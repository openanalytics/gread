#' @title Quick and easy reading of gtf/gff/bed/bam files
#'
#' @description \code{read_gtf}, \code{read_gff}, \code{read_bed} and 
#' \code{read_bam} are specific functions corresponding to respective file 
#' formats. 
#' 
#' \code{read_format} is a top-level convenience function that detects and 
#' reads automatically from the extension and is sufficient in most cases.
#'
#' @param file Complete path to the input file.
#' @param format Type of annotation. Currently allowed values are 
#' \code{"gtf"}, \code{"gff"}, \code{"bed"} and \code{"bam"}. If no input is 
#' provided, then it will be guessed from input file's extension.
#' @param filter Filter rows where \code{start/stop} are invalid, and where 
#' \code{start > end}. Default (\code{FALSE}) is to not filter, i.e., retain 
#' all invalid rows as well.
#' @param chromosomes Argument to \code{read_bam}. If specified only reads 
#' from those chromosomes will be loaded.
#' @param tags Additional (optional) tags to load from the bam file. By 
#' default, the tags \code{"NM"} and \code{"MD"} are loaded.
#' @param verbose If \code{TRUE}, sends useful status messages to the console. 
#' Default is \code{FALSE}.
#' @aliases read_format read_gtf read_gff read_bed read_bam
#' @return An object of class \code{gtf}, \code{gff}, \code{bed} or \code{bam}, 
#' corresponding to the input file format, that inherits from \code{GRanges}.
#' @seealso \code{\link{supported_formats}} 
#' \code{\link{extract}} \code{\link{construct_introns}}
#' @keywords file
#' @export
#' @examples
#' path <- system.file("tests", package="gread")
#' gff_file <- file.path(path, "sample.gff")
#' gtf_file <- file.path(path, "sample.gtf")
#' bed_file <- file.path(path, "sample.bed")
#' bam_file <- file.path(path, "sample.bam")
#' 
#' read_format(gff_file) # read gff
#' read_gff(gff_file)    # same as above
#' read_format(gtf_file) # read gtf
#' read_gtf(gtf_file)    # same as above
#' read_format(bed_file) # read bed
#' read_bed(bed_file)    # same as above
#' read_format(bam_file) # read bam
#' read_bam(bam_file)    # same as above
#' 
#' gtf_filter_file <- file.path(path, "sample_filter.gtf")
#' read_format(gtf_filter_file, filter=TRUE) # filter invalid rows
#' read_gtf(gtf_filter_file, filter=TRUE)    # same as above
read_format <- function(file, format=detect_format(file), filter=FALSE, 
            chromosomes=NULL, tags=c("NM", "MD"), verbose=FALSE) {
    
    if (!file.exists(file))
        stop(toupper(format), " file '", file, "' not found.")
    types = format_types(format=format)
    names = format_names(format=format)
    token = tokenise(file, format, filter, chromosomes, tags, verbose, 
                        types=types, names=names)
    ans = gread(token)
    tidy_cols(ans, verbose=verbose)
    if (token$filter && !identical(format, "bam")) {
        before = nrow(ans)
        # TODO: optimise this subset using internal any in data.table
        ans = ans[start <= end] # will take care of NAs as well
        after = nrow(ans)
        if (token$verbose && before-after) 
            cat(before-after, " invalid rows filtered. These include rows",
                " where start/end coordinates were NA/NaN and ", 
                " where start > end.\n", sep="")
    }
    # Returning 'GRanges' following Hervé's feedback
    # See https://github.com/Bioconductor/Contributions/issues/25
    new(class(ans)[1L], as(setDF(ans), "GRanges"))
}

#' @rdname read_format
#' @export
read_gtf <- function(file, filter=FALSE, verbose=FALSE) {
    read_format(file, "gtf", filter, verbose)
}

#' @rdname read_format
#' @export
read_gff <- function(file, filter=FALSE, verbose=FALSE) {
    read_format(file, "gff", filter, verbose)
}

#' @rdname read_format
#' @export
read_bed <- function(file, filter=FALSE, verbose=FALSE) {
    read_format(file, "bed", filter, verbose)
}

#' @rdname read_format
#' @export
read_bam <- function(file, filter=FALSE, chromosomes=NULL, 
                        tags=c("NM", "MD"), verbose=FALSE) {
    read_format(file, "bam", filter, chromosomes, tags, verbose)
}

# Internal helper functions --------------------------------------------------

detect_format <- function(file) {

    type = substr(tolower(tools::file_ext(file)), 1L, 3L)
    if (!type %in% c("gtf", "gff", "bam", "bed"))
        stop("File does not have gtf/gff/bed/bam extension.")
    type
}

format_types <- function(format) {
    chr = "character"
    int = "integer"
    num = "numeric"
    switch(format, gtf=, gff=c(rep(chr, 3L), rep(int, 2L), num, chr, int, chr),
                   bed=c(chr, rep(int, 2L), rep(chr, 3L), rep(int, 2L), 
                            chr, int, rep(chr, 2L)),
                   bam=c())
}

format_names <- function(format) {
    switch(format, gtf=c("seqnames", "source", "feature", "start", 
                    "end", "score", "strand", "frame", "attributes"), 
                   gff=c("seqnames", "source", "feature", "start", 
                    "end", "score", "strand", "phase", "attributes"), 
                   # all 15 columns for bed, only first 3 are compulsory
                   bed=c("seqnames", "start", "end", "name", "score", 
                    "strand", "thickStart", "thickEnd", "itemRgb", 
                    "blockCount", "blockSizes", "blockStarts", "expCount", 
                    "expIds", "expScores"), 
                   bam=c("seqnames", "start", "end", "strand", "cigar", 
                    "qwidth", "width", "njunc", "flag"))
}

tokenise <- function(file, format, filter, chromosomes, 
                       tags, verbose, ...) {
    arglist = as.list(match.call(expand.dots=TRUE))[-1L]
    argnames = setdiff(names(arglist), "format")
    if (!length(tags)) tags = character(0) # ScanBamParam needs this
    tokens = structure(c(list(file, filter, chromosomes, tags, verbose), 
                list(...)), class=paste(format, "format", sep="_"))
    data.table::setattr(tokens, 'names', argnames)
    tokens
}

`.` <- function(...) {
    stop(".() is only meant for use in `j` argument of 
        data.table, and is an alias to list()")
}
