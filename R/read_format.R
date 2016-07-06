#' @title Quickly and easily read GTF/GFF/BED/BAM files
#'
#' @description \code{read_gtf}, \code{read_gff}, \code{read_bed} and 
#' \code{read_bam} are specific functions corresponding to respective file 
#' formats. 
#' 
#' \code{read_format} is a top-level convenience function that detects and 
#' reads automatically from the extension and is sufficient in most cases.
#'
#' @details Note that \code{GFF1} uses \code{group} instead of 
#' \code{attributes} as the column name, but \code{gread} always names it 
#' as \code{attributes}. Similarly the first three columns of \code{bed} 
#' format are named \code{seqname, start, end} instead of \code{chrom, 
#' chromStart, chromEnd} for consistency.
#' 
#' The argument \code{tidy} (\code{TRUE} by default) tidies up the 
#' \code{attributes} column in case of \code{GTF/GFF} format files. The 
#' \code{attributes} column itself is removed since it is tidied up into 
#' multiple columns. If this is not desirable, use \code{tidy = FALSE} to 
#' load the file \emph{as-is} and then use the \code{tidy} function with 
#' \code{remove_cols = NULL}.
#' 
#' In case of \code{GFF} format, when `tidy=TRUE`, generation of columns 
#' \code{"transcript_id"} and \code{"gene_id"} will be attempted as these 
#' columns are essential in most cases for downstream analyses. If possible, 
#' the columns \code{"transcript_name"} and \code{"gene_name"} will also be 
#' created.
#' 
#' @param file Complete path to the input file.
#' @param format Type of annotation. Currently allowed values are 
#' \code{"GTF"}, \code{"GFF"}, \code{"BED"} and \code{"BAM"}. If no input is 
#' provided, then it will be guessed from input file's extension.
#' @param filter Filter rows where \code{start/stop} are invalid, and where 
#' \code{start > end}. Default is \code{FALSE} (not to filter, but retain all 
#' invalid rows).
#' @param tidy If \code{TRUE} (default), returns only essential columns for 
#' further analysis (for e.g, \code{score}, \code{itemRgb} etc. are removed), 
#' \code{attributes} column is cleaned up with separate columns for 
#' \code{gene}, \code{transcript} id etc. Default is \code{TRUE}.
#' @param chromosomes Argument to \code{read_bam}. If specified only reads 
#' from those chromosomes will be loaded.
#' @param tags Additional (optional) tags to load from the bam file. By 
#' default, the tags \code{"NM"} and \code{"MD"} are loaded.
#' @param verbose If \code{TRUE}, sends useful status messages to the console. 
#' Default is \code{FALSE}.
#' @aliases read_format read_gtf read_gff read_bed read_bam
#' @return An object of class \code{gtf}, \code{gff}, \code{bed} or \code{bam}, 
#' corresponding to the input file format, that inherits from \code{data.table}.
#' @seealso \code{\link{supported_formats}} \code{\link{tidy}} 
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
#' read_format(gff_file) # read GFF
#' read_gff(gff_file)    # same as above
#' read_format(gtf_file) # read GTF
#' read_gtf(gtf_file)    # same as above
#' read_format(bed_file) # read BED
#' read_bed(bed_file)    # same as above
#' read_format(bam_file) # read BAM
#' read_bam(bam_file)    # same as above
#' 
#' read_format(gtf_file, tidy=FALSE) # load as is, don't tidy
#' 
#' gtf_filter_file <- file.path(path, "sample_filter.gtf")
#' read_format(gtf_filter_file, filter=TRUE) # filter invalid rows
#' read_gtf(gtf_filter_file, filter=TRUE)    # same as above
read_format <- function(file, format=detect_format(file), filter=FALSE, 
            tidy=TRUE, chromosomes=NULL, tags=c("NM", "MD"),  verbose=FALSE) {
    
    if (!file.exists(file))
        stop(toupper(format), " file '", file, "' not found.")
    types = format_types(format=format)
    names = format_names(format=format)
    token = tokenise(file, format, filter, tidy, chromosomes, tags, verbose, 
                        types=types, names=names)
    ans = gread(token)
    if (tidy) tidy(ans, verbose=verbose)
    if (token$filter && !identical(format, "bam")) {
        before = nrow(ans)
        # TO DO: optimise this subset using internal any in data.table
        ans = ans[start <= end] # will take care of NAs as well
        after = nrow(ans)
        if (token$verbose && before-after) 
            cat(before-after, " invalid rows filtered. These include rows",
                " where start/end coordinates were NA/NaN and ", 
                " where start > end.\n", sep="")
    }
    ans
}

#' @rdname read_format
#' @export
read_gtf <- function(file, filter=FALSE, tidy=TRUE, verbose=FALSE) {
    read_format(file, "gtf", filter, tidy, verbose)
}

#' @rdname read_format
#' @export
read_gff <- function(file, filter=FALSE, tidy=TRUE, verbose=FALSE) {
    read_format(file, "gff", filter, tidy, verbose)
}

#' @rdname read_format
#' @export
read_bed <- function(file, filter=FALSE, tidy=TRUE, verbose=FALSE) {
    read_format(file, "bed", filter, tidy, verbose)
}

#' @rdname read_format
#' @export
read_bam <- function(file, filter=FALSE, tidy=TRUE, chromosomes=NULL, 
                        tags=c("NM", "MD"), verbose=FALSE) {
    read_format(file, "bam", filter, tidy, chromosomes, tags, verbose)
}

# Helper/Internal functions for read_format ---------------------

detect_format <- function(file) {

    type = substr(tolower(tools::file_ext(file)), 1L, 3L)
    if (!type %in% c("gtf", "gff", "bam", "bed"))
        stop("File does not have GTF/GFF/BED/BAM extension.")
    type
}

tokenise <- function(file, format, filter, tidy, chromosomes, 
                       tags, verbose, ...) {
    arglist = as.list(match.call(expand.dots=TRUE))[-1L]
    argnames = setdiff(names(arglist), "format")
    if (!length(tags)) tags = character(0) # ScanBamParam needs this
    tokens = structure(c(list(file, filter, tidy, chromosomes, tags, verbose), 
                list(...)), class=paste(format, "format", sep="_"))
    data.table::setattr(tokens, 'names', argnames)
    tokens
}

# Helper/Internal functions for read_gtf/read_gff/read_bed -------

format_names <- function(format) {
    switch(format, gtf=c("seqname", "source", "feature", "start", 
                    "end", "score", "strand", "frame", "attributes"), 
                   gff=c("seqname", "source", "feature", "start", 
                    "end", "score", "strand", "phase", "attributes"), 
                   # all 15 columns for bed, only first 3 are compulsory
                   bed=c("seqname", "start", "end", "name", "score", 
                    "strand", "thickStart", "thickEnd", "itemRgb", 
                    "blockCount", "blockSizes", "blockStarts", "expCount", 
                    "expIds", "expScores"), 
                   bam=c("seqname", "start", "end", "strand", "cigar", 
                    "qwidth", "width", "njunc", "flag"))
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

read_table <- function(file, ...) 
    setDT(read.table(file, stringsAsFactors=FALSE, as.is=TRUE, ...))

gread <- function(token) {
    UseMethod("gread")
}

gread.default <- function(token) {
    stop("No method available for format ", gsub("(.*)_format$", 
            "", class(token)))
}

gread.gtf_format <- function(token) {
    # to please R CMD CHECK
    score=phase=NULL
    # rtracklayer::readGFF reads GFF v1 'attributes' column as 'group',
    # but we want it to be always 'attributes'. Also, it tidies up the 
    # attributes column by default
    rtracklayer_fun <- function() {
        ans = setDT(rtracklayer::readGFF(token$file, version = 2L, 
                columns = rtracklayer::GFFcolnames()))
        if (!token$tidy) set(ans, j = tail(names(ans), -9L), value = NULL)
        ans[, (1:3) := lapply(.SD, as.character), .SDcols=1:3][]
        # remove additional attributes
        setattr(ans, 'pragmas', NULL)
        setattr(ans, 'attrcol_fmt', NULL)
        setattr(ans, 'ncol0', NULL)
        setattr(ans, 'ntag', NULL)
        setattr(ans, 'raw_data', NULL)
        ans
    }
    fread_fun <- function() fread(token$file, colClasses = token$types, 
        showProgress = token$verbose)
    read_table_fun <- function() read_table(token$file, sep = "\t", 
        header = FALSE, comment.char = "#", nrows = -1L, colClasses = 
        token$types, quote = "")
    ans = tryCatch(rtracklayer_fun(), error = function(o) { 
                if (token$verbose) cat("rtracklayer::readGFF failed to read", 
                    " the GTF file. Reverting to data.table::fread.\n", sep="")
                tryCatch(fread_fun(), error = function(o) {
                    if (token$verbose) cat("data.table::fread failed to read", 
                        " as well. Reverting to base::read.table.\n", sep="")
                    read_table_fun()
                })
          })
    setnames(ans, head(names(ans), 9L), token$names)
    if (is.character(ans[["score"]])) 
        ans[, "score" := suppressWarnings(as.numeric(score))]
    if (is.character(ans[["phase"]])) 
        ans[, "phase" := suppressWarnings(as.integer(phase))]
    setattr(ans, 'class', c("gtf", "data.table", "data.frame"))
}

gread.gff_format <- function(token) {
    # to please R CMD CHECK
    score=phase=NULL
    # rtracklayer::readGFF reads GFF v1 'attributes' column as 'group',
    # but we want it to be always 'attributes'. Also, it tidies up the 
    # attributes column by default
    rtracklayer_fun <- function() {
        ans = setDT(as.data.frame(rtracklayer::readGFF(token$file, columns = 
                rtracklayer::GFFcolnames())))
        if (!token$tidy) set(ans, j = tail(names(ans), -9L), value = NULL)
        ans[, (1:3) := lapply(.SD, as.character), .SDcols=1:3][]
        list_cols = which(vapply(ans, is.list, TRUE))
        if (length(list_cols)) {
            list_paste <- function(x) sapply(x, paste, collapse=",")
            # paste all values of list cols together
            ans[, (list_cols) := lapply(.SD, list_paste), .SDcols=list_cols]
        }
        # as.data.frame drops all extra attributes
        ans
    }
    fread_fun <- function() fread(token$file, colClasses = token$types, 
        showProgress = token$verbose)
    read_table_fun <- function() read_table(token$file, sep = "\t", 
        header = FALSE, comment.char = "#", nrows = -1L, colClasses = 
        token$types, quote = "")
    ans = tryCatch(rtracklayer_fun(), error = function(o) { 
                if (token$verbose) cat("rtracklayer::readGFF failed to read",
                    " the GTF file. Reverting to data.table::fread.\n", sep="")
                tryCatch(fread_fun(), error = function(o) {
                    if (token$verbose) cat("data.table::fread failed to read", 
                        " as well. Reverting to base::read.table.\n", sep="")
                    read_table_fun()
                })
          })
    setnames(ans, head(names(ans), 9L), token$names)
    if (is.character(ans[["score"]])) 
        ans[, "score" := suppressWarnings(as.numeric(score))]
    if (is.character(ans[["phase"]])) 
        ans[, "phase" := suppressWarnings(as.integer(phase))]
    setattr(ans, 'class', c("gff", "data.table", "data.frame"))
}

gread.bed_format <- function(token) {
    ans = tryCatch(fread(token$file, colClasses = token$types, sep="\t", 
                showProgress = token$verbose), error = function(o) {
                    if (token$verbose) cat("data.table::fread failed to read", 
                        " the bed file. Reverting to 'read.table'", 
                        " from base R.\n", sep="")
                    read_table(token$file, sep = "\t", header = FALSE, 
                        comment.char = "#", nrows = -1L, colClasses = 
                        token$types, quote = "")
                    })
    # setnames will error if bed file has < 3 cols
    setnames(ans, head(token$names, max(3L, ncol(ans))))
    # ans[, start := start+1L] # TODO: should we automatically set start+1L?
    setattr(ans, 'class', c("bed", "data.table", "data.frame"))
}

gread.bam_format <- function(token) {
    if (token$verbose) 
        cat("Loading bam file: ",token$file,", please wait...\n", sep="")
    stopifnot(length(token$file)==1L)
    
    what  = c("flag")
    flags = Rsamtools::scanBamFlag(isUnmappedQuery=FALSE)
    if (!length(token$chromosomes)) {
        param = Rsamtools::ScanBamParam(what=what, tag=token$tags, flag=flags)
    } else {
            hdr = Rsamtools::scanBamHeader(token$file)[[1]]$targets
            hdr = hdr[names(hdr) %in% token$chromosomes]
            if (!all(names(hdr) %in% token$chromosomes)) {
                stop("chromosomes [", paste(setdiff(token$chromosomes, 
                    names(hdr)), collapse=", "), "] not found in bam file.")
            }
            which = as_granges(data.table(seqname = names(hdr), 
                        start=1L, end = hdr, strand="*"), ignore_strand=TRUE)
            param = Rsamtools::ScanBamParam(what = what, which = which, 
                                tag = token$tags, flag = flags)
    }
    # load bam with corresponding param
    ans = as_bam(GenomicAlignments::readGAlignments(token$file, param=param))
    setattr(ans, 'class', c("bam", "data.table", "data.frame"))
}

`.` <- function(...) {
    stop(".() is only meant for use in `j` argument of 
        data.table, and is an alias to list()")
}
