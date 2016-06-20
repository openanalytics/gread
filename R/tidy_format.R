#' @title Tidy up GTF/GFF/BED/BAM objects
#'
#' @description Tidy up data depending on the type of input format. 
#'
#' Note that this operation is performed \emph{by reference}, to avoid 
#' unnecessary copies, for efficiency. Therefore there is no need to assign 
#' the result back to a variable. 
#'
#' @details In case of \code{GTF} and \code{GFF} formats, the 
#' \code{attributes} column is extracted into separate columns after which 
#' it is removed (by default). Set \code{remove_cols} to \code{NULL} if you 
#' would like to retain the column even after tidying up.
#' 
#' In case of \code{BED} and \code{BAM} files, no columns are removed by 
#' default. However, you can use \code{remove_cols} argument if necessary.
#' 
#' In case of \code{BAM} files, the function \code{GAlignments} from the 
#' \code{GenomicAlignments} Bioconductor package is used to read in the 
#' file, with the additional arguments for \code{ScanBamParam()} function 
#' using \code{NM} and \code{MD} tags by default.
#'
#' @seealso \code{\link{supported_formats}} \code{\link{read_format}} 
#' \code{\link{extract}} \code{\link{construct_introns}} 
#' \code{\link{as_granges}}
#' @param x Input object, of class \code{GTF}, \code{GFF}, \code{BED} or 
#' \code{BAM}.
#' @param remove_cols A character vector of column names to be removed.
#' @param verbose Logical. Default is \code{FALSE}. If \code{TRUE}, helpful 
#' status messages are printed on to the console. 
#' @param ... Arguments that are ignored at the moment.
#' @aliases tidy tidy.gtf tidy.gff tidy.bed tidy.bam
#' @return A tidied object of class \code{GTF}, \code{GFF}, \code{BED} or 
#' \code{BAM}, that inherits from \code{data.table}.
#' @export
#' @examples
#' path = system.file("tests", package="gread")
#' gff_file = file.path(path, "sample.gff")
#' gtf_file = file.path(path, "sample.gtf")
#' bed_file = file.path(path, "sample.bed")
#' bam_file = file.path(path, "sample.bam")
#' 
#' gtf = read_format(gtf_file, tidy=FALSE)
#' gff = read_format(gtf_file, tidy=FALSE)
#' bed = read_format(bed_file, tidy=FALSE)
#' bam = read_format(bam_file, tidy=FALSE)
#' 
#' tidy(gtf, remove_cols=NULL)[]          # tidy attributes col, but not remove
#' tidy(gff, remove_cols=NULL)[]          # same as above, but for gff
#' tidy(gtf, remove_cols="attributes")[]  # tidy *and* remove attributes col
#' tidy(gff, remove_cols="attributes")[]  # same as above, but for gff
#' 
#' tidy(bed, remove_cols="name")[]        # remove name column
#' tidy(bam, remove_cols=c("NM", "MD"))[] # remove additional loaded tags 
tidy <- function(x, ...) {
    UseMethod("tidy")
}

#' @rdname tidy
#' @export
tidy.default <- function(x, ...) {
    stop("No default method available.")
}

#' @rdname tidy
#' @export
tidy.gtf <- function(x, remove_cols="attributes", verbose=FALSE, ...) {
    stopifnot(identical(head(names(x), 9L), format_names("gtf")))
    meta_fun <- function(vec, attr) {
        data.table::setattr(as.list(vec), 'names', attr)
    }
    if (ncol(x) > 9L) { # rtracklayer::readGFF output, already tidied
        for (col in which(vapply(x, is.character, TRUE)))
            set(x, i=which(x[[col]] == ""), j=col, value=NA)
    } else if ("attributes" %in% names(x)) {
        if (verbose) cat("Tidying up 'attributes' column.\n")
        # base R's strsplit with perl=TRUE is faster than without, 
        # but for 'gtf', is twice as slow as stringi (20 vs 10 s)
        prep = stringi::stri_split_regex(x$attributes, 
                  "[ ]\"|\";[ ]*", omit_empty=TRUE)
        cols = lapply(prep, `[`, c(TRUE, FALSE))
        vals = lapply(prep, `[`, c(FALSE, TRUE))
        len1 = vapply(cols, length, 0L)
        len2 = vapply(vals, length, 0L)
        if (!identical(len1, len2)) {
            stop("Ill-formed gff file. Attributes column in gff file must be 
                of the form 'id1=val1;id2=val2;...")
        }
        prep = setDT(list(id=rep.int(seq_along(len1), len1), 
                          cols=unlist(cols), vals=unlist(vals)))
        meta = dcast(prep, id ~ cols, value.var="vals")
        set(meta, j="id", value=NULL)
        # x = data.table:::shallow(x) # TODO: use when exported
        x[, names(meta) := meta]
    }
    remove_cols(x, remove_cols, verbose) # updates by reference
    if (!all(c("transcript_id", "gene_id") %in% names(x)))
        warning("Columns 'transcript_id' and/or 'gene_id' are not found in ", 
            "the GTF file or has been removed while tidying. Note that ", 
            "these columns are essential for downstream analysis in ", 
            "most cases.")
    invisible(x)
}

#' @rdname tidy
#' @export
tidy.gff <- function(x, remove_cols="attributes", verbose=FALSE, ...) {
    stopifnot(identical(head(names(x), 9L), format_names("gff")))
    meta_fun <- function(vec, attr) {
        data.table::setattr(as.list(vec), 'names', attr)
    }
    if (ncol(x) > 9L) { # rtracklayer::readGFF output, already tidied
        for (col in which(vapply(x, is.character, TRUE)))
            set(x, i=which(x[[col]] == ""), j=col, value=NA)
    } else if ("attributes" %in% names(x)) {
        if (verbose) cat("Tidying up 'attributes' column.\n")
        # base R's strsplit with perl=TRUE is faster than without, 
        # but for 'gtf', is twice as slow as stringi (20 vs 10 s)
        prep = stringi::stri_split_regex(x$attributes, "=|;", omit_empty=TRUE)
        cols = lapply(prep, `[`, c(TRUE, FALSE))
        vals = lapply(prep, `[`, c(FALSE, TRUE))
        len1 = vapply(cols, length, 0L)
        len2 = vapply(vals, length, 0L)
        if (!identical(len1, len2)) {
            stop("Ill-formed gff file. Attributes column in gff file must be 
                    of the form 'id1=val1;id2=val2;...")
        }
        prep = setDT(list(id=rep.int(seq_along(len1), len1), 
                          cols=unlist(cols), vals=unlist(vals)))
        meta = dcast(prep, id ~ cols, value.var="vals")
        set(meta, j="id", value=NULL)
        # x = data.table:::shallow(x) # TODO: use when exported
        x[, names(meta) := meta]
    }
    gff_gene_transcript_cols(x) # add transcript_id, gene_id cols
    # and if possible transcript_name and gene_name cols as well
    remove_cols(x, remove_cols, verbose) # updates by reference
    invisible(x)
}

#' @rdname tidy
#' @export
tidy.bed <- function(x, remove_cols=NULL, verbose=FALSE, ...) {
    stopifnot(identical(names(x), head(format_names("bed"), length(x))))
    remove_cols(x, remove_cols, verbose) # updates by reference
    invisible(x)
}

#' @rdname tidy
#' @export
tidy.bam <- function(x, remove_cols=NULL, verbose=FALSE, ...) {
    remove_cols(x, remove_cols, verbose) # updates by reference
    invisible(x)
}

## internal function used in tidy methods -------------------

remove_cols <- function(x, cols, verbose) {
    if (is.null(cols)) {
        if (verbose) cat("remove_cols is NULL. No columns are removed.")
    } else {
        invalid = cols[!cols %in% names(x)]
        if (length(invalid)) {
            warning("Skipping column(s) [", paste(invalid, 
                collapse=", "), "] as they are not present in 'x'.")
            cols = setdiff(cols, invalid)
        }
        if (length(cols)) {
            if (verbose) cat("Removing [", paste(cols, collapse=", "), 
                "] columns.\n", sep="")
            set(x, j=cols, value=NULL)
        }
    }
    x[]
}

gff_gene_transcript_cols <- function(x) {
    # to please R CMD CHECK
    gene_id=transcript_id=NULL
    cols = c("ID", "Parent")
    if (any(check <- !cols %in% names(x))) {
        warning("GFF file doesn't contain ", paste(cols[check], collapse=","), 
        "column(s) in its attributes column. They are necessary for ", 
        " creating 'transcript_id' and 'gene_id' columns which will be", 
        " required for most downstream analysis. Please fix the GFF file.")
        cols = cols[!check]
        x[, c(cols) := lapply(.SD, function(x) gsub("^.*:", "", x)), 
                    .SDcols=cols]
    } else {
        ID=Parent=i.Parent=Name=i.Name=NULL # suppress R CMD CHECK warnings
        x[, "ID" := gsub("^.*:", "", ID)]
        x[, "Parent" := gsub("^.*:", "", Parent)]
        # transcript_id
        rnaids = grepl("rna", x$feature, ignore.case=TRUE)
        x[, "transcript_id" := Parent][rnaids, "transcript_id" := ID]
        # gene_id
        tmp = x[rnaids, list(transcript_id=ID, Parent)]
        if (anyDuplicated(tmp)) tmp = unique(tmp, by=names(tmp))
        x[, gene_id := ID][tmp, gene_id := i.Parent, on="transcript_id"]
        if ("Name" %in% names(x)) {
            # transcript_name
            tmp = x[rnaids, list(transcript_id, Name)]
            if (anyDuplicated(tmp)) tmp = unique(tmp, by=names(tmp))
            i.Name = NULL
            x[tmp, "transcript_name" := i.Name, on="transcript_id"]
            # gene_name
            geneids = grepl("gene", x$feature, ignore.case=TRUE)
            tmp = x[geneids, list(gene_id, Name)]
            if (anyDuplicated(tmp)) tmp = unique(tmp, by=names(tmp))
            x[tmp, "gene_name" := i.Name, on="gene_id"]
        }
    }
}
