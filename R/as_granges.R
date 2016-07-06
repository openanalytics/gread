## internal functions ---------------------------------

#' @title Convert to GRanges object
#' 
#' @description Convert \code{gtf}, \code{gff}, \code{bed}, \code{bam} 
#' or a valid \code{data.table} to a GRanges object, preserving all 
#' additional columns.
#'
#' @param x An object of class \code{gtf}, \code{gff}, \code{bed} or 
#' \code{bam} or a valid \code{data.table} object.
#' @param ignore_strand Logical argument to pass to \code{GRanges} function. 
#' Indicates whether \code{strand} should be ignored when constructing 
#' \code{GRanges} object or not. Default is \code{FALSE}.
#' @return A \code{GRanges} object.
#' @export
#' @aliases as_granges
#' @examples
#' path <- system.file("tests", package="gread")
#' gff_file <- file.path(path, "sample.gff")
#' gtf_file <- file.path(path, "sample.gtf")
#' bed_file <- file.path(path, "sample.bed")
#' bam_file <- file.path(path, "sample.bam")
#' 
#' gff <- read_format(gff_file)
#' gtf <- read_format(gtf_file)
#' bed <- read_format(bed_file)
#' bam <- read_format(bam_file)
#' 
#' as_granges(gff)
#' as_granges(gtf)
#' as_granges(bed)
#' as_granges(bam)
#' 
#' as_granges(gff, ignore_strand=FALSE)
#' as_granges(gtf, ignore_strand=FALSE)
#' as_granges(bed, ignore_strand=FALSE)
#' as_granges(bam, ignore_strand=FALSE)
#' @seealso \code{\link{read_format}} 
#' \code{\link{tidy}} \code{\link{extract}} \code{\link{construct_introns}}
as_granges <- function(x, ignore_strand=FALSE) {

    stopifnot(is.gtf(x)||is.gff(x)||is.bed(x)||is.bam(x)||is.data.table(x))
    x = shallow(x)
    if (ignore_strand) x[, "strand" := NULL]
    as(setDF(x), "GRanges")
}

#' @title Return only those rows where rows per group is > 1.
#'
#' @description For internal use only. It returns only rows corresponding to 
#' those groups where the number of rows > 1L. This is equivalent to doing 
#' \code{DT[, if (.N>1L) .SD, by=c(...)]}.
#' @param x A \code{data.table}.
#' @param cols Character vector of column names (from \code{x}).
#' @return A \code{data.table}.
#' @examples
#' \dontrun{
#' require(data.table)
#' dt <- data.table(x=c(1,1,1,2,2), y=c(3,3,4,5,6))
#' strictly_nonunique(dt) # Only 1,3 occurs more than once
#' strictly_nonunique(dt, "x") # all values occur more than once
#' }
strictly_nonunique <- function(x, cols=names(x)) {
    # to please R CMD CHECK
    `..N..`=NULL
    stopifnot(is.data.table(x), all(cols %in% names(x)))
    shallow(x)[, "..N.." := .N, by=c(cols)
             ][`..N..` > 1L
             ][, "..N.." := NULL][] 
             # last [] is to ensure printing on first time
}

#' @title Find overlapping indices of two gtf/gff/bed/bam objects
#'
#' @description For internal use only. Function for finding overlaps between 
#' two objects of class \code{gtf/gff/bed/bam} using 
#' \code{GenomicRanges::findOverlaps}.
#' 
#' @param x,y An object of class \code{gtf}, \code{gff}, \code{bed} or 
#' \code{bam}.
#' @param ignore_redundant Should redundant overlaps be ignored?
#' @param ignore_strand Logical argument to pass to \code{GRanges} function. 
#' Indicates whether \code{strand} should be ignored when constructing 
#' \code{GRanges} object or not. Default is \code{FALSE}.
#' @param ... Additional arguments passed to 
#' \code{GenomicRanges::findOverlaps}.
#' @return A \code{data.table} containing overlapping indices.
find_overlaps <- function(x, y, ignore_redundant=FALSE, 
                    ignore_strand=FALSE, ...) {
    stopifnot(is.gtf(x) || is.gff(x) || is.bam(x) || is.bed(x), 
                is.gtf(y) || is.gff(y) || is.bam(y) || is.bed(y), 
                ignore_redundant %in% c(FALSE, TRUE), 
                ignore_strand %in% c(FALSE, TRUE))
    x = as_granges(x, ignore_strand)
    y = as_granges(y, ignore_strand)
    olaps = GenomicRanges::findOverlaps(x, y, ...)
    olaps = setDT(list(queryHits = queryHits(olaps), 
                    subjectHits = subjectHits(olaps)))
    # findOverlaps for GRanges objects doesn't seem to have ignoreRedundant
    # argument. so mimicing that functionality below.
    if (ignore_redundant) {
        olaps = olaps[, `:=`(queryHits = pmin(queryHits, subjectHits), 
                            subjectHits = pmax(queryHits, subjectHits))]
        olaps = unique(olaps, by=names(olaps))
    }
    olaps[]
}

#' @title Compute reduced ranges on a gtf/gff/bed/bam object.
#' 
#' @description For internal use only. Computes reduced ranges on a 
#' \code{gtf/gff/bed/bam} object by converting the input object to a 
#' \code{GRangesList} object and calling \code{reduce()} (from 
#' \code{GenomicRanges} package) on it. Returns an object of same class 
#' as input.
#'
#' @param x An object of class \code{gtf/gff/bed/bam}.
#' @param by Character vector of column names in \code{x} to \emph{group by}.
#' @param ignore_strand Logical argument to pass to \code{GRanges} function. 
#' Indicates whether \code{strand} should be ignored when constructing 
#' \code{GRanges} object or not. Default is \code{FALSE}.
#' @return A \code{data.table} with reduced ranges.
reduce_overlaps <- function(x, by="gene_id", ignore_strand=FALSE) {
    # to please R CMD CHECK
    group=seqname=i.strand=`..N..`=NULL
    stopifnot(length(by) == 1L, by %in% names(x), is.gtf(x) || is.gff(x) 
        || is.bed(x) || is.bam(x), ignore_strand %in% c(FALSE, TRUE))
    red = as.data.frame(GenomicRanges::reduce(GenomicRanges::split(as_granges(
                            x, ignore_strand), x[[by]])))
    setDT(red)[, group := NULL]
    setcolorder(red, c("seqnames", "start", "end", "width", "strand", 
                        "group_name"))
    setnames(red, c("seqnames", "width", "group_name"), c("seqname", 
                        "length", by))
    red[, `:=`(seqname=as.character(seqname), strand=as.character(strand))]
    # restore original order, nomatch = "errror" would be great to have here!
    ux = unique(shallow(x, c(by, "strand")), by=c(by))
    red = red[ux, on=c(by)]
    # if ignore_strand, replace strand
    if (ignore_strand) red[, strand := i.strand]
    red[, i.strand := NULL]
    setattr(red, 'class', class(x))
    red[]
}

#' @title Compute disjoint ranges on a gtf/gff/bed/bam object.
#' 
#' @description For internal use only. Computes disjoint ranges on a 
#' \code{gtf/gff/bed/bam} object by converting the input object to a 
#' GRangesList object and calling \code{disjoin()} (from the 
#' \code{GenomicRanges} package) on it. Returns an object of same class 
#' as input.
#' 
#' @param x An object of class \code{gtf/gff/bed/bam}.
#' @param by Character vector of column names in \code{x} to \emph{group by}.
#' @param ignore_strand Logical argument to pass to \code{GRanges} function. 
#' Indicates whether \code{strand} should be ignored when constructing 
#' \code{GRanges} object or not. Default is \code{FALSE}.
#' @return A \code{data.table} with disjoint ranges.
disjoin_overlaps <- function(x, by="gene_id", ignore_strand=FALSE) {
    stopifnot(length(by) == 1L, by %in% names(x), is.gtf(x) || is.gff(x) || 
        is.bed(x) || is.bam(x), ignore_strand %in% c(FALSE, TRUE))

    # to please R CMD CHECK
    group=seqname=i.strand=NULL
    dj = as.data.frame(GenomicRanges::disjoin(GenomicRanges::split(as_granges(
                        x, ignore_strand), x[[by]])))
    setDT(dj)[, group := NULL]
    setcolorder(dj, c("seqnames", "start", "end", "width", "strand", 
                        "group_name"))
    setnames(dj, c("seqnames", "width", "group_name"), c("seqname", 
                        "length", by))
    dj[, `:=`(seqname=as.character(seqname), strand=as.character(strand))]
    # restore original order
    # TODO: revisit when nomatch = "error" is implemented
    ux = unique(shallow(x, c(by, "strand")), by=c(by))
    dj = dj[ux, on=c(by)]
    # if ignore_strand, replace strand
    if (ignore_strand) dj[, strand := i.strand]
    dj[, i.strand := NULL]
    setattr(dj, 'class', class(x))
    dj[]
}

#' @title Compute intersecting ranges on a gtf/gff/bed/bam object with itself.
#' 
#' @description For internal use only. Computes intersecting ranges on a 
#' \code{gtf/gff} object by converting the input object to a GRangesList 
#' object and calling \code{disjoin()} (from the \code{GenomicRanges} 
#' package) on it. Returns an object of same class as input.
#' 
#' @param x An object of class \code{gtf} or \code{gff}.
#' @param by Character vector of column names in \code{x} to \emph{group by}.
#' @param ignore_strand Logical argument to pass to \code{GRanges} function. 
#' Indicates whether \code{strand} should be ignored when constructing 
#' \code{GRanges} object or not. Default is \code{FALSE}.
#' @return A \code{data.table} with intersecting ranges.
intersect_overlaps <- function(x, by="gene_id", ignore_strand=FALSE) {
    # to please R CMD CHECK
    group=seqname=i.strand=NULL
    stopifnot(length(by) == 1L, by %in% names(x), is.gtf(x) || is.gff(x) || 
        is.bed(x) || is.bam(x), ignore_strand %in% c(FALSE, TRUE))
    s_gr = GenomicRanges::split(as_granges(x, ignore_strand), x[[by]])
    isect = as.data.frame(GenomicRanges::intersect(s_gr, s_gr))
    setDT(isect)[, group := NULL]
    setcolorder(isect, c("seqnames", "start", "end", "width", "strand", 
                        "group_name"))
    setnames(isect, c("seqnames", "width", "group_name"), c("seqname", 
                        "length", by))
    isect[, `:=`(seqname=as.character(seqname), strand=as.character(strand))]
    # restore original order, nomatch = "errror" would be great to have here!
    ux = unique(shallow(x, c(by, "strand")), by=c(by))
    isect = isect[ux, on=c(by)]
    # if ignore_strand, replace strand
    if (ignore_strand) isect[, strand := i.strand]
    isect[, i.strand := NULL]
    setattr(isect, 'class', class(x))
    isect[]
}

#' @title shallow copy a \code{data.table}
#' 
#' @description Convenience function to shallow copy a \code{data.table} 
#' (until this function is exported in the \code{data.table} package). For 
#' internal use only.
#' @param x A \code{data.table}.
#' @param cols Character vector of column names (from \code{x}).
#' @return A shallow copied \code{data.table}.
#' @examples
#' \dontrun{
#' # For internal use only
#' library(data.table)
#' x <- data.table(a=1:2, b=3:4)
#' setattr(x, 'class', c("tmp", class(x)))
#'
#' y <- gread:::shallow(x) # only copies column pointers
#' class(y) # class(x) is retained
#' }
shallow <- function(x, cols = names(x)) {
    stopifnot(is.data.table(x), all(cols %in% names(x)))
    ans = vector("list", length(cols))
    setattr(ans, 'names', data.table::copy(cols))
    for (col in cols)
        ans[[col]] = x[[col]]
    setDT(ans)
    setattr(ans, 'class', data.table::copy(class(x)))
    ans[]
}

#' @title convert \code{GAlignments} object to \code{data.table}
#' 
#' @description For internal use only. Converts a \code{GAlignments} object to 
#' \code{bam} object.
#' 
#' @param x A \code{GAlignments} object.
#' @return An object of class \code{bam} that inherits from \code{data.table}.
#' @examples
#' \dontrun{
#' For internal use only
#' library(GenomicAlignments)
#' path <- system.file("tests", package="gread")
#' bam_file <- file.path(path, "sample.bam")
#' 
#' # with no metadata
#' bam <- GenomicAlignments::readGAlignments(bam_file)
#' gread:::as_bam(bam)
#' 
#' # with some metadata
#' bam <- as(as_granges(read_format(bam_file)), "GAlignments")
#' gread:::as_bam(bam)
#' }
as_bam <- function(x) {
    stopifnot(inherits(x, "GAlignments"))
    ans = list(
            seqname=as.character(seqnames(x)), 
            start=start(x), 
            end=end(x), 
            strand=as.character(strand(x)),
            cigar=as.character(cigar(x)),
            qwidth=qwidth(x),
            width=width(x),
            njunc=njunc(x))
    setDT(ans)
    if (ncol(mcols(x))) {
        mdata = as.data.table(mcols(x))
        ans = ans[, names(mdata) := mdata]
    }
    setattr(ans, 'class', c("bam", "data.table", "data.frame"))
    ans[]
}

# #' @title Convert to TxDb object
# #' 
# #' @description This is another helper function which allows extracting 
# #' features -- \code{genes}, \code{transcripts}, \code{exonsBy} genes, 
# #' transcripts etc. For internal use only. 
# #' 
# #' See \code{\link{transcriptsBy}} from the \code{GenomicFeatures} package for 
# #' more. 
# #' 
# #' @param x An object of class \code{gtf}, \code{gff}, \code{bed} or 
# #' \code{bam}.
# #' @param ignore_strand Logical argument to pass to \code{GRanges} function. 
# #' Indicates whether \code{strand} should be ignored when constructing 
# #' \code{GRanges} object or not. Default is \code{FALSE}.
# #' @return A \code{TxDb} object.
# #' @aliases as_txdb
# #' @seealso \code{\link{as_granges}} \code{\link{read_format}} 
# #' \code{\link{tidy}} \code{\link{extract}}
# #' @export
# #' @examples
# #' path <- system.file("tests", package="gread")
# #' gff_file <- file.path(path, "sample.gff")
# #' gtf_file <- file.path(path, "sample.gtf")
# #' 
# #' gff <- read_format(gff_file)
# #' gtf <- read_format(gtf_file)
# #' 
# #' as_txdb(gff)
# #' as_txdb(gtf)
# #' 
# #' as_txdb(gff, ignore_strand=FALSE)
# #' as_txdb(gtf, ignore_strand=FALSE)
# as_txdb <- function(x, ignore_strand=FALSE) {
#     stopifnot(is.gtf(x)||is.gff(x))
#     # for 'makeTxDbFromGRanges' if necessary to extract features via 
#     # bioc functions like exonsBy
#     if ("feature" %in% names(x)) x[, "type" := feature] 
#     else x[, c("type", "feature") := "reads"]
#     GenomicFeatures::makeTxDbFromGRanges(as_granges(x, ignore_strand))
# }
