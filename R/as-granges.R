#' @title Convert data.table to GRanges object
#' 
#' @description Convert \code{gtf}, \code{gff}, \code{bed}, \code{bam} 
#' or a valid \code{data.table} to a GRanges object.
#'
#' @param x An object of class \code{gtf}, \code{gff}, \code{bed} or 
#' \code{bam} or a valid \code{data.table} object.
#' @param ignore_strand Logical argument to pass to \code{GRanges} function. 
#' Indicates whether \code{strand} should be ignored when constructing 
#' \code{GRanges} object or not. Default is \code{FALSE}.
#' @return A \code{GRanges} object.
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
#' @seealso \code{\link{read_format}} \code{\link{extract}} 
#' \code{\link{construct_introns}}
as_granges <- function(x, ignore_strand=FALSE) {

    stopifnot(is.gtf(x)||is.gff(x)||is.bed(x)||is.bam(x)||is.data.table(x))
    x = shallow(x)
    if (ignore_strand) x[, "strand" := NULL]
    as(setDF(x), "GRanges")
}
