#' @title Fast Reading and Processing  of Common Gene Annotation and Next 
#' Generation Sequencing Format Files.
#' 
#' @description The \code{gread} package allows to read and process common 
#' gene annotation and next generation sequencing file formats easily. It 
#' currently supports:
#' 
#' \itemize{
#'  \item Gene Feature Format or \code{GFF} using \code{rtracklayer} package,
#'  \item Gene Transfer Format or \code{GTF} (also referred to as 
#'     \code{GFF 2.5 at times}) using \code{rtracklayer} package,
#'  \item Browser Extensible Data or \code{BED} format using the 
#'     \code{data.table} package.
#'  \item Binary sequence Alignment/Map Format or \code{BAM} format using 
#'     \code{GenomicAlignments} package.
#' }
#' 
#' In addition, it also provides additional utility functions such as:
#' 
#' \itemize{
#'  \item Filling missing intron coordinates - \code{construct_introns}, 
#'  \item Extracting particular features from GFF/GTF files - \code{extract}, 
#'  \item Extract non-overlapping features (e.g., \code{genes}, 
#'     \code{introns}) etc. - \code{non_overlaps}
#' }
#' 
#' @references \url{http://gmod.org/wiki/GFF3} \cr
#' \url{http://www.ensembl.org/info/website/upload/gff.html} \cr
#' \url{https://github.com/rdatatable/data.table} \cr
#' \url{https://bioconductor.org} \cr
#' @aliases gread gread-package
#' @docType package
#' @name gread
#' @importFrom tools file_ext
#' @importFrom stringi stri_split_regex
#' @import data.table
#' @import methods
#' @importFrom utils capture.output head packageVersion read.table tail
#' @importFrom stats na.omit
#' @importFrom Rsamtools ScanBamParam scanBamFlag scanBamHeader
#' @importFrom GenomicRanges GRanges split
#' @importFrom IRanges IRanges
#' @importMethodsFrom GenomicRanges disjoin reduce intersect  findOverlaps 
#' @importMethodsFrom GenomicRanges countOverlaps seqnames start end strand
#' @import GenomicFeatures
#' @importFrom GenomicAlignments findOverlaps summarizeOverlaps cigar qwidth
#' @importFrom GenomicAlignments readGAlignments readGAlignmentPairs 
#' @importFrom GenomicAlignments width njunc rname
#' @importFrom S4Vectors queryHits subjectHits
#' @importMethodsFrom S4Vectors Rle elementMetadata
#' @importFrom rtracklayer readGFF GFFcolnames
#' @seealso \code{\link{read_format}} \code{\link{read_gff}} 
#' \code{\link{read_gtf}} \code{\link{read_bed}} \code{\link{read_bam}} 
#' \code{\link{extract}} \code{\link{construct_introns}}
#' @keywords file
NULL
