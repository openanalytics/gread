#' @title Extract gtf coordinates intersecting input bed file
#' 
#' @description Given a \code{gtf} file, and a \code{bed} file, this function 
#' extracts all \emph{transcripts} that overlap with coordinates from 
#' \code{bed} file and returns a \emph{filtered} \code{gtf} file. This 
#' function is at the moment only for internal purposes.
#' 
#' @param gtf_file Full path to a \code{gtf} file.
#' @param bed_file Full path to a \code{bed} file. The file must contain at 
#' least three columns, with them corresponding to \code{seqnames} (or 
#' \code{chr}), \code{start} and \code{end}.
#' @param select_features A 1-column \code{data.table} or a named list of 
#' length=1. The name of the list indicates the column in the \code{gtf} file 
#' to filter on. The column/value should be a \emph{character} vector 
#' containing the values to \emph{retain}.
#' @seealso \code{\link{read_format}}, \code{\link{non_overlaps}}, 
#' \code{\link{construct_introns}}, \code{\link{extract}}
#' @return A \code{gtf} object that inherits from \code{GRanges} containing 
#' just those transcripts that overlap with the provided bed file.
intersect_bed <- function(gtf_file, bed_file, select_features) {

    # to please R CMD CHECK
    feature=seqnames=transcript_id=NULL
    # load seqnames, start and end columns
    bed = read_format(bed_file)
    setnames(bed, names(bed)[1:3], c("seqnames", "start", "end"))
    bed[, "start" := start+1L]
    setkey(bed)

    gtf = read_format(gtf_file)
    if (!missing(select_features)) {
        select_features = unique(as.data.table(select_features))
        gtf = gtf[select_features, on=names(select_features)[1L], nomatch=0L]
    }
    bed_chrs = unique(bed$seqnames)
    gtf_chrs = unique(gtf$seqnames)
    if (length(diff_chrs <- setdiff(bed_chrs, gtf_chrs))) {
        if (length(diff_chrs) == length(bed_chrs))
            stop("Chromosomes in bed file are completely different from 
                that of gtf file.")
        else warning("Chromosomes ", paste(diff_chrs, collapse=","), " are 
                found in bed but not in gtf file.")
    }

    construct_transcripts <- function(x) {
        x[feature %chin% "exon", 
            .(seqnames = seqnames[1L], start = min(start), 
                end = max(end), strand = strand[1L]), 
        by = "transcript_id"]
    }
    transcripts = construct_transcripts(gtf)

    olaps = foverlaps(transcripts, bed, type="any", which=TRUE, nomatch=0L)
    tids  = transcripts[olaps$xid, unique(transcript_id)]
    ans = gtf[.(tids), on="transcript_id", nomatch=0L]
    new("gtf", as(setDF(ans), "GRanges"))
}
