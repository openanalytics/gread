#' @title Extract GTF coordinates intersecting input BED file
#' 
#' @description Given a \code{GTF} file, and a \code{BED} file, this function 
#' extracts all \emph{transcripts} that overlap with coordinates from 
#' \code{BED} file and returns a \emph{filtered} \code{GTF} file. This 
#' function is at the moment only for internal purposes.
#' 
#' @param gtf_file Full path to a \code{GTF} file.
#' @param bed_file Full path to a \code{BED} file. The file must contain at 
#' least three columns, with them corresponding to \code{seqname} (or 
#' \code{chr}), \code{start} and \code{end}.
#' @param select_features A 1-column \code{data.table} or a named list of 
#' length=1. The name of the list indicates the column in the \code{GTF} file 
#' to filter on. The column/value should be a \emph{character} vector 
#' containing the values to \emph{retain}.
#' @seealso \code{\link{read_format}}, \code{\link{non_overlaps}}, 
#' \code{\link{construct_introns}}, \code{\link{extract}}
#' @return A \code{gtf} object containing just those transcripts that overlap 
#' with the provided bed file.
intersect_bed <- function(gtf_file, bed_file, select_features) {

    # to please R CMD CHECK
    feature=seqname=transcript_id=NULL
    # load seqname, start and end columns
    bed = read_format(bed_file)
    setnames(bed, names(bed)[1:3], c("seqname", "start", "end"))
    bed[, "start" := start+1L]
    setkey(bed)

    gtf = read_format(gtf_file)
    if (!missing(select_features)) {
        select_features = unique(as.data.table(select_features))
        gtf = gtf[select_features, on=names(select_features)[1L], nomatch=0L]
    }
    bed_chrs = unique(bed$seqname)
    gtf_chrs = unique(gtf$seqname)
    if (length(diff_chrs <- setdiff(bed_chrs, gtf_chrs))) {
        if (length(diff_chrs) == length(bed_chrs))
            stop("Chromosomes in BED file are completely different from 
                that of GTF file.")
        else warning("Chromosomes ", paste(diff_chrs, collapse=","), " are 
                found in BED but not in GTF file.")
    }

    construct_transcripts <- function(x) {
        x[feature %chin% "exon", 
            .(seqname = seqname[1L], start = min(start), 
                end = max(end), strand = strand[1L]), 
        by = "transcript_id"]
    }
    transcripts = construct_transcripts(gtf)

    olaps = foverlaps(transcripts, bed, type="any", which=TRUE, nomatch=0L)
    tids  = transcripts[olaps$xid, unique(transcript_id)]
    gtf[.(tids), on="transcript_id", nomatch=0L]
}
