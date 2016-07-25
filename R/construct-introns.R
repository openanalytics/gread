#' @title Construct introns from gtf/gff objects
#'
#' @description This function generates intronic coordinates by extracting 
#' all the \code{exons} from a \code{gtf} or \code{gff} object.
#' 
#' @param x An object of class \code{gtf/gff}. It has to have \code{exon} 
#' feature present.
#' @param update If \code{TRUE} (default), \code{x} is updated by reference 
#' and returned invisibily, else just the \code{intron} coordinates are 
#' returned.
#' @return An object of class \code{gtf/gff} with updated \code{intron} 
#' coordinates or just the intron coordinates depending on \code{update}. 
#' They all inherit from \code{GRanges}.
#' @export
#' @examples
#' path <- system.file("tests", package="gread")
#' gtf_file <- file.path(path, "sample.gtf")
#' gtf <- read_format(gtf_file)
#' # update gtf with intron coordinats from the exon gaps 
#' # and return the updated object
#' ans <- construct_introns(gtf, update=TRUE)[] # default
#' # same as above, but returns only constructed introns
#' introns <- construct_introns(gtf, update=FALSE)
#' @seealso \code{\link{supported_formats}} \code{\link{read_format}} 
#' \code{\link{extract}} \code{\link{as_granges}} 
construct_introns <- function(x, update=TRUE) {
    stopifnot(is.gtf(x) || is.gff(x))
    x = as_data_table(x)
    stopifnot("feature" %chin% names(x), 
                "transcript_id" %chin% names(x), 
                update %in% c(TRUE, FALSE))

    # to please R CMD CHECK
    feature=seqnames=transcript_id=NULL
    x_class = copy(class(x))
    exons = x[feature == "exon"]
    if (!nrow(exons)) stop("feature == 'exon' returned 0 rows.")
    setorderv(exons, c("transcript_id", "seqnames", "start", "end", "strand"))
    introns = exons[, .(seqnames = seqnames[1L], 
                        start = end[seq_len(.N-1L)]+1L, 
                        end = start[seq_len(.N-1L)+1L]-1L,
                        feature = "intron",
                        strand = strand[1L]), by=transcript_id]
    introns = na.omit(introns, cols = c("start", "end"))
    check = introns[, any(start > end), by=transcript_id]
    if (length(ids <- which(check[["V1"]])))
        stop("Exons for transcript ids [", paste(ids, collaspse=" "), 
            "] have start > end")
    exons[, c(setdiff(names(introns), "transcript_id")) := NULL]
    exons = unique(shallow(exons, reset_class=TRUE), by="transcript_id")
    ecols = names(exons)
    introns[exons, (ecols) := mget(ecols), on="transcript_id"]
    # reset all exon related cols
    introns[, grep("^exon", names(introns), value=TRUE) := NA]
    colorder = names(x)
    if (update) {
        x = rbind(x, introns)
        setorderv(x, c("seqnames", "start", "end"))
        x = x[, .SD, by="transcript_id"]
        setattr(x, 'class', x_class)
        setcolorder(x, colorder)
        ans = x
    } else {
        setcolorder(introns, colorder)
        setattr(introns, 'class', c("intron", "data.table", "data.frame"))
        ans = introns
    }
    new(class(ans)[1L], as(setDF(ans), "GRanges"))
}
