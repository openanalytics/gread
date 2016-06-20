#' @title Construct introns from \code{gtf} or \code{gff} objects
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
#' @export
#' @examples
#' path = system.file("tests", package="gread")
#' gtf_file = file.path(path, "sample.gtf")
#' gtf = read_format(gtf_file)
#' # update gtf with intron coordinats from the exon gaps 
#' # and return the updated object
#' ans = construct_introns(gtf, update=TRUE)[] # default
#' # same as above, but returns only constructed introns
#' introns = construct_introns(gtf, update=FALSE)
#' @seealso \code{\link{supported_formats}} \code{\link{read_format}} 
#' \code{\link{extract}} \code{\link{tidy}} \code{\link{as_granges}} 
construct_introns <- function(x, update=TRUE) {
    stopifnot(inherits(x, "gtf") || inherits(x, "gff"), 
                "feature" %chin% names(x), 
                "transcript_id" %chin% names(x), 
                update %in% c(TRUE, FALSE))

    # to please R CMD CHECK
    feature=seqname=transcript_id=NULL
    x_class = copy(class(x))
    exons = x[feature == "exon"]
    if (!nrow(exons)) stop("feature == 'exon' returned 0 rows.")
    setorderv(exons, c("transcript_id", "seqname", "start", "end", "strand"))
    introns = exons[, .(seqname = seqname[1L], 
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
    exons = unique(exons, by="transcript_id")
    ecols = names(exons)
    introns[exons, (ecols) := mget(ecols), on="transcript_id"]
    # reset all exon related cols
    introns[, grep("^exon", names(introns), value=TRUE) := NA]
    colorder = names(x)
    if (update) {
        x = rbind(x, introns)
        setorderv(x, c("seqname", "start", "end"))
        x = x[, .SD, by="transcript_id"]
        setattr(x, 'class', x_class)
        setcolorder(x, colorder)
        return(x)
    } else {
        setcolorder(introns, colorder)
            setattr(introns, 'class', unique(c("intron", x_class)))
        return(introns[])
    }
}
