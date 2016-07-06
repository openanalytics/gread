#' @title Currently supported formats
#'
#' @description Maintains a vector of all the formats supported currently.
#'
#' @return A character vector listing all currently supported formats.
#' @aliases supported_formats
#' @seealso \code{\link{read_format}} \code{\link{tidy}}
#' @examples
#' supported_formats()
#' @export
supported_formats <- function() {
    c("gtf", "gff", "bed", "bam")
}

error_format <- function() {
    stop("'format' must be one of ", paste(supported_formats(), collapse="/"))
}
