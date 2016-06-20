.onLoad <- function(libname, pkgname) {
    gread_verbose = "FALSE"
    # global options
    opts = c(
              "gread_verbose" = gread_verbose
            )
    for (i in setdiff(names(opts), names(options())) ) {
        text = paste('options(', i, '=', opts[i], ')', sep="")
        eval(parse(text=text))
    }
    invisible()
}
