.onAttach <- function(libname, pkgname) {
    # Runs when attached to search() path such as by library() or require()
    if (interactive()) {
        packageStartupMessage('v', as.character(packageVersion("gread")), 
          ', type vignette("gread-vignette", package="gread") to get started.')
    }
}
