#' @title Run a set of tests
#' 
#' @description Runs a set of tests to check if \code{gread} is working 
#' correctly.
#'
#' @details Runs a series of tests. These can be used to see features and 
#' examples of usage, too. Running \code{test_gread} will tell you the full 
#' location of the test file(s) to open.
#' @param verbose If \code{TRUE} sets \code{gread_verbose} to \code{TRUE} for 
#' the duration of the tests.
#' @param pkg Root directory name under which all package content 
#' (ex: \code{DESCRIPTION}, \code{src/}, \code{R/}, \code{inst/} etc..) 
#' resides.
#' @param silent Logical, default \code{FALSE}. When \code{TRUE} it will not 
#' raise error on in case of test fails.
#' @return When silent equals to \code{TRUE} it will return \code{TRUE} if 
#' all tests were successful. \code{FALSE} otherwise. If silent equals to 
#' \code{FALSE} it will return \code{TRUE} if all tests were successful. Error 
#' otherwise.
#' @aliases test_gread
#' @seealso \code{\link{read_format}} \code{\link{extract}} 
#' \code{\link{tidy_cols}}
#' @examples
#' \dontrun{
#' gread:::test_gread()
#' }
test_gread <- function(verbose=FALSE, pkg="pkg", silent=FALSE) {
    if (exists("test_gread", .GlobalEnv, inherits=FALSE)) { # package developer
        if ("package:gread" %in% search()) stop("gread package loaded")
        if (.Platform$OS.type == "unix") {
            d = path.expand("~/Documents/oa-git/gread/inst/tests")
        }
    } else { # user
        d = paste(getNamespaceInfo("gread", "path"), "/tests", sep="")
    }
    olddir = setwd(d)
    on.exit(setwd(olddir))
    envirs <- list()
    for (fn in file.path(d, 'tests.Rraw')) {    # not testthat
        cat("Running", fn, "\n")
        oldverbose = getOption("gread_verbose")
        if (verbose) options(gread_verbose=TRUE)
        envirs[[fn]] = new.env(parent=.GlobalEnv)
        if(isTRUE(silent)){
            try(sys.source(fn, envir=envirs[[fn]]), silent=silent)
        } else {
            sys.source(fn, envir=envirs[[fn]])
        }
        options(gread_verbose=oldverbose)
    }
    invisible(sum(sapply(envirs, `[[`, "nfail"))==0)
}

# Define test() and its globals here, for use in dev
# But primarily called by test_gread() calling inst/tests/tests.Rraw
# Initialized at the top of tests.Raw ...
# nfail = ntest = lastnum = 0
# whichfail = NULL
# .devtesting = TRUE

compactprint <- function(DT, topn=2) {
    cn = paste(" [Key=", paste(key(DT), collapse=","), " Types=", 
         paste(substring(gsub("integer64", "i64", sapply(DT, class)), 1, 3), 
           collapse=","), "]", sep="")
    print(copy(DT)[, (cn):=""], topn=topn)
    invisible()
}

test <- function(num, x, y, error=NULL, warning=NULL, output=NULL) {
    # Usage:
    # 1. tests that x equals y when both x and y are supplied, the most 
    #    common usage
    # 2. tests that x is TRUE when y isn't supplied
    # 3. if error is supplied, y should be missing and x is tested to result 
    #    in an error message matching the pattern
    # 4. if warning is supplied, y (if supplied) is checked to equal x, and 
    #    x should result in a warning message matching the pattern
    # 5. if output is supplied, x is evaluated and printed and the output is 
    #    checked to match the pattern
    # At most one of error, warning or output may be supplied, all single 
    # character strings (passed to grep)
    # num just needs to be numeric and unique. We normally increment integers 
    # at the end, but inserts can be made using decimals e.g. 10, 11, 11.1, 
    # 11.2, 12, 13,...
    # Motivations:
    # 1. we'd like to know all tests that fail not just stop at the first. 
    #    This often helps by revealing a common feature across a set of
    #    failing tests
    # 2. test() tests more deeply than a diff on console output and uses a 
    #    gread appropriate definition of "equals" different
    #    from all.equal and different to identical related to row.names and 
    #    unused factor levels
    # 3. each test has a unique id which we refer to in commit messages, 
    #    emails etc.

    # to cater for both test_gread() and stepping through tests in dev
    nfail = get("nfail", parent.frame()) 
    whichfail = get("whichfail", parent.frame())
    assign("ntest", get("ntest", parent.frame())+1L, parent.frame(), 
            inherits=TRUE)   # bump number of tests run
    assign("lastnum", num, parent.frame(), inherits=TRUE)
    v = getOption("gread_verbose")
    i = interactive()
    if (v || i) {
        cat(if (i) "\r" else "\n\n", "Running test id ", num, "     ", sep="")
    }
    # TODO: every line that could possibly fail should ideally be 
    # inside test()
    xsub = substitute(x)
    ysub = substitute(y)
    if (is.null(output)) err <<- try(x, TRUE)
    else out = gsub("NULL$", "", paste(capture.output(print(err<<-try(x, 
               TRUE))), collapse=""))
    if (!is.null(output)) {
        output = gsub("[[]", "<LBRACKET>", output)
        output = gsub("[]]", "<RBRACKET>", output)
        output = gsub("<LBRACKET>", "[[]", output)
        output = gsub("<RBRACKET>", "[]]", output)
        output = gsub("[(]", "[(]", output)
        output = gsub("[)]", "[)]", output)
        if (!length(grep(output, out))) {
            cat("Test", num, "didn't produce correct output:\n")
            cat(">", deparse(xsub), "\n")
            cat("Expected: '", output, "'\n", sep="")
            cat("Observed: '", out, "'\n", sep="")
            assign("nfail", nfail+1, parent.frame(), inherits=TRUE)
            assign("whichfail", c(whichfail, num), parent.frame(), 
                   inherits=TRUE)
            return()
        }
    }
    if (!is.null(error) || !is.null(warning)) {
        type = ifelse(!is.null(error), "error", "warning")
        patt = txt = ifelse(!is.null(error), error, warning)
        patt = gsub("[(]", "[(]", patt)
        patt = gsub("[)]", "[)]", patt)
        patt = gsub("\\^", "\\\\^", patt)
        observedtype = ifelse(length(grep("converted from warning", err)), 
                       "warning", "error")
        if (! (inherits(err, "try-error") &&
               length(grep(patt, err)) &&
               type==observedtype)) {
            cat("Test", num, "didn't produce correct", type, ":\n")
            cat(">", deparse(xsub), "\n")
            cat("Expected ", type, ": '", txt, "'\n", sep="")
            if (!inherits(err, "try-error"))
                cat("Observed: no error or warning\n")
            else
                cat("Observed ", observedtype, ": '", 
                    gsub("^[(]converted from warning[)] ", "", 
                        gsub("\n$", "", 
                            gsub("^Error.* : \n  ", "", as.character(err)))), 
                    "'\n", sep="")
            assign("nfail", nfail+1L, parent.frame(), inherits=TRUE)
            # Not the same as nfail <<- nfail + 1, it seems 
            # (when run via R CMD check)
            assign("whichfail", c(whichfail, num), parent.frame(), 
                            inherits=TRUE)
            return()
        }
        if (type=="warning")
            err <- if (is.null(output)) x<-try(suppressWarnings(x), TRUE) 
                   else out <- paste(capture.output(x<-try(suppressWarnings(
                               x), TRUE)), collapse="")
        else return()
    }
    if (inherits(err, "try-error") || (!missing(y) && inherits(err<-try(y, 
        TRUE), "try-error"))) {
        cat("Test", num, err)
        assign("nfail", nfail+1, parent.frame(), inherits=TRUE)  
        assign("whichfail", c(whichfail, num), parent.frame(), inherits=TRUE)
        return()
    }
    if (missing(y)) {
        if (!is.null(output)) return()
        if (isTRUE(as.vector(x))) return()
        cat("Test", num, "expected TRUE but observed:\n")
        cat(">", deparse(xsub), "\n")
        if (is.data.table(x)) compactprint(x) else print(x)
        assign("nfail", nfail+1, parent.frame(), inherits=TRUE)
        assign("whichfail", c(whichfail, num), parent.frame(), inherits=TRUE)
        return()
    } else {
        if (identical(x, y)) return()
        if (is.data.table(x) && is.data.table(y)) {
            # TO DO:  test 166 doesn't pass with these :
            # if (!selfrefok(x)) stop("x selfref not ok")
            # if (!selfrefok(y)) stop("y selfref not ok")
            xc=copy(x)
            yc=copy(y)  # so we don't affect the original data which 
                        # may be used in the next test
            # drop unused levels in factors
            if (length(x)) for (i in which(sapply(x, is.factor))) {
                .xi=x[[i]]
                xc[, (i):=factor(.xi)]
            }
            if (length(y)) {
                for (i in which(sapply(y, is.factor))) {
                    .yi=y[[i]]
                    yc[, (i):=factor(.yi)]
                }
            }
            setattr(xc, "row.names", NULL)  
            setattr(yc, "row.names", NULL)
            setattr(xc, "index", NULL)
            setattr(yc, "index", NULL)
            if (identical(xc, yc) && identical(key(x), key(y))) return()
            if (isTRUE(all.equal(xc, yc)) && identical(key(x), key(y))) 
                return()
        }
        if (is.factor(x) && is.factor(y)) {
            x = factor(x)
            y = factor(y)
            if (identical(x, y)) return()
        }
        if (is.atomic(x) && is.atomic(y) && isTRUE(all.equal(x, y))) 
            return()
    }
    cat("Test", num, "ran without errors but failed check that x equals y:\n")
    cat("> x =", deparse(xsub), "\n")
    if (is.data.table(x)) compactprint(x) else if (length(x)>6) {
        cat("First 6 of", length(x), ":")
        print(head(x))
    } else print(x)
    cat("> y =", deparse(ysub), "\n")
    if (is.data.table(y)) compactprint(y) else if (length(y)>6) {
        cat("First 6 of", length(y), ":")
        print(head(y))
    } else print(y)
    assign("nfail", nfail+1L, parent.frame(), inherits=TRUE)
    assign("whichfail", c(whichfail, num), parent.frame(), inherits=TRUE)
    invisible()
}
