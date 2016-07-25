gread <- function(token) {
    UseMethod("gread")
}

gread.default <- function(token) {
    stop("No method available for format ", gsub("(.*)_format$", 
            "", class(token)))
}

gread.gtf_format <- function(token) {
    # to please R CMD CHECK
    score=phase=NULL
    # rtracklayer::readGFF reads gff v1 'attributes' column as 'group',
    # but we want it to be always 'attributes'. Also, it tidies up the 
    # attributes column by default
    rtracklayer_fun <- function() {
        ans = setDT(rtracklayer::readGFF(token$file, version = 2L, 
                columns = rtracklayer::GFFcolnames()))
        ans[, (1:3) := lapply(.SD, as.character), .SDcols=1:3][]
        # remove additional attributes
        setattr(ans, 'pragmas', NULL)
        setattr(ans, 'attrcol_fmt', NULL)
        setattr(ans, 'ncol0', NULL)
        setattr(ans, 'ntag', NULL)
        setattr(ans, 'raw_data', NULL)
        ans
    }
    fread_fun <- function() fread(token$file, colClasses = token$types, 
        showProgress = token$verbose)
    read_table_fun <- function() read_table(token$file, sep = "\t", 
        header = FALSE, comment.char = "#", nrows = -1L, colClasses = 
        token$types, quote = "")
    ans = tryCatch(rtracklayer_fun(), error = function(o) { 
                if (token$verbose) cat("rtracklayer::readGFF failed to read", 
                    " the gtf file. Reverting to data.table::fread.\n", sep="")
                tryCatch(fread_fun(), error = function(o) {
                    if (token$verbose) cat("data.table::fread failed to read", 
                        " as well. Reverting to base::read.table.\n", sep="")
                    read_table_fun()
                })
          })
    setnames(ans, head(names(ans), 9L), token$names)
    if (is.character(ans[["score"]])) 
        ans[, "score" := suppressWarnings(as.numeric(score))]
    if (is.character(ans[["phase"]])) 
        ans[, "phase" := suppressWarnings(as.integer(phase))]
    setattr(ans, 'class', c("gtf", "data.table", "data.frame"))
}

gread.gff_format <- function(token) {
    # to please R CMD CHECK
    score=phase=NULL
    # rtracklayer::readGFF reads gff v1 'attributes' column as 'group',
    # but we want it to be always 'attributes'. Also, it tidies up the 
    # attributes column by default
    rtracklayer_fun <- function() {
        ans = setDT(as.data.frame(rtracklayer::readGFF(token$file, columns = 
                rtracklayer::GFFcolnames())))
        if (!token$tidy_cols) set(ans, j = tail(names(ans), -9L), value = NULL)
        ans[, (1:3) := lapply(.SD, as.character), .SDcols=1:3][]
        list_cols = which(vapply(ans, is.list, TRUE))
        if (length(list_cols)) {
            list_paste <- function(x) sapply(x, paste, collapse=",")
            # paste all values of list cols together
            ans[, (list_cols) := lapply(.SD, list_paste), .SDcols=list_cols]
        }
        # as.data.frame drops all extra attributes
        ans
    }
    fread_fun <- function() fread(token$file, colClasses = token$types, 
        showProgress = token$verbose)
    read_table_fun <- function() read_table(token$file, sep = "\t", 
        header = FALSE, comment.char = "#", nrows = -1L, colClasses = 
        token$types, quote = "")
    ans = tryCatch(rtracklayer_fun(), error = function(o) { 
                if (token$verbose) cat("rtracklayer::readGFF failed to read",
                    " the gtf file. Reverting to data.table::fread.\n", sep="")
                tryCatch(fread_fun(), error = function(o) {
                    if (token$verbose) cat("data.table::fread failed to read", 
                        " as well. Reverting to base::read.table.\n", sep="")
                    read_table_fun()
                })
          })
    setnames(ans, head(names(ans), 9L), token$names)
    if (is.character(ans[["score"]])) 
        ans[, "score" := suppressWarnings(as.numeric(score))]
    if (is.character(ans[["phase"]])) 
        ans[, "phase" := suppressWarnings(as.integer(phase))]
    setattr(ans, 'class', c("gff", "data.table", "data.frame"))
}

gread.bed_format <- function(token) {
    ans = tryCatch(fread(token$file, colClasses = token$types, sep="\t", 
                showProgress = token$verbose), error = function(o) {
                    if (token$verbose) cat("data.table::fread failed to read", 
                        " the bed file. Reverting to 'read.table'", 
                        " from base R.\n", sep="")
                    read_table(token$file, sep = "\t", header = FALSE, 
                        comment.char = "#", nrows = -1L, colClasses = 
                        token$types, quote = "")
                    })
    # setnames will error if bed file has < 3 cols
    setnames(ans, head(token$names, max(3L, ncol(ans))))
    # ans[, start := start+1L] # TODO: should we automatically set start+1L?
    setattr(ans, 'class', c("bed", "data.table", "data.frame"))
}

gread.bam_format <- function(token) {
    if (token$verbose) 
        cat("Loading bam file: ",token$file,", please wait...\n", sep="")
    stopifnot(length(token$file)==1L)
    
    what  = c("flag")
    flags = Rsamtools::scanBamFlag(isUnmappedQuery=FALSE)
    if (!length(token$chromosomes)) {
        param = Rsamtools::ScanBamParam(what=what, tag=token$tags, flag=flags)
    } else {
            hdr = Rsamtools::scanBamHeader(token$file)[[1]]$targets
            hdr = hdr[names(hdr) %in% token$chromosomes]
            if (!all(names(hdr) %in% token$chromosomes)) {
                stop("chromosomes [", paste(setdiff(token$chromosomes, 
                    names(hdr)), collapse=", "), "] not found in bam file.")
            }
            which = as_granges(data.table(seqnames = names(hdr), 
                        start=1L, end = hdr, strand="*"), ignore_strand=TRUE)
            param = Rsamtools::ScanBamParam(what = what, which = which, 
                                tag = token$tags, flag = flags)
    }
    # load bam with corresponding param
    ans = as_bam(GenomicAlignments::readGAlignments(token$file, param=param))
    setattr(ans, 'class', c("bam", "data.table", "data.frame"))
}

# Internal helper functions --------------------------------------------------
read_table <- function(file, ...) {
    setDT(read.table(file, stringsAsFactors=FALSE, as.is=TRUE, ...))
}
