tidy_cols <- function(x, ...) {
    UseMethod("tidy_cols")
}

tidy_cols.default <- function(x, ...) {
    stop("No default method available.")
}

tidy_cols.gtf <- function(x, remove_cols="attributes", verbose=FALSE, ...) {
    stopifnot(identical(head(names(x), 9L), format_names("gtf")))
    meta_fun <- function(vec, attr) {
        data.table::setattr(as.list(vec), 'names', attr)
    }
    if (ncol(x) > 9L) { # rtracklayer::readgff output, already tidied
        for (col in which(vapply(x, is.character, TRUE)))
            set(x, i=which(x[[col]] == ""), j=col, value=NA)
    } else if ("attributes" %in% names(x)) {
        if (verbose) cat("Tidying up 'attributes' column.\n")
        # base R's strsplit with perl=TRUE is faster than without, 
        # but for 'gtf', is twice as slow as stringi (20 vs 10 s)
        prep = stringi::stri_split_regex(x$attributes, 
                  "[ ]\"|\";[ ]*", omit_empty=TRUE)
        cols = lapply(prep, `[`, c(TRUE, FALSE))
        vals = lapply(prep, `[`, c(FALSE, TRUE))
        len1 = vapply(cols, length, 0L)
        len2 = vapply(vals, length, 0L)
        if (!identical(len1, len2)) {
            stop("Ill-formed gff file. Attributes column in gff file must be 
                of the form 'id1=val1;id2=val2;...")
        }
        prep = setDT(list(id=rep.int(seq_along(len1), len1), 
                          cols=unlist(cols), vals=unlist(vals)))
        meta = dcast(prep, id ~ cols, value.var="vals")
        set(meta, j="id", value=NULL)
        # x = data.table:::shallow(x) # TODO: use when exported
        x[, names(meta) := meta]
    }
    remove_cols(x, remove_cols, verbose) # updates by reference
    if (!all(c("transcript_id", "gene_id") %in% names(x)))
        warning("Columns 'transcript_id' and/or 'gene_id' are not found in ", 
            "the gtf file or has been removed while tidying. Note that ", 
            "these columns are essential for downstream analysis in ", 
            "most cases.")
    invisible(x)
}

tidy_cols.gff <- function(x, remove_cols="attributes", verbose=FALSE, ...) {
    stopifnot(identical(head(names(x), 9L), format_names("gff")))
    meta_fun <- function(vec, attr) {
        data.table::setattr(as.list(vec), 'names', attr)
    }
    if (ncol(x) > 9L) { # rtracklayer::readgff output, already tidied
        for (col in which(vapply(x, is.character, TRUE)))
            set(x, i=which(x[[col]] == ""), j=col, value=NA)
    } else if ("attributes" %in% names(x)) {
        if (verbose) cat("Tidying up 'attributes' column.\n")
        # base R's strsplit with perl=TRUE is faster than without, 
        # but for 'gtf', is twice as slow as stringi (20 vs 10 s)
        prep = stringi::stri_split_regex(x$attributes, "=|;", omit_empty=TRUE)
        cols = lapply(prep, `[`, c(TRUE, FALSE))
        vals = lapply(prep, `[`, c(FALSE, TRUE))
        len1 = vapply(cols, length, 0L)
        len2 = vapply(vals, length, 0L)
        if (!identical(len1, len2)) {
            stop("Ill-formed gff file. Attributes column in gff file must be 
                    of the form 'id1=val1;id2=val2;...")
        }
        prep = setDT(list(id=rep.int(seq_along(len1), len1), 
                          cols=unlist(cols), vals=unlist(vals)))
        meta = dcast(prep, id ~ cols, value.var="vals")
        set(meta, j="id", value=NULL)
        # x = data.table:::shallow(x) # TODO: use when exported
        x[, names(meta) := meta]
    }
    gff_gene_transcript_cols(x) # add transcript_id, gene_id cols
    # and if possible transcript_name and gene_name cols as well
    remove_cols(x, remove_cols, verbose) # updates by reference
    invisible(x)
}

tidy_cols.bed <- function(x, remove_cols=NULL, verbose=FALSE, ...) {
    stopifnot(identical(names(x), head(format_names("bed"), length(x))))
    remove_cols(x, remove_cols, verbose) # updates by reference
    invisible(x)
}

tidy_cols.bam <- function(x, remove_cols=NULL, verbose=FALSE, ...) {
    remove_cols(x, remove_cols, verbose) # updates by reference
    invisible(x)
}

## internal helper functions -------------------------------------------------

remove_cols <- function(x, cols, verbose) {
    if (is.null(cols)) {
        if (verbose) cat("remove_cols is NULL. No columns are removed.")
    } else {
        invalid = cols[!cols %in% names(x)]
        if (length(invalid)) {
            warning("Skipping column(s) [", paste(invalid, 
                collapse=", "), "] as they are not present in 'x'.")
            cols = setdiff(cols, invalid)
        }
        if (length(cols)) {
            if (verbose) cat("Removing [", paste(cols, collapse=", "), 
                "] columns.\n", sep="")
            set(x, j=cols, value=NULL)
        }
    }
    x[]
}

gff_gene_transcript_cols <- function(x) {
    # to please R CMD CHECK
    gene_id=transcript_id=NULL
    cols = c("ID", "Parent")
    if (any(check <- !cols %in% names(x))) {
        warning("gff file doesn't contain ", paste(cols[check], collapse=","), 
        "column(s) in its attributes column. They are necessary for ", 
        " creating 'transcript_id' and 'gene_id' columns which will be", 
        " required for most downstream analysis. Please fix the gff file.")
        cols = cols[!check]
        x[, c(cols) := lapply(.SD, function(x) gsub("^.*:", "", x)), 
                    .SDcols=cols]
    } else {
        ID=Parent=i.Parent=Name=i.Name=NULL # suppress R CMD CHECK warnings
        x[, "ID" := gsub("^.*:", "", ID)]
        x[, "Parent" := gsub("^.*:", "", Parent)]
        # transcript_id
        rnaids = grepl("rna", x$feature, ignore.case=TRUE)
        x[, "transcript_id" := Parent][rnaids, "transcript_id" := ID]
        # gene_id
        tmp = x[rnaids, list(transcript_id=ID, Parent)]
        if (anyDuplicated(tmp)) tmp = unique(tmp, by=names(tmp))
        x[, gene_id := ID][tmp, gene_id := i.Parent, on="transcript_id"]
        if ("Name" %in% names(x)) {
            # transcript_name
            tmp = x[rnaids, list(transcript_id, Name)]
            if (anyDuplicated(tmp)) tmp = unique(tmp, by=names(tmp))
            i.Name = NULL
            x[tmp, "transcript_name" := i.Name, on="transcript_id"]
            # gene_name
            geneids = grepl("gene", x$feature, ignore.case=TRUE)
            tmp = x[geneids, list(gene_id, Name)]
            if (anyDuplicated(tmp)) tmp = unique(tmp, by=names(tmp))
            x[tmp, "gene_name" := i.Name, on="gene_id"]
        }
    }
}
