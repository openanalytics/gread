#' @title Extract features from gtf/gff objects
#'
#' @description Provides functions for further post processing on objects of 
#' class \code{gtf} and \code{gff}.
#'
#' @details Extract features based on various criteria (usually intended for  
#' obtaining read counts using \code{gcount} for a given \code{bam} file.
#' 
#' @param x Input object of class \code{gtf} or \code{gff} which inherits 
#' from \code{GRanges}.
#' @param feature A character vector of (usually related) features to extract 
#' from. One of \code{"gene_exon"}, \code{"gene"}, \code{"gene_intron"}, 
#' \code{"exon"}, \code{"intron"}. NB: \code{"exon"} feature must be present 
#' in \code{x}.
#' @param type \code{default} just extracts the features and returns it as 
#' such.
#' 
#' \code{union} merges all overlapping intervals into one. For e.g., with 
#' intervals \code{[a,b], [c,d], [e,f]} where \code{c < a < e < d < b < f}, 
#' the \code{union} is \code{[c, f]}. NB: There may be more than one row per 
#' \code{feature}.
#' 
#'  \code{intersect} returns only the intersecting part. Using the same 
#' intervals as before, the intersection is \code{[e,d]}. NB: If there is an 
#' intersection, exactly one row is returned, else the \code{feature} is 
#' skipped entirely (0-rows).
#'            
#' \code{disjoin} splits intervals into non-overlapping pieces. Using 
#' the same interval as before, the pieces would be \code{[c,a-1]} and 
#' \code{[b+1,f]}. NB: it could result in multiple rows for each a given 
#' \code{feature}.
#' 
#' \code{longest} retains only the longest interval.
#' 
#' \code{shortest} retains only the shortest interval.
#' 
#' \code{overlap} is a special case. Of the overlapping intervals, only 
#' the shortest interval is retained iff they all have identical \code{start}, 
#' \code{end}, or both. If not, all overlapping intervals are retained. For 
#' e.g., with intervals \code{[a,b], [c,d], [e,f]} where \code{a == c, b == f, 
#' d > b,f and e > a,c}, the interval \code{[e,f]} will be retained.
#'
#' @param ignore_strand Logical argument to pass to \code{GRanges} function. 
#' Indicates whether \code{strand} should be ignored when constructing 
#' \code{GRanges} object or not. Default is \code{FALSE}.
#' @param transcript_id Column name in \code{x} corresponding to transcript 
#' id. Default value is \code{"transcript_id"}.
#' @param gene_id Column name in \code{x} corresponding to gene id. Default 
#' value is \code{"gene_id"}.
#' @param ... Arguments passed to other functions. Ignored at the moment.
#' @aliases extract extract_feature
#' @return An object of class \code{"gene"} when \code{feature} is 
#' \code{"gene"}, \code{"gene_exon"} or \code{"gene_intron"}, and of class 
#' \code{"exon"} and \code{"intron"} when \code{feature} is \code{"exon"} or 
#' \code{"intron"} respectively. They all inherit from \code{GRanges}.
#' @seealso \code{\link{read_format}} \code{\link{as_granges}} 
#' \code{\link{extract}} \code{\link{construct_introns}}
#' @export
#' @examples
#' path <- system.file("tests", package="gread")
#' gtf_file <- file.path(path, "sample.gtf")
#' gtf <- read_format(gtf_file)
#' # extract exons, combine coordinates of overlapping exons
#' exons <- extract(gtf, feature="exon", type="union")
#' # extract all exons within the gene, but combine overlapping exons
#' exons <- extract(gtf, feature="gene_exon", type="union")
#' ## extract gene span (uses exon coordinates if feature='gene' doesn't exist)
#' genes <- extract(gtf, feature="gene", type="default")
extract <- function(x, feature=c("gene_exon", "gene", "gene_intron", "exon", 
    "intron"), type=c("default", "union", "disjoin", "intersect", "longest", 
    "shortest", "overlap"), ignore_strand=FALSE, 
    transcript_id="transcript_id", gene_id="gene_id", ...) {
    feature = match.arg(feature)
    splits  = unlist(strsplit(feature, "_", fixed=TRUE))
    feature = splits[1L]
    subfeature = if (is.na(splits[2L])) NULL else splits[2L]
    stopifnot(is.gtf(x) || is.gff(x))
    x <- as_data_table(x)
    stopifnot("feature" %in% names(x), 
                ignore_strand %in% c(FALSE, TRUE), 
                transcript_id %in% names(x), 
                gene_id %in% names(x), nrow(x)>0L)
    setattr(x, 'class', c(match.arg(feature), class(x)))
    setnames(x, c(transcript_id, gene_id), 
                c("transcript_id", "gene_id"))
    ans = extract_feature(x, unique(x$feature), subfeature, 
                match.arg(type), ignore_strand, ...)
    new(class(ans)[1L], as(setDF(ans), "GRanges"))
}

# ----- internal functions for extracting specific features ----- #

extract_feature <- function(x, uniq_features, feature, type, 
                        ignore_strand=FALSE, ...) {
    UseMethod("extract_feature")
}

extract_feature.default <- function(x, ...) {
    stop("No default method implemented.")
}

extract_feature.gene <- function(x, uniq_features, feature, type, 
                        ignore_strand=FALSE, ...) {

    # to pleae R CMD CHECK
    seqnames=gene_id=transcript_id=NULL
    if (is.null(feature)) {
        if (!"exon" %in% uniq_features)
            stop("'exon' must be a feature in input object 'x'.")
        construct_transcripts <- function(x) {
            x[feature %chin% "exon", .(seqnames=seqnames[1L], start=min(start), 
                    end=max(end), strand=strand[1L], gene_id=gene_id[1L]), 
            by="transcript_id"]
        }
        transcripts = construct_transcripts(x)
    } else {
        cols=c("transcript_id", "seqnames", "start", "end", "strand", "gene_id")
        .feature = feature
        transcripts = x[feature %chin% .feature, cols, with=FALSE]
    }
    ans = switch (type, 
                    default = {
                        genes = transcripts[, .(seqnames=seqnames[1L], 
                            start=min(start), end=max(end), strand=strand[1L], 
                        transcript_id=paste(unique(transcript_id), 
                            collapse=";")), by="gene_id"]
                    },
                    union = {
                        genes = suppressWarnings(reduce_overlaps(transcripts, 
                                by="gene_id", ignore_strand=ignore_strand))
                        olaps = find_overlaps(genes, transcripts, 
                                    ignore_strand=FALSE)
                        olaps[, "transcript_id" := 
                                    transcripts$transcript_id[subjectHits]]
                        olaps = olaps[, .(transcript_id = 
                                paste(unique(transcript_id), collapse=";")), 
                                by="queryHits"]
                        genes[olaps$queryHits, 
                                "transcript_id" := olaps$transcript_id]
                    },
                    disjoin = {
                        genes = suppressWarnings(disjoin_overlaps(transcripts, 
                                by="gene_id", ignore_strand=ignore_strand))
                        olaps = find_overlaps(genes, transcripts, 
                                    ignore_strand=FALSE)
                        olaps[, "transcript_id" := 
                                    transcripts$transcript_id[subjectHits]]
                        olaps = olaps[, .(transcript_id = 
                                paste(unique(transcript_id), collapse=";")), 
                                by="queryHits"]
                        genes[olaps$queryHits, 
                                "transcript_id" := olaps$transcript_id]
                    }, 
                    intersect = {
                        genes = transcripts[, if (max(start) <= min(end)) 
                                .(seqnames=seqnames[1L], start=max(start), 
                                end=min(end), strand=strand[1L], 
                                transcript_id=paste(unique(transcript_id), 
                                    collapse=";")), by="gene_id"]
                    }, 
                    longest = {
                        max_len_idx = transcripts[, .I[which.max(end-start+1L
                            )], by="gene_id"]$V1
                        genes = transcripts[max_len_idx]
                    }, 
                    shortest = {
                        min_len_idx = transcripts[, .I[which.min(end-start+1L
                            )], by="gene_id"]$V1
                        genes = transcripts[min_len_idx]
                    }, 
                    overlap = {
                        stop("Not yet implemented.")
                    }                   
                )
    genes[, "length" := end-start+1L]
    olaps = find_overlaps(genes, genes, type="any", select="all", 
                ignore_strand=FALSE, ignore_redundant=FALSE)
    olaps = strictly_nonunique(setorder(olaps), "queryHits")
    olaps = olaps[, .(gene_id = paste(genes$gene_id[subjectHits], 
                collapse=";")), by="queryHits"]
    genes[, "overlaps":=gene_id][olaps$queryHits, "overlaps":=olaps$gene_id]

    colorder = c("seqnames", "start", "end", "length", "strand", 
                    "transcript_id", "gene_id", "overlaps")
    setcolorder(genes, colorder)
    genes[]
}

extract_feature.exon <- function(x, uniq_features, feature, type, 
                        ignore_strand=FALSE, ...) {
    stopifnot("exon" %in% uniq_features)

    # to please R CMD CHECK
    transcript_id=NULL
    cols = c("transcript_id", "seqnames", "start", "end", "strand", "gene_id")
    .feature = if (is.null(feature)) "exon" else feature
    exons = x[feature %chin% .feature, cols, with=FALSE]
    ans = switch (type, 
                    default = {
                        # exons are already extracted
                        ;
                    },
                    union = {
                        # gets all columns except 'transcript_id'
                        # warnings are from BiocGenerics..
                        exons_red = suppressWarnings(reduce_overlaps(exons, 
                                    by="gene_id", ignore_strand=ignore_strand))
                        olaps = find_overlaps(exons_red, exons, 
                                    ignore_strand=ignore_strand)
                        olaps[, "transcript_id" := 
                                    exons$transcript_id[subjectHits]]
                        olaps = olaps[, .(transcript_id = paste(unique(
                                transcript_id), collapse=";")), by="queryHits"]
                        exons_red[olaps$queryHits, 
                                "transcript_id" := olaps$transcript_id]
                        exons = exons_red
                    },
                    disjoin = {
                        exons_red = suppressWarnings(disjoin_overlaps(exons, 
                                    by="gene_id", ignore_strand=ignore_strand))
                        olaps = find_overlaps(exons_red, exons, 
                                    ignore_strand=ignore_strand)
                        olaps[, "transcript_id" := 
                                    exons$transcript_id[subjectHits]]
                        olaps = olaps[, .(transcript_id = paste(unique(
                                transcript_id), collapse=";")), by="queryHits"]
                        exons_red[olaps$queryHits, 
                                "transcript_id" := olaps$transcript_id]
                        exons = exons_red                       
                    }, 
                    intersect = {
                        # TODO: could we exclude the filtering step by 
                        # computing overlaps within 'gene_id'?
                        # should be straightforward with data. 
                        # table::foverlaps, as it can have multiple identifiers
                        # not sure how to do this with GRangesList object as I 
                        # don't get the same result.
                        # exons_isect = find_overlaps(exons, exons, 
                        #       ignore_strand=ignore_strand)
                        # exons_isect = exons_isect[exons$gene_id[queryHits] 
                        #       == exons$gene_id[subjectHits]]
                        stop("Not yet implemented.")
                    }, 
                    longest = {
                        stop("Not yet implemented.")                        
                    }, 
                    shortest = {
                        stop("Not yet implemented.")
                    }, 
                    overlap = {
                        stop("Not yet implemented.")
                    }                   
                )
    exons[, "length" := end-start+1L]
    # # TODO: mark overlapping genes
    # olaps = find_overlaps(ans, ans, type="any", select="all", 
    #                               ignore_strand=ignore_strand, 
    #                               ignore_redundant=FALSE)
    # olaps = strictly_nonunique(setorder(olaps), "queryHits")
    # olaps = olaps[, .(gene_id = paste(ans$gene_id[subjectHits], 
    #                       collapse=",")), by=queryHits]
    # ans[, "overlaps" := gene_id][olaps$queryHits, "overlaps" := 
    #                       olaps$gene_id]

    colorder = c("seqnames", "start", "end", "length", "strand", 
                    "transcript_id", "gene_id") #, "overlaps")
    setcolorder(exons, colorder)
    exons[]
}

extract_feature.intron <- function(x, uniq_features, feature, type, 
                            ignore_strand=FALSE, ...) {
    stopifnot("exon" %in% uniq_features)
    # to please R CMD CHECK
    transcript_id=NULL
    cols = c("transcript_id", "seqnames", "start", "end", "strand", "gene_id")
    .feature = if (is.null(feature)) "intron" else feature
    introns = x[feature %chin% .feature, cols, with=FALSE]
    if (nrow(introns) == 0L) {
        # most likely the object doesn't have introns, so construct them
        warning("No introns found in input object. Attempting to construct 
            introns using construct_introns().")
        xx = new(class(x)[2L], as_granges(x))
        introns = as_data_table(construct_introns(xx, update=FALSE))
        setattr(introns, 'class', data.table::copy(class(x)))
    }
    if (nrow(introns) == 0L)
        stop("construct_introns() resulted in 0 introns as well. 
                Nothing to do.")
    ans = switch (type, 
                    default = {
                        # introns are already extracted
                        ;
                    },
                    union = {
                        # gets all columns except 'transcript_id'
                        # warnings are from BiocGenerics..
                        introns_red = suppressWarnings(reduce_overlaps(introns
                                , by="gene_id", ignore_strand=ignore_strand))
                        olaps = find_overlaps(introns_red, introns, 
                                    ignore_strand=ignore_strand)
                        olaps[, "transcript_id" := 
                                    introns$transcript_id[subjectHits]]
                        olaps = olaps[, .(transcript_id = paste(unique(
                                transcript_id), collapse=";")), by="queryHits"]
                        introns_red[olaps$queryHits, 
                                "transcript_id" := olaps$transcript_id]
                        introns = introns_red
                    },
                    disjoin = {
                        introns_red = suppressWarnings(disjoin_overlaps(introns
                                , by="gene_id", ignore_strand=ignore_strand))
                        olaps = find_overlaps(introns_red, introns, 
                                ignore_strand=ignore_strand)
                        olaps[, "transcript_id" := 
                                introns$transcript_id[subjectHits]]
                        olaps = olaps[, .(transcript_id = paste(unique(
                                transcript_id), collapse=";")), by="queryHits"]
                        introns_red[olaps$queryHits, 
                                "transcript_id" := olaps$transcript_id]
                        introns = introns_red                       
                    }, 
                    intersect = {
                        # TODO: could we exclude the filtering step by 
                        # computing overlaps within 'gene_id'?
                        # should be straightforward with data.
                        # table::foverlaps, as it can have multiple identifiers
                        # not sure how to do this with GRangesList object as I 
                        # don't get the same result.
                        # introns_isect = find_overlaps(introns, introns, 
                        #                   ignore_strand=ignore_strand)
                        # introns_isect = introns_isect[introns$gene_id[
                        #         queryHits] == introns$gene_id[subjectHits]]
                        stop("Not yet implemented.")
                    }, 
                    longest = {
                        stop("Not yet implemented.")                        
                    }, 
                    shortest = {
                        stop("Not yet implemented.")
                    }, 
                    overlap = {
                        stop("Not yet implemented.")
                    }                   
                )
    introns[, "length" := end-start+1L]
    # # TODO: mark overlapping genes
    # olaps = find_overlaps(ans, ans, type="any", select="all", 
    #                               ignore_strand=ignore_strand, 
    #                               ignore_redundant=FALSE)
    # olaps = strictly_nonunique(setorder(olaps), "queryHits")
    # olaps = olaps[, .(gene_id = paste(ans$gene_id[subjectHits], 
    #                       collapse=",")), by=queryHits]
    # ans[, "overlaps" := gene_id][olaps$queryHits, "overlaps" := 
    #                               olaps$gene_id]

    colorder = c("seqnames", "start", "end", "length", "strand", 
                    "transcript_id", "gene_id") #, "overlaps")
    setcolorder(introns, colorder)
    introns[]
}

# ----- Helper functions for extract_features --------------------------------

is.gtf <- function(x) inherits(x, 'gtf')
is.gff <- function(x) inherits(x, 'gff')
is.bed <- function(x) inherits(x, 'bed')
is.bam <- function(x) inherits(x, 'bam')

#' @title Convert gtf/gff/bam/bed S4 objects to data.table and preserve class
#' 
#' @description \code{as.data.table} converts the input object to 
#' \code{data.table}, but strips away all other classes except 
#' \code{data.table} and \code{data.frame}. 
#' 
#' This internal function is a simple wrapper to \code{as.data.table} but 
#' also retains the extra class information of the input object.
#' 
#' @param x An object that inherits from \code{GRanges}.
#' @param reset_class logical (default \code{FALSE}). If \code{TRUE}, resets 
#' the class to \code{"data.table", "data.frame"}.
#' @aliases as_data_table
#' @return A data.table object with input class preserved.
#' @examples
#' \dontrun{
#' setClass("foo", contains="GRanges")
#' x = new("foo", GRanges(letters[1:5], IRanges(1:5, 6:10)))
#' class(x) # [1] "foo"
#' class(gread:::as_data_table(x)) # [1] "foo" "data.table" "data.frame"
#' }
as_data_table <- function(x, reset_class=FALSE) {
    stopifnot(inherits(x, "GRanges"))
    cl = if (identical(class(x), "GRanges") || reset_class) NULL else class(x)
    ans = as.data.table(x)
    setattr(ans, 'class', c(cl, "data.table", "data.frame"))
    ans[]
}
