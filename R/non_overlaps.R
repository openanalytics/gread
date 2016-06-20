#' @title Generate non-overlapping intron coordinates
#'
#' @description This function takes an \code{gtf/gff} object. If the 
#' features \code{exon} and \code{intron} are present, then, it attempts 
#' to find a set of exon and intron coordinates that do NOT overlap with each 
#' other. If the intron coordinates are not present, it constructs the introns 
#' first using \code{gread::construct_introns} and then attempts to extract 
#' non-overlapping exonic and intronic coordinates. If \code{exon} feature is 
#' not present, it returns an error.
#' 
#' Note: This function is not exported at the moment.
#' 
#' @param x An object of class \code{gtf/gff} object which has to have 
#' features \code{intron} and \code{exon}.
#' @param transcript_id Column name in \code{x} corresponding to transcript id.
#' @param gene_id Column name in \code{x} corresponding to gene id.
#' @return Returns a list of data.table containing \emph{non overlapping 
#' exons} and \emph{non overlapping introns}.
#' @examples
#' \dontrun{
#' path = system.file("tests", package="gread")
#' gtf_file = file.path(path, "sample.gtf")
#' gtf = read_format(gtf_file)
#' non_overlaps(gtf)
#' }
non_overlaps <- function(x, transcript_id="transcript_id", gene_id="gene_id"){

    # to please R CMD CHECK
    feature=idx=epos=ecnt=keep=edge=check=V1=seqname=pos=e.gid=i.gid=N=NULL
    # taken from gffutils::non_overlaps    
    stopifnot(is.gtf(x) || is.gff(x), "feature" %chin% names(x), 
        transcript_id %chin% names(x), gene_id %chin% names(x))
    uniq_features = unique(x[["feature"]])
    if (!"intron" %in% uniq_features) x = construct_introns(x)
    uniq_features = unique(x[["feature"]])
    stopifnot(all(c("intron", "exon") %in% uniq_features))
    i = x[feature == "intron"]
    e = x[feature == "exon"]
    e[, `:=`(ecnt = .N, epos = 1:.N), by=c(transcript_id)]

    i_gr = as_granges(i)

    # non-overlapping intron
    # "within" contains the potential IR events. So we shouldn't delete them
    a = find_overlaps(i, e, type="any", ignore_strand=TRUE)    
    w = find_overlaps(i, e, type="within", ignore_strand=TRUE)
    
    # idea: remove all `any` events that are not as well `within`. Also, 
    # equate the corresponding intronic and exonic gene_id and obtain a 
    # logical vector idx
    a_minus_w = a[!w, on=names(a)
                ][, `:=`(idx = i[[gene_id]][queryHits] == 
                                    e[[gene_id]][subjectHits], 
                        epos = e$epos[subjectHits], 
                        ecnt = e$ecnt[subjectHits])]
    setkey(a_minus_w, queryHits)

    # if the logical vector idx sums to `0`, when grouped by `queryHits` then, 
    # this means that this intron coordinate overlaps only with an overlapping 
    # gene and so we'll have to keep it. Hence the `keep != 0`
    a_minus_w = a_minus_w[, list(N = .N, keep = sum(idx), edge = 
                    (all(epos == 1) | all(epos == ecnt))), by=queryHits]
    a_minus_w = a_minus_w[keep != 0][, keep := NULL]

    # Now resolve end-cases and special "both start,end unequal overlap" cases 
    # scenario A: we'll remove this intron and it'll be falsely detected as a 
    # CI event... We'll have to retain them!
    # ==========-------------========= (iso1) (will be removed usually)
    #                   ============== (iso2) 
    a_minus_w = a_minus_w[!(edge)][, edge := NULL]
    
    # scenario B: we'll remove introns of this fashion:
    # ============--------------============= (iso1)
    # ===============----------------======== (iso2) (will be removed usually)
    # In this case we won't identify the isoforms associated with these events 
    # here. So, we'll have to recover them as well. How? by finding intron 
    # overlaps with itself. Can't use find_overlaps here as it requires x,y 
    # arguments.
    io = findOverlaps(i_gr, drop.self=TRUE, type="any", drop.redundant=FALSE)
    io = setDT(list(queryHits=queryHits(io), subjectHits=subjectHits(io)))
    setkey(io)
    
    # you'll have to subtract "equal" or ("exact") overlaps:
    io.eq = findOverlaps(i_gr, drop.self=TRUE, type="equal", 
                drop.redundant=FALSE)
    io.eq = setDT(list(queryHits=queryHits(io.eq), 
                        subjectHits=subjectHits(io.eq)))
    setkey(io.eq)
    io = io[!io.eq]
    setkey(io, queryHits)
    io[, check := i$start[queryHits] != i$start[subjectHits] & 
                    i$end[queryHits] != i$end[subjectHits]]
    io.fil = io[, sum(check) == .N, by=queryHits][V1 == TRUE]

    # all these io.fil must be retained. So, we've to remove them from 
    # a_minus_w.
    a_minus_w = a_minus_w[!list(io.fil$queryHits)] 
    # old code: a_minus_w[!queryHits %in% io.fil$queryHits]
    # The above lines will take care of scenario B.
    # changes for 7th August ends here

    # Now, we have in a_minus_w all the indices of introns that we'd like 
    # to remove as, they overlap with at least one exon of their own gene, 
    # meaning there is a smaller intron at the same position for this gene 
    # that we've already accounted for
    ni = i[!a_minus_w$queryHits]
    setkey(ni, "seqname", "start", "end")
    tmp__ = ni[[transcript_id]]
    ni[, c(transcript_id) := as.list(tmp__)]
    ni = ni[, c("seqname", "start", "end", "strand", 
                    transcript_id, gene_id), with=FALSE]
    count = ni[, .N, by = list(seqname, start, end)]
    # code here is moved to the end... July 25 (due to overlapping genes 
    # interpretation difference between intron and exon)

    # non_overlapping exon
    # edited May 30th 2013, redid the whole thing: July 17 2013
    # data.table ABSOLUTELY ROCKS!
    ne = e[, list(seqname=seqname[1L], pos=c(start, end)), by=c(gene_id)]
    setkey(ne)
    ne = ne[, list(seqname=seqname[1L], start=pos[1:(.N-1)], end=pos[2:.N]), 
                    by=c(gene_id)]

    nolaps = find_overlaps(ne, ni, type = "any")
    nolaps[, `:=`(e.gid=ne[[gene_id]][queryHits], 
                    i.gid=ni[[gene_id]][subjectHits])]
    nolaps = nolaps[e.gid == i.gid]

    ne = ne[!unique(nolaps$queryHits)]
    ne = ne[start != end]
    ne = ne[, list(seqname=seqname[1L], pos=c(start, end)), by=c(gene_id)]
    setkeyv(ne, c(gene_id, "pos"))
    ne = ne[ne[, .N, by=key(ne)][N == 1]][, N := NULL]
    ne = ne[, list(seqname=seqname[1L], start=pos[c(TRUE, FALSE)], 
                    end=pos[c(FALSE, TRUE)]), by=c(gene_id)]
    setkeyv(ne, c(gene_id))

    setnames(e, transcript_id, "transcript_id")
    on.exit(setnames(e, "transcript_id", transcript_id), add=TRUE)
    ne.rest = e[, list(strand=strand[1L], 
                    transcript_id=list(unique(transcript_id))), by=c(gene_id)]
    setnames(ne.rest, "transcript_id", transcript_id)
    ne = ne[ne.rest]
    setcolorder(ne, names(ni))
    setkey(ne, seqname, start, end)
    
    # now get overlapping intron
    ni[count[N > 1], `:=`(strand=paste(unique(strand), sep="", collapse=""), 
                        transcript_id=list(unlist(transcript_id)), 
                        gene_id=list(unique(unlist(gene_id)))), 
       by=.EACHI] ## version 1.9.3+
    setnames(ni, c("transcript_id", "gene_id"), c(transcript_id, gene_id))
    ni = unique(ni)
    # calling unique on 'ni' is more efficient than this, since 'ni' is key'd. 
    # changed July 17 
    # .gff[["non_overlapping_intron"]] = ni[count, mult = "first"][, N := NULL]
    setkey(ne, NULL)
    setkey(ni, NULL)
    list(non_overlapping_exon = ne, non_overlapping_intron = ni)
}
