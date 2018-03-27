#' Extracts short UniProtKB gene accessions from the long version identifiers,
#' e.g. 'sp|Q8GWW7|AGUA_ARATH' will become 'Q8GWW7'.
#'
#' @param ukb.accs A character vector of long UniProtKB gene accessions
#'
#' @export
#' @return The short gene accessions extracted from the argument
#' \code{ukb.accs}.
sanitizeAccession <- function(ukb.accs) {
    sub("^[^|]+\\|", "", sub("\\|[^|]+$", "", ukb.accs))
}

#' Extracts the unique and sorted Gene Ontology Term Annotations for a given
#' gene.
#'
#' @param gene.id The unique identifier of the gene to lookup GOA for
#' @param goa.tbl An instance of \code{data.frame} holding GOAs. Default is
#' \code{getOption('MapMan2GO.goa.tbl', ukb.goa.hits)}
#' @param gene.col The column of \code{goa.tbl} number or name in which to
#' lookup the genes. Default is \code{getOption('MapMan2GO.goa.tbl.gene.col',
#' 3)}.
#' @param go.col The column of \code{goa.tbl} number or name in which to lookup
#' the GO terms. Default is \code{getOption('MapMan2GO.goa.tbl.go.col',2)}
#' @param extend.goas.with.ancestors boolean indicating whether to extend each
#' proteins' GOA with the ancestors of the respective GO Terms. Default is
#' \code{TRUE}.
#'
#' @export
#' @return A character holding the GO terms for \code{gene.id}
compoundGoAnnotation <- function(gene.id, goa.tbl = getOption("MapMan2GO.goa.tbl", 
    ukb.goa.hits), gene.col = getOption("MapMan2GO.goa.tbl.gene.col", 3), 
    go.col = getOption("MapMan2GO.goa.tbl.go.col", 2), extend.goas.with.ancestors = TRUE) {
    res.goa <- unique(sort(goa.tbl[which(goa.tbl[, gene.col] == gene.id), 
        go.col]))
    if (extend.goas.with.ancestors) {
        addAncestors(res.goa)
    } else {
        res.goa
    }
}

#' Adds all GO Terms that are ancestral to the argument \code{go.terms}.
#'
#' @param go.terms a character vector of GO identifier, e.g.
#' \code{'GO:006969'}.
#' @param ancestors a list of keys ontologies and values the respective
#' ancestral GO terms. Default is \code{list(MF = GOMFANCESTOR, BP =
#' GOBPANCESTOR, CC = GOCCANCESTOR)}.
#' @param exclude.root boolean indicating wether to exclude the ROOT GO Term
#' \code{'all'}. Default is \code{TRUE}.
#' @param root.go String indicating the ROOT GO Term. Default is \code{'all'}.
#'
#' @return A character vector including the argument \code{go.terms} and all
#' found ancestors.
#' @export
addAncestors <- function(go.terms, ancestors = list(MF = GOMFANCESTOR, 
    BP = GOBPANCESTOR, CC = GOCCANCESTOR), exclude.root = TRUE, root.go = "all") {
    sort(unique(unlist(lapply(go.term, function(g.id) {
        g.t <- GOTERM[[g.id]]
        if (!is.null(g.t) && length(g.t) > 0) {
            g.name <- attr(g.t, "Term")
            g.ont <- attr(g.t, "Ontology")
            g.incl.anc <- c(g.id, get(g.id, go.term.dist.to.root.funks[[g.ont]]))
            if (exclude.root) 
                g.incl.anc <- setdiff(g.incl.anc, root.go)
            g.incl.anc
        } else {
            g.id
        }
    }))))
}

#' Infers the compound GO annotations found for the genes related to the
#' argument MapMan Bin and measures the resulting Shannon Entropy of these
#' compound GO annotations.
#'
#' @param map.man.bin The identifier of the MapMan Bin to assign GO terms to.
#' @param mm.bins.vs.genes An instance of \code{data.frame} holding
#' MapManBin-Gene-Relations. Default is
#' \code{getOption('MapMan2GO.seq.sim.tbl', mm.bins.vs.sprot)}.
#' @param mm.bin.col The column of \code{mm.bins.vs.genes} in which to lookup
#' the MapMan-Bins. Default is \code{getOption('MapMan2GO.seq.sim.tbl.bin.col',
#' 'MapManBin')}.
#' @param mm.gene.col The column of \code{mm.bins.vs.genes} in which to lookup
#' the gene identifiers. Default is
#' \code{getOption('MapMan2GO.seq.sim.tbl.gene.col', 'Swissprot.Short.ID')}.
#' @param mm.bins.vs.genes An instance of \code{data.frame} with at least two
#' columns. It must hold mappings of \code{map.man.bin} to genes at least
#' partially found in \code{goa.tbl}.
#' @param goa.tbl An instance of \code{data.frame} holding GOAs. Default is
#' \code{getOption('MapMan2GO.goa.tbl', ukb.goa.hits)}
#' @param gene.col The column of \code{goa.tbl} number or name in which to
#' lookup the genes. Default is \code{getOption('MapMan2GO.goa.tbl.gene.col',
#' 3)}.
#' @param go.col The column of \code{goa.tbl} number or name in which to lookup
#' the GO terms. Default is \code{getOption('MapMan2GO.goa.tbl.go.col',2)}
#'
#' @export
#' @return An instance of \code{list} with the following named entries:
#' 'Shannon.Entropy' is the measured entropy of compound GO annotations
#' retrieved for the genes related to \code{map.man.bin}. The second entry
#' 'genes.goa' is a list with each genes' compound GO Annotations,
#' 'mutual.information' holds the mutual information between the Bin's GOA and
#' those found for the Bin's reference proteins, 'MapManBin.GO' is the
#' intersection of the genes' compound GO Annotations to be used as the
#' MapMan-Bin's compound GO Annotation, n.GO is the number of GO Terms in the
#' MapMan-Bin's compound GO Annotation, median.n.GO is the median of the number
#' of GO Terms in the genes' GOAs, and n.genes is the number of genes related
#' to \code{map.man.bin}.
compoundGoAnnotationEntropy <- function(map.man.bin, mm.bins.vs.genes = getOption("MapMan2GO.seq.sim.tbl", 
    mm.bins.vs.sprot), mm.bin.col = getOption("MapMan2GO.seq.sim.tbl.bin.col", 
    "MapManBin"), mm.gene.col = getOption("MapMan2GO.seq.sim.tbl.gene.col", 
    "Swissprot.Short.ID"), goa.tbl = getOption("MapMan2GO.goa.tbl", ukb.goa.hits), 
    gene.col = getOption("MapMan2GO.goa.tbl.gene.col", 3), go.col = getOption("MapMan2GO.goa.tbl.go.col", 
        2)) {
    gene.ids <- mm.bins.vs.genes[which(mm.bins.vs.genes[, mm.bin.col] == 
        map.man.bin), mm.gene.col]
    genes.goa <- setNames(lapply(gene.ids, function(g.id) {
        compoundGoAnnotation(g.id, goa.tbl, gene.col, go.col)
    }), gene.ids)
    s.e <- entropy(table(as.character(unlist(lapply(genes.goa, paste, collapse = ",")))))
    m.i <- mutualInformationBinGoaGenesGoas(genes.goa)
    bin.goa <- sort(Reduce(intersect, genes.goa))
    list(Shannon.Entropy = s.e, genes.goa = genes.goa, mutual.information = m.i, 
        MapManBin.GO = paste(bin.goa, collapse = ","), n.GO = length(bin.goa), 
        median.n.GO = median(unlist(lapply(genes.goa, length)), na.rm = TRUE), 
        n.genes = length(gene.ids))
}

#' Computes the empirical counts of events in a sample. Events are GO Term
#' annotations.
#'
#' @param go.sample.space A vector representing the sample space.
#' @param go.annos A vector of samples
#'
#' @export
#' @return A named numeric of empirical counts, names are the events in
#' \code{go.sample.space} and values are the empirical counts as found in
#' \code{go.annos}.
goCounts <- function(go.sample.space, go.annos) {
    setNames(as.numeric(lapply(go.sample.space, function(go.t) length(which(go.annos == 
        go.t)))), go.sample.space)
}

#' Compute the mutual information between a MapMan Bin's GOA and the GOAs found
#' for it's reference genes. In this, regard each single GO Term as a event in
#' the sample space. Because each GO Term can only have binary count values in
#' the Bin's GOA the Bin's counts are normalized by multiplying with them with
#' the number of reference genes.
#'
#' @param genes.goa A list of character vectors. Names are reference gene
#' identifiers and values, the character vectors, represent each gene's GO Term
#' Annotations. See \code{compoundGoAnnotationEntropy} for more details.
#'
#' @export
#' @return A numeric value, the computed mutual information in bits (log2).
mutualInformationBinGoaGenesGoas <- function(genes.goa) {
    genes.gos <- Reduce(c, genes.goa)
    bin.gos <- Reduce(intersect, genes.goa)
    go.sample.space <- unique(genes.gos)
    genes.go.counts <- goCounts(go.sample.space, genes.gos)
    bin.go.counts.norm <- length(genes.goa) * goCounts(go.sample.space, 
        bin.gos)
    mi.plugin(rbind(genes.go.counts, bin.go.counts.norm), unit = "log2")
}

#' Generates a two row plot with the first one being a Histogram and the secons
#' row a horizontal Boxplot.
#'
#' @param x The values passed into \code{hist} and \code{boxplot}
#' @param The main title of the resulting plot
#'
#' @export
#' @return TRUE if and only if no error has occurred
plotDistAsHistAndBox <- function(x, main) {
    def.mar <- par("mar")
    m.1 <- def.mar
    m.1[[1]] <- 0
    op <- par(mfcol = 2:1, mar = m.1)
    hist(x, col = "lightgrey", main = main, xlab = NULL)
    m.2 <- def.mar
    m.2[[3]] <- 0
    par(mar = m.2)
    boxplot(x, col = "lightgrey", horizontal = TRUE, frame = FALSE, pch = "|")
    par(op)
    TRUE
}
