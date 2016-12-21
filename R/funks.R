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
#'
#' @export
#' @return A character holding the GO terms for \code{gene.id}
compoundGoAnnotation <- function(gene.id, goa.tbl = getOption("MapMan2GO.goa.tbl", 
    ukb.goa.hits), gene.col = getOption("MapMan2GO.goa.tbl.gene.col", 3), go.col = getOption("MapMan2GO.goa.tbl.go.col", 
    2)) {
    unique(sort(goa.tbl[which(goa.tbl[, gene.col] == gene.id), go.col]))
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
#' 'genes.goa' is a list with each genes' compound GO Annotations, and
#' 'MapManBin.GO' is the intersection of the genes' compound GO Annotations to
#' be used as the MapMan-Bin's compound GO Annotation, n.GO is the number of GO
#' Terms in the MapMan-Bin's compound GO Annotation, median.n.GO is the median
#' of the number of GO Terms in the genes' GOAs, and n.genes is the number of
#' genes related to \code{map.man.bin}.
compoundGoAnnotationEntropy <- function(map.man.bin, mm.bins.vs.genes = getOption("MapMan2GO.seq.sim.tbl", 
    mm.bins.vs.sprot), mm.bin.col = getOption("MapMan2GO.seq.sim.tbl.bin.col", "MapManBin"), 
    mm.gene.col = getOption("MapMan2GO.seq.sim.tbl.gene.col", "Swissprot.Short.ID"), 
    goa.tbl = getOption("MapMan2GO.goa.tbl", ukb.goa.hits), gene.col = getOption("MapMan2GO.goa.tbl.gene.col", 
        3), go.col = getOption("MapMan2GO.goa.tbl.go.col", 2)) {
    gene.ids <- mm.bins.vs.genes[which(mm.bins.vs.genes[, mm.bin.col] == map.man.bin), 
        mm.gene.col]
    genes.goa <- setNames(lapply(gene.ids, function(g.id) {
        compoundGoAnnotation(g.id, goa.tbl, gene.col, go.col)
    }), gene.ids)
    s.e <- entropy(table(as.character(unlist(lapply(genes.goa, paste, collapse = ",")))))
    bin.goa <- sort(Reduce(intersect, genes.goa))
    list(Shannon.Entropy = s.e, genes.goa = genes.goa, MapManBin.GO = paste(bin.goa, 
        collapse = ","), n.GO = length(bin.goa), median.n.GO = median(unlist(lapply(genes.goa, 
        length)), na.rm = TRUE), n.genes = length(gene.ids))
}
