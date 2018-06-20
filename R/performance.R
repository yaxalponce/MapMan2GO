#' Compute the F1-Score for the prediction \code{predicted.annos} as the
#' harmonic mean of precision and recall.
#'
#' @param true.annos A vector of atomic annotations representing the TRUTH.
#' @param predicted.annos A vector of atomic annotations representing the
#' PREDICTIONS.
#'
#' @return An instance of \code{base::data.frame} with the following numeric
#' columns: precision, recall, true.pos, false.pos, f1.score
#' @export
f1Score <- function(true.annos, predicted.annos) {
    true.pos <- sum(true.annos %in% predicted.annos)
    false.pos <- length(predicted.annos) - true.pos
    
    if (identical(true.annos, predicted.annos)) {
        return(data.frame(precision = 1, recall = 1, true.pos = true.pos, false.pos = 0, 
            f1.score = 1, stringsAsFactors = FALSE))
    } else if (true.pos == 0) {
        return(data.frame(precision = 0, recall = 0, true.pos = true.pos, false.pos = 0, 
            f1.score = 0, stringsAsFactors = FALSE))
    }
    
    precision <- true.pos/length(predicted.annos)
    recall <- true.pos/length(true.annos)
    f1.score <- 2 * precision * recall/(precision + recall)
    
    data.frame(precision = precision, recall = recall, true.pos = true.pos, false.pos = false.pos, 
        f1.score = f1.score, stringsAsFactors = FALSE)
}

#' Tests the function \code{f1Score}.
#'
#' @return \code{TRUE} if and only if all tests pass successfully.
#' @export
testF1Score <- function() {
    t.1 <- f1Score(c(), c())$f1.score == 1
    t.2 <- f1Score(c(), LETTERS[1:3])$f1.score == 0
    t.3 <- f1Score(LETTERS[1:3], LETTERS[1:3])$f1.score == 1
    t.4 <- f1Score(LETTERS[1:4], LETTERS[3:6])$f1.score == 0.5
    all(c(t.1, t.2, t.3, t.4))
}

#' For each row in \code{mercator.annos} measures the respective F1-Score of the
#' MapMan-Bin annotation based on the reference annotations.
#'
#' @param mercator.annos an instance of \code{base::data.frame} produced by
#' invoking function \code{MapMan2GO::readMercatorResultTable}.
#' @param reference.gene.annos An instance of \code{base::data.frame} holding the
#' reference annotations for the referene genes, i.e. the TRUTH the mercator
#' predictions are evaluated against. Default is
#' \code{getOption('MapMan2GO.reference.gene.annos', ukb.ref.goas)}.
#' @param rga.gene.col An integer or string indicating in which column of
#' \code{reference.gene.annos} to lookup the reference gene identifiers. Default
#' is \code{getOption('MapMan2GO.rga.gene.col', 'V5')}.
#' @param rga.anno.col An integer or string indicating in which column of
#' \code{reference.gene.annos} to lookup the reference gene annotations. Default
#' is \code{getOption('MapMan2GO.rga.anno.col', 'V2')}.
#' @param process.ref.annos.funk Function applied to each reference gene's
#' annotation vector. Set to \code{base::identity} if you whish to leave them
#' unchanged. Because we are dealing with Gene Ontology Term annotions the
#' default here is \code{getOption( 'MapMan2GO.process.ref.annos.funk',
#' MapMan2GO::addAncestors )}.
#'
#' @return An instance of \code{base::data.frame} with the following numeric
#' columns: precision, recall, true.pos, false.pos, f1.score
#' @export
mercatorF1Score <- function(mercator.annos, reference.gene.annos = getOption("MapMan2GO.reference.gene.annos", 
    ukb.ref.goas), rga.gene.col = getOption("MapMan2GO.rga.gene.col", "V5"), rga.anno.col = getOption("MapMan2GO.rga.anno.col", 
    "V2"), process.ref.annos.funk = getOption("MapMan2GO.process.ref.annos.funk", 
    MapMan2GO::addAncestors)) {
    Reduce(rbind, mclapply(1:nrow(mercator.annos), function(i) {
        x <- mercator.annos[i, ]
        f1 <- data.frame(precision = NA, recall = NA, true.pos = NA, false.pos = NA, 
            f1.score = NA, stringsAsFactors = FALSE)
        if (!is.na(x$TYPE) && x$TYPE) {
            rga.i <- which(reference.gene.annos[, rga.gene.col] == x$IDENTIFIER)
            if (length(rga.i) > 0) {
                predicted.annos <- strsplit(x$MapManBin.GO, ",")[[1]]
                reference.annos <- process.ref.annos.funk(reference.gene.annos[rga.i, 
                  rga.anno.col])
                f1 <- f1Score(reference.annos, predicted.annos)
            }
        }
        f1
    }))
}

#' Computes the F1-scores of InterProScan protein function predictions.
#'
#' @param ipr.annos An instance of \code{base::data.frame} with at least the
#' following required columns: 'V2' holds the function annotions, and 'V4'
#' holds the gene identifiers.
#' @param reference.genes A character vector of gene identifiers forming the
#' complete set of reference genes ('gold standard').
#' @param reference.gene.annos An instance of \code{base::data.frame} holding the
#' reference annotations for the referene genes, i.e. the TRUTH the mercator
#' predictions are evaluated against. Default is
#' \code{getOption('MapMan2GO.reference.gene.annos', ukb.ref.goas)}.
#' @param rga.gene.col An integer or string indicating in which column of
#' \code{reference.gene.annos} to lookup the reference gene identifiers. Default
#' is \code{getOption('MapMan2GO.rga.gene.col', 'V4')}.
#' @param rga.anno.col An integer or string indicating in which column of
#' \code{reference.gene.annos} to lookup the reference gene annotations. Default
#' is \code{getOption('MapMan2GO.rga.anno.col', 'V2')}.
#' @param process.ref.annos.funk Function applied to each reference gene's
#' annotation vector. Set to \code{base::identity} if you whish to leave them
#' unchanged. Because we are dealing with Gene Ontology Term annotions the
#' default here is \code{getOption( 'MapMan2GO.process.ref.annos.funk',
#' MapMan2GO::addAncestors )}.
#'
#' @return An instance of \code{base::data.frame} with the following numeric
#' columns: precision, recall, true.pos, false.pos, f1.score, gene.
#' @export
interProScanF1Scores <- function(ipr.annos, reference.genes = getOption("MapMan2GO.reference.genes", 
    ref.gene.ids), reference.gene.annos = getOption("MapMan2GO.reference.gene.annos", 
    ukb.ref.goas), rga.gene.col = getOption("MapMan2GO.rga.gene.col", "V4"), rga.anno.col = getOption("MapMan2GO.rga.anno.col", 
    "V2"), process.ref.annos.funk = getOption("MapMan2GO.process.ref.annos.funk", 
    MapMan2GO::addAncestors)) {
    Reduce(rbind, mclapply(reference.genes, function(gene) {
        ref.annos <- process.ref.annos.funk(unique(reference.gene.annos[which(reference.gene.annos[, 
            rga.gene.col] == gene), rga.anno.col]))
        pred.annos <- unique(ipr.annos[which(ipr.annos$V4 == gene), "V2"])
        ipr.f1.df <- f1Score(ref.annos, pred.annos)
        ipr.f1.df$gene <- gene
        ipr.f1.df
    }))
}
