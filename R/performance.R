#' Compute the F1-Score for the prediction \code{predicted.annos} as the
#' harmonic mean of precision and recall.
#'
#' @param true.annos A vector of atomic annotations representing the TRUTH.
#' @param predicted.annos A vector of atomic annotations representing the
#' PREDICTIONS.
#'
#' @return Number - The weighted harmonic mean of precision and recall.
#' @export
f1Score <- function(true.annos, predicted.annos) {
    true.pos <- sum(true.annos %in% predicted.annos)
    
    if (identical(true.annos, predicted.annos)) {
        return(1)
    } else if (true.pos == 0) {
        return(0)
    }
    
    precision <- true.pos/length(predicted.annos)
    recall <- true.pos/length(true.annos)
    
    2 * precision * recall/(precision + recall)
}

#' Tests the function \code{f1Score}.
#'
#' @return \code{TRUE} if and only if all tests pass successfully.
#' @export
testF1Score <- function() {
    t.1 <- f1Score(c(), c()) == 1
    t.2 <- f1Score(c(), LETTERS[1:3]) == 0
    t.3 <- f1Score(LETTERS[1:3], LETTERS[1:3]) == 1
    t.4 <- f1Score(LETTERS[1:4], LETTERS[3:6]) == 0.5
    all(c(t.1, t.2, t.3, t.4))
}

#' For each row in \c{mercator.annos} measures the respective F1-Score of the
#' MapMan-Bin annotation based on the reference annotations.
#'
#' @param mercator.annos an instance of \c{base::data.frame} produced by
#' invoking function \c{MapMan2GO::readMercatorResultTable}.
#' @param reference.gene.annos An instance of \c{base::data.frame} holding the
#' reference annotations for the referene genes, i.e. the TRUTH the mercator
#' predictions are evaluated against. Default is
#' \c{getOption('MapMan2GO.reference.gene.annos', ukb.ref.goas)}.
#' @param rga.gene.col An integer or string indicating in which column of
#' \c{reference.gene.annos} to lookup the reference gene identifiers. Default
#' is \c{getOption('MapMan2GO.rga.gene.col', 'V4')}.
#' @param rga.anno.col An integer or string indicating in which column of
#' \c{reference.gene.annos} to lookup the reference gene annotations. Default
#' is \c{getOption('MapMan2GO.rga.anno.col', 'V2')}.
#'
#' @return A numeric vector holding the respective F1-Scores.
#' @export
mercatorF1Score <- function(mercator.annos, reference.gene.annos = getOption("MapMan2GO.reference.gene.annos", 
    ukb.ref.goas), rga.gene.col = getOption("MapMan2GO.rga.gene.col", "V4"), rga.anno.col = getOption("MapMan2GO.rga.anno.col", 
    "V2")) {
    unlist(mclapply(1:nrow(mercator.annos), function(i) {
        x <- mercator.annos[i, ]
        f1 <- NA
        if (!is.na(x$TYPE) && x$TYPE) {
            rga.i <- which(reference.gene.annos[, rga.gene.col] == toupper(x$IDENTIFIER))
            if (length(rga.i) > 0) {
                predicted.annos <- strsplit(x$MapManBin.GO, ",")[[1]]
                reference.annos <- addAncestors(reference.gene.annos[rga.i, rga.anno.col])
                f1 <- f1Score(reference.annos, predicted.annos)
            }
        }
        f1
    }))
}
