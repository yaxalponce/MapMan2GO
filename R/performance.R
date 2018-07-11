#' Compute the Matthew's correlation coefficient (mcc) and F1-Score for the
#' prediction \code{predicted.annos}. F1-Score as the harmonic mean of
#' precision and recall. mcc using \code{mccr::mccr}.
#'
#' @param true.annos A vector of atomic annotations representing the TRUTH.
#' @param predicted.annos A vector of atomic annotations representing the
#' PREDICTIONS.
#' @param universe.annos A vector of atomic annotations representing all
#' possible annotations that could theoretically be assigned, i.e. the
#' 'annotation universe'. Default is
#' \code{getOption('MapMan2GO.performance.universe.annotations',
#' ukb.ref.universe.gos)}.
#' @param na.for.empty.references boolean indicating whether to return all NA
#' values in case no references are given. In the case of protein function
#' prediction missing references can be interpreted as missing knowledge,
#' because there are no proteins with no function. Default is
#' \code{getOption('MapMan2GO.f1.na.for.empty.references', TRUE)}.
#'
#' @return An instance of \code{base::data.frame} with the following numeric
#' columns: n.truth, n.pred, true.pos, false.pos, precision, recall,
#' false.pos.rate, f1.score, mcc
#' @export
performanceScores <- function(true.annos, predicted.annos, universe.annos = getOption("MapMan2GO.performance.universe.annotations", 
    ukb.ref.universe.gos), na.for.empty.references = getOption("MapMan2GO.f1.na.for.empty.references", 
    TRUE)) {
    t.a <- unique(true.annos)
    p.a <- unique(predicted.annos)
    if ((is.null(t.a) || all(is.na(t.a)) || length(t.a) == 0) && na.for.empty.references) {
        return(data.frame(n.truth = 0, n.pred = length(p.a), true.pos = NA, false.pos = NA, 
            precision = NA, recall = NA, false.pos.rate = NA, f1.score = NA, mcc = NA, 
            stringsAsFactors = FALSE))
    }
    true.pos <- sum(t.a %in% p.a)
    false.pos <- length(p.a) - true.pos
    negatives <- setdiff(universe.annos, t.a)
    fp.rate <- if (length(negatives) == 0) {
        0
    } else {
        false.pos/length(negatives)
    }
    mcc <- mccr::mccr(universe.annos %in% t.a, universe.annos %in% p.a)
    
    precision <- if (length(p.a) == 0) 
        0 else true.pos/length(p.a)
    recall <- if (length(t.a) == 0) 
        0 else true.pos/length(t.a)
    f1.score <- if (precision + recall == 0) 
        0 else 2 * precision * recall/(precision + recall)
    
    data.frame(n.truth = length(t.a), n.pred = length(p.a), true.pos = true.pos, 
        false.pos = false.pos, precision = precision, recall = recall, false.pos.rate = fp.rate, 
        f1.score = f1.score, mcc = mcc, stringsAsFactors = FALSE)
}

#' Tests the function \code{performanceScores}.
#'
#' @return \code{TRUE} if and only if all tests pass successfully.
#' @export
testPerformanceScores <- function() {
    t.1 <- is.na(performanceScores(c(), c(), c(), TRUE)$f1.score)
    t.11 <- performanceScores(c(), c(), c(), FALSE)$f1.score == 0
    t.2 <- performanceScores(c(), LETTERS[1:3], LETTERS[1:3], FALSE)$f1.score == 
        0
    t.3 <- performanceScores(LETTERS[1:3], LETTERS[1:3], c(), FALSE)$f1.score == 
        1
    t.4 <- performanceScores(LETTERS[1:4], LETTERS[3:6], c(), FALSE)$f1.score == 
        0.5
    t.5 <- is.na(performanceScores(character(0), LETTERS[1:3], c(), TRUE)$f1.score)
    t.6 <- performanceScores(c(), c(), c(), FALSE)$mcc == 0
    t.7 <- performanceScores(c(), LETTERS[1:3], c(), FALSE)$mcc == 0
    t.8 <- performanceScores(LETTERS[1:3], LETTERS[1:3], LETTERS[1:3], FALSE)$mcc == 
        0
    t.9 <- performanceScores(LETTERS[1:3], LETTERS[4:6], LETTERS[1:6], FALSE)$mcc == 
        -1
    t.10 <- is.na(performanceScores(character(0), LETTERS[1:3], c(), TRUE)$mcc)
    all(c(t.1, t.2, t.3, t.4, t.5, t.6, t.7, t.8, t.9, t.10, t.11))
}


#' Computes the performance scores of protein function predictions.
#'
#' @param pred.annos An instance of \code{base::data.frame} with at least two
#' columns, one holding the protein identifier and the other the predicted
#' annotations.
#' @param pa.gene.col An integer or string indicating in which column of
#' \code{pred.annos} to lookup the gene identifiers.
#' @param pa.anno.col An integer or string indicating in which column of
#' \code{pred.annos} to lookup the predicted annotations.
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
#' @param process.predicted.annos.funk Function applied to the predicted
#' annotations before being passed to argument \code{process.annos.funk}. Set
#' to \code{base::identity} (the default), if you want to leave them unchanged.
#' Set to \code{MapMan2GO::splitMapManBinGOAs} if the function prediction
#' method was Mercator / MapMan2GO. Default is
#' \code{getOption('MapMan2GO.process.predicted.annos.funk', identity)}.
#' @param process.annos.funk Function applied to each
#' annotation vector. Set to \code{base::identity} if you whish to leave them
#' unchanged. Because we are dealing with Gene Ontology Term annotions the
#' default here is \code{getOption( 'MapMan2GO.process.annos.funk',
#' MapMan2GO::addAncestors )}.
#'
#' @return An instance of \code{base::data.frame} with the columns as returned
#' by \code{MapMan2GO::performanceScores} and the additional column \code{gene}
#' (the reference ientifier).
#' @export
predictorPerformance <- function(pred.annos, pa.gene.col, pa.anno.col, reference.genes = getOption("MapMan2GO.reference.genes", 
    ref.gene.ids), reference.gene.annos = getOption("MapMan2GO.reference.gene.annos", 
    ukb.ref.goas), rga.gene.col = getOption("MapMan2GO.rga.gene.col", "V4"), rga.anno.col = getOption("MapMan2GO.rga.anno.col", 
    "V2"), process.predicted.annos.funk = getOption("MapMan2GO.process.predicted.annos.funk", 
    identity), process.annos.funk = getOption("MapMan2GO.process.annos.funk", MapMan2GO::addAncestors)) {
    Reduce(rbind, mclapply(reference.genes, function(gene) {
        ref.annos <- process.annos.funk(reference.gene.annos[which(reference.gene.annos[, 
            rga.gene.col] == gene), rga.anno.col])
        pred.annos <- process.annos.funk(process.predicted.annos.funk(pred.annos[which(pred.annos[, 
            pa.gene.col] == gene), pa.anno.col]))
        pred.performance.df <- performanceScores(ref.annos, pred.annos)
        pred.performance.df$gene <- gene
        pred.performance.df
    }))
}



#' Function to parse a Blastp output table, retain only the best non self Hits,
#' where query and hit are identical, also discarding any hits that match
#' protein identifiers given one per line in argument
#' \code{exclude.accessions.lines.path}. This function expects protein
#' identifiers to be in the UniprotKB format, eg. 'sp|Q54IF9|MYBG_DICDI'. 
#'
#' @param path.2.blast.tbl Valid file path to the Blastp output table.
#' @param exclude.accessions.lines.path Valid file path to the file in which
#' protein identifiers are stored that are to be excluded from the hits.
#' Expected format is the short protein ID from uniprot, eg. 'Q54IF9'.
#'
#' @return An instance of \code{data.frame} with 14 columns, two additional to
#' the original blast tabular output: \code{'query.ukb.short.id'} and
#' \code{'hit.ukb.short.id'}, which hold the short Uniprot protein accessions
#' of the Query and Hit proteins respectively, eg. 'Q54IF9'. The other 12
#' columns are obtained from the Blast tabular output.
#' @export
extractBestBlastHits <- function(path.2.blast.tbl, exclude.accessions.lines.path) {
    blast.tbl <- read.table(path.2.blast.tbl, sep = "\t", stringsAsFactors = F)
    excl.ids <- readLines(exclude.accessions.lines.path)
    blast.tbl$query.ukb.short.id <- sub("^[^|]+\\|", "", sub("\\|[^|]+$", "", blast.tbl$V1))
    blast.tbl$hit.ukb.short.id <- sub("^[^|]+\\|", "", sub("\\|[^|]+$", "", blast.tbl$V2))
    blast.tbl.filtered <- blast.tbl[which(blast.tbl$V1 != blast.tbl$V2 & !blast.tbl$hit.ukb.short.id %in% 
        excl.ids), ]
    blast.tbl.sort <- blast.tbl.filtered[order(blast.tbl.filtered$V1, blast.tbl.filtered$V11), 
        ]
    blast.tbl.sort[!duplicated(blast.tbl.sort$V1), ]
}

#' Function to generate the annotations (GO Term predictions) for Best Blast.
#'
#' @param best.blast.tbl The result of invoking function
#' \code{MapMan2GO::extractBestBlastHits}.
#' @param blast.hit.goa The result of reading in a pre-processed form of the
#' UniProt GOA tables. See Vignette for more details. Is expected to have at
#' least two columns, 'V3' holding the Blast Hit identifier matching
#' 'query.ukb.short.id' in \code{best.blast.tbl} and 'V2' holding the Gene
#' Ontology Terms assigned to these respective Hits.
#'
#' @return An instance of \code{base::data.frame} with two columns: 'query'
#' holds the respective query proteins and 'GO' the GO predictions obtained by
#' the Best Blast method.
#' @export
bestBlastPredictions <- function(best.blast.tbl, blast.hit.goa) {
    Reduce(rbind, mclapply(1:nrow(best.blast.tbl), function(i) {
        query.id <- best.blast.tbl[i, "query.ukb.short.id"]
        hit.id <- best.blast.tbl[i, "hit.ukb.short.id"]
        goa.i <- which(blast.hit.goa$V3 == hit.id)
        if (length(goa.i) > 0) {
            gos <- blast.hit.goa[goa.i, "V2"]
            data.frame(query = query.id, GO = gos, stringsAsFactors = FALSE)
        }
    }))
}
