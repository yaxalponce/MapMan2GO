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
#' columns: precision, recall, false.pos.rate, mcc, f1.score
#' @export
f1Score <- function(true.annos, predicted.annos, universe.annos = getOption("MapMan2GO.performance.universe.annotations", 
    ukb.ref.universe.gos), na.for.empty.references = getOption("MapMan2GO.f1.na.for.empty.references", 
    TRUE)) {
    t.a <- unique(true.annos)
    p.a <- unique(predicted.annos)
    if ((is.null(t.a) || all(is.na(t.a)) || length(t.a) == 0) && na.for.empty.references) {
        return(data.frame(precision = NA, recall = NA, false.pos.rate = NA, f1.score = NA, 
            mcc = NA, stringsAsFactors = FALSE))
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
    
    # if (identical(t.a, p.a)) { return(data.frame(precision = 1, recall = 1,
    # false.pos.rate = 0, f1.score = 1, mcc=1, stringsAsFactors = FALSE)) } else
    # if (true.pos == 0) { return(data.frame(precision = 0, recall = 0,
    # false.pos.rate = fp.rate, f1.score = 0, mcc = -1, stringsAsFactors = FALSE))
    # }
    
    precision <- true.pos/length(p.a)
    recall <- true.pos/length(t.a)
    f1.score <- 2 * precision * recall/(precision + recall)
    
    data.frame(precision = precision, recall = recall, false.pos.rate = fp.rate, 
        f1.score = f1.score, mcc = mcc, stringsAsFactors = FALSE)
}

#' Tests the function \code{f1Score}.
#'
#' @return \code{TRUE} if and only if all tests pass successfully.
#' @export
testF1Score <- function() {
    t.1 <- is.na(f1Score(c(), c(), c(), FALSE)$f1.score)
    t.2 <- is.na(f1Score(c(), LETTERS[1:3], c(), FALSE)$f1.score)
    t.3 <- f1Score(LETTERS[1:3], LETTERS[1:3], c(), FALSE)$f1.score == 1
    t.4 <- f1Score(LETTERS[1:4], LETTERS[3:6], c(), FALSE)$f1.score == 0.5
    t.5 <- is.na(f1Score(character(0), LETTERS[1:3], c(), TRUE)$f1.score)
    t.6 <- f1Score(c(), c(), c(), FALSE)$mcc == 0
    t.7 <- f1Score(c(), LETTERS[1:3], c(), FALSE)$mcc == 0
    t.8 <- f1Score(LETTERS[1:3], LETTERS[1:3], LETTERS[1:3], FALSE)$mcc == 0
    t.9 <- f1Score(LETTERS[1:3], LETTERS[4:6], LETTERS[1:6], FALSE)$mcc == -1
    t.10 <- is.na(f1Score(character(0), LETTERS[1:3], c(), TRUE)$mcc)
    all(c(t.1, t.2, t.3, t.4, t.5, t.6, t.7, t.8, t.9, t.10))
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
#' @return An instance of \code{base::data.frame} with the folque gold standard
#' proteinslowing numeric columns: precision, recall, true.pos, false.pos,
#' f1.score, mcc, gene (the reference)
#' @export
mercatorF1Score <- function(mercator.annos, reference.gene.annos = getOption("MapMan2GO.reference.gene.annos", 
    ukb.ref.goas), rga.gene.col = getOption("MapMan2GO.rga.gene.col", "V5"), rga.anno.col = getOption("MapMan2GO.rga.anno.col", 
    "V2"), process.ref.annos.funk = getOption("MapMan2GO.process.ref.annos.funk", 
    MapMan2GO::addAncestors)) {
    # Iterate over unique gold standard proteins
    unique.reference.prot.ids <- unique(reference.gene.annos$V5)
    Reduce(rbind, mclapply(1:length(unique.reference.prot.ids), function(i) {
        x <- unique.reference.prot.ids[i]
        f1 <- data.frame(precision = NA, recall = NA, true.pos = NA, false.pos = NA, 
            f1.score = NA, stringsAsFactors = FALSE)
        rga.i <- which(mercator.annos$IDENTIFIER == x)
        if (length(rga.i) > 0) {
            predicted.annos <- unique(unlist(strsplit(mercator.annos$MapManBin.GO[rga.i], 
                ",")))
            ref.ann.coor <- which(reference.gene.annos[, rga.gene.col] == x)
            reference.annos <- process.ref.annos.funk(reference.gene.annos[ref.ann.coor, 
                rga.anno.col])
            f1 <- f1Score(reference.annos, predicted.annos)
        } else {
            predicted.annos <- character(0)
            f1 <- f1Score(reference.annos, predicted.annos)
        }
        f1$gene <- x
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
#' columns: precision, recall, true.pos, false.pos, f1.score, mcc, gene.
#' @export
interProScanF1Scores <- function(ipr.annos, reference.genes = getOption("MapMan2GO.reference.genes", 
    ref.gene.ids), reference.gene.annos = getOption("MapMan2GO.reference.gene.annos", 
    ukb.ref.goas), rga.gene.col = getOption("MapMan2GO.rga.gene.col", "V4"), rga.anno.col = getOption("MapMan2GO.rga.anno.col", 
    "V2"), process.ref.annos.funk = getOption("MapMan2GO.process.ref.annos.funk", 
    MapMan2GO::addAncestors)) {
    Reduce(rbind, mclapply(reference.genes, function(gene) {
        ref.annos <- process.ref.annos.funk(unique(reference.gene.annos[which(reference.gene.annos[, 
            rga.gene.col] == gene), rga.anno.col]))
        pred.annos <- process.ref.annos.funk(unique(ipr.annos[which(ipr.annos$V4 == 
            gene), "V2"]))
        ipr.f1.df <- f1Score(ref.annos, pred.annos)
        ipr.f1.df$gene <- gene
        ipr.f1.df
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

#' Function to compute the performance F1-scores of Best Blast protein function
#' predictions using the Gene Ontology. The Best Blast method simply passes the
#' GO Term annotations found for the best Blast Hit to the respective query
#' protein. 
#'
#' @param best.blast.res.tbl An instance of \code{data.frame} result of
#' invoking \code{MapMan2GO::extractBestBlastHits}.
#' @param blast.hit.annos An instance of \code{data.frame} holding the GO
#' annotations for the proteins in the Blast database. Must have at least two
#' columns, one with the Blast Hit protein identifier (expected to be in short
#' Uniprot format), and another holding the GO Term identifiers. See this
#' packages Vignette, section 'Prepare UniprotKB Gene Ontology Annotations
#' (GOA) for processing in R.'.
#' @param bbrt.query.col An integer or string indicating the column of argument
#' \code{best.blast.res.tbl} in which to find the query protein identifiers.
#' Default is \code{getOption('MapMan2GO.best.blast.tbl.query.col',
#' 'query.ukb.short.id')}.
#' @param bbrt.hit.col An integer or string indicating the column of
#' argument \code{best.blast.res.tbl} in which to find the hit protein
#' identifiers. Default is \code{getOption('MapMan2GO.best.blast.tbl.hit.col',
#' 'hit.ukb.short.id')}.
#' @param bha.prot.col An integer or string indicating the column of argument
#' \code{blast.hit.annos} in which to find the hit protein identifiers.
#' Default is \code{getOption('MapMan2GO.blast.hit.goas.prot.col', 3)}.
#' @param bha.go.col An integer or string indicating the column of argument
#' \code{blast.hit.annos} in which to find the annotated GO Term
#' identifiers. Default is \code{getOption('MapMan2GO.blast.hit.goas.go.col',
#' 2)}.
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
#' MapMan2GO::addAncestors)}.
#'
#' @return An instance of \code{base::data.frame} with the following numeric
#' columns: precision, recall, true.pos, false.pos, f1.score, mcc, gene.
#' @export
bestBlastF1Score <- function(best.blast.res.tbl, blast.hit.annos, bbrt.query.col = getOption("MapMan2GO.best.blast.tbl.query.col", 
    "query.ukb.short.id"), bbrt.hit.col = getOption("MapMan2GO.best.blast.tbl.hit.col", 
    "hit.ukb.short.id"), bha.prot.col = getOption("MapMan2GO.blast.hit.goas.prot.col", 
    3), bha.go.col = getOption("MapMan2GO.blast.hit.goas.go.col", 2), reference.gene.annos = getOption("MapMan2GO.reference.gene.annos", 
    ukb.ref.goas), rga.gene.col = getOption("MapMan2GO.rga.gene.col", "V4"), rga.anno.col = getOption("MapMan2GO.rga.anno.col", 
    "V2"), process.ref.annos.funk = getOption("MapMan2GO.process.ref.annos.funk", 
    MapMan2GO::addAncestors)) {
    unique.reference.prot.ids <- unique(reference.gene.annos$V4)
    Reduce(rbind, mclapply(1:length(unique.reference.prot.ids), function(i) {
        x <- unique.reference.prot.ids[i]
        query.prot.id.coord <- which(best.blast.res.tbl[, bbrt.query.col] == x)
        query.prot.id <- best.blast.res.tbl[query.prot.id.coord, bbrt.query.col]
        hit.prot.id.coord <- which(best.blast.res.tbl[, bbrt.query.col] == x)
        hit.prot.id <- best.blast.res.tbl[hit.prot.id.coord, bbrt.hit.col]
        reference.annos.coord <- which(reference.gene.annos[, rga.gene.col] == 
            x)
        reference.annos <- process.ref.annos.funk(reference.gene.annos[reference.annos.coord, 
            rga.anno.col])
        if (length(hit.prot.id) > 0 & length(query.prot.id) > 0) {
            predicted.annos.coord <- which(blast.hit.annos[, bha.prot.col] == hit.prot.id)
            predicted.annos <- process.ref.annos.funk(blast.hit.annos[predicted.annos.coord, 
                bha.go.col])
        } else {
            predicted.annos <- character(0)
        }
        f1.df <- f1Score(reference.annos, predicted.annos)
        f1.df$gene <- x
        f1.df
    }))
}


#' Function testing \code{MapMan2GO::bestBlastF1Score}.
#'
#' @return \code{TRUE} if and only if all tests pass.
#' @export
testBestBlastF1Score <- function() {
    bb.res <- data.frame(query.ukb.short.id = LETTERS[1:3], hit.ukb.short.id = letters[10:12], 
        stringsAsFactors = FALSE)
    ref.annos <- data.frame(V5 = bb.res$query.ukb.short.id, V2 = paste("GO:", 1:3, 
        sep = ""), stringsAsFactors = FALSE)
    hit.annos <- data.frame(V1 = rep(NA, 3), V2 = paste("GO:", 1:3, sep = ""), 
        V3 = bb.res$hit.ukb.short.id, stringsAsFactors = FALSE)
    f1.a <- bestBlastF1Score(bb.res, hit.annos, reference.gene.annos = ref.annos, 
        process.ref.annos.funk = identity)
    t.1 <- sum(f1.a$f1.score) == 3
    hit.annos.b <- data.frame(V1 = rep(NA, 3), V2 = paste("GO:", 4:6, sep = ""), 
        V3 = bb.res$hit.ukb.short.id, stringsAsFactors = FALSE)
    f1.b <- bestBlastF1Score(bb.res, hit.annos.b, reference.gene.annos = ref.annos, 
        process.ref.annos.funk = identity)
    t.2 <- sum(f1.b$f1.score) == 0
    all(c(t.1, t.2))
}
