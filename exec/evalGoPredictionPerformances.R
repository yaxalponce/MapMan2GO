require(MapMan2GO)
options(mc.cores = detectCores())

message("USAGE: Rscript path/2/MapMan2GO/exec/evalGoPredictionPerformances.R path/2/mercator_result.txt path/2/blast_result_table path/2/protein_identifiers_2_exclude path/2/preprocessed_uniprot_goa_table path/2/MapMan2GO")

input.args <- commandArgs(trailingOnly = TRUE)


#' Compute F1-Scores and Matthews correlation coefficient of competing GO predictors.
#' ##################################################################################

#' First consider the universe of all possible GO term annotations 'as is',
#' that is without adding ancestral GO terms.

ipr.annos <- ukb.ref.goas[!is.na(ukb.ref.goas$V3), ]
ipr.annos.no.ref.anc <- predictorPerformance(ipr.annos, pa.gene.col = "V4", pa.anno.col = "V2", 
    process.annos.funk = identity)

best.blast.tbl <- extractBestBlastHits(input.args[[2]], input.args[[3]])
goa.tbl <- fread(input.args[[4]], sep = "\t", data.table = FALSE, stringsAsFactors = FALSE, 
    quote = "", na.strings = "", header = FALSE)
blast.hit.goa <- goa.tbl[which(goa.tbl$V3 %in% best.blast.tbl$hit.ukb.short.id), 
    ]
rm(goa.tbl)
best.blast.pred <- bestBlastPredictions(best.blast.tbl, blast.hit.goa)
bb.annos.no.ref.anc <- predictorPerformance(best.blast.pred, pa.gene.col = "query", 
    pa.anno.col = "GO", process.annos.funk = identity)


#' Now consider the universe of all possible GO term annotations to include
#' ancestral terms. Also extend both reference GO term annotations as well as
#' predicted GO term annotations with ancestral termsc

options(MapMan2GO.performance.universe.annotations = ukb.ref.universe.gos.w.anc)

mercator.annos <- readMercatorResultTable(input.args[[1]], sanitize.accession = TRUE)
mercator.annos.f1 <- predictorPerformance(mercator.annos, pa.gene.col = "IDENTIFIER", 
    pa.anno.col = "MapManBin.GO", process.predicted.annos.funk = splitMapManBinGOAs, 
    reference.genes = tolower(ref.gene.ids), rga.gene.col = "V5")

ipr.annos.w.ref.anc <- predictorPerformance(ipr.annos, pa.gene.col = "V4", pa.anno.col = "V2")

bb.annos.w.ref.anc <- predictorPerformance(best.blast.pred, pa.gene.col = "query", 
    pa.anno.col = "GO")

#' Set the respective option back to default value:
options(MapMan2GO.performance.universe.annotations = NULL)


#' Plot histograms of F1-Scores and Matthews correlation coefficient:
scores <- c("precision", "recall", "false.pos.rate", "f1.score", "mcc")

#' - for mercator
for (score.i in scores) {
    pdf(file.path(input.args[[length(input.args)]], "inst", paste("mercator_", 
        score.i, "_Hist.pdf", sep = "")))
    plotDistAsHistAndBox(mercator.annos.f1[, score.i], main = paste("Mercator", 
        score.i, "distribution"), summary.as.title = TRUE)
    dev.off()
}


#' - for InterProScan
for (score.i in scores) {
    pdf(file.path(input.args[[length(input.args)]], "inst", paste("interProScanNoAncRefGoTerms_", 
        score.i, "_Hist.pdf", sep = "")))
    plotDistAsHistAndBox(ipr.annos.no.ref.anc[, score.i], main = paste("InterProScan", 
        score.i, "distribution"), summary.as.title = TRUE)
    dev.off()
}
for (score.i in scores) {
    pdf(file.path(input.args[[length(input.args)]], "inst", paste("interProScanWithAncRefGoTerms_", 
        score.i, "_Hist.pdf", sep = "")))
    plotDistAsHistAndBox(ipr.annos.w.ref.anc[, score.i], main = paste("InterProScan", 
        score.i, "distribution"), summary.as.title = TRUE)
    dev.off()
}

#' - for BestBlast
for (score.i in scores) {
    pdf(file.path(input.args[[length(input.args)]], "inst", paste("bestBlastWithAncRefGoTerms_", 
        score.i, "_Hist.pdf", sep = "")))
    plotDistAsHistAndBox(bb.annos.w.ref.anc[, score.i], main = paste("BestBlast", 
        score.i, "distribution"), summary.as.title = TRUE)
    dev.off()
}
for (score.i in scores) {
    pdf(file.path(input.args[[length(input.args)]], "inst", paste("bestBlastNoAncRefGoTerms_", 
        score.i, "_Hist.pdf", sep = "")))
    plotDistAsHistAndBox(bb.annos.no.ref.anc[, score.i], main = paste("BestBlast", 
        score.i, "distribution"), summary.as.title = TRUE)
    dev.off()
}


#' Plot Blast sequence similarity distribution among the queries for which Best
#' Blast achieved a non zero False Positive Rate (FPR):
bb.i <- which(bb.annos.w.ref.anc$false.pos.rate > 0)
bb.queries <- bb.annos.w.ref.anc[bb.i, "gene"]
pdf(file.path(input.args[[length(input.args)]], "inst", "bestBlastWithAncestorsSeqSimOnQueriesWithFPRgtZero_Hist.pdf"))
plotDistAsHistAndBox(best.blast.tbl[which(best.blast.tbl$query.ukb.short.id %in% 
    bb.queries), "V3"], main = "Blast Sequence Similarity for queries with FPR > 0", 
    summary.as.title = TRUE)
dev.off()



#' Scatterplots of of n.truth againt false positives, true positives and n.pred:
scatter.y <- c("n.pred", "false.pos", "true.pos")

#' - for mercator
for (scatter.i in scatter.y) {
    pdf(file.path(input.args[[length(input.args)]], "inst", paste("mercator_", 
        scatter.i, "_Scatter.pdf", sep = "")))
    plot(mercator.annos.f1$n.truth, mercator.annos.f1[, scatter.i], xlab = "n.truth", ylab = scatter.i, pch = 20, 
         main = paste("Mercator", scatter.i, "distribution"))
    dev.off()
}

#' - for Best Blast with ancestral GO Terms
for (scatter.i in scatter.y) {
    pdf(file.path(input.args[[length(input.args)]], "inst", paste("bestBlastWithAncRefGoTerms_", 
        scatter.i, "_Scatter.pdf", sep = "")))
    plot(bb.annos.w.ref.anc$n.truth, bb.annos.w.ref.anc[, scatter.i], xlab = "n.truth", ylab = scatter.i, pch = 20, 
         main = paste("Mercator", scatter.i, "distribution"))
    dev.off()
}

#' - for InterProScan with ancestral GO Terms
for (scatter.i in scatter.y) {
    pdf(file.path(input.args[[length(input.args)]], "inst", paste("interProScanWithAncRefGoTerms_", 
        scatter.i, "_Scatter.pdf", sep = "")))
    plot(ipr.annos.w.ref.anc$n.truth, ipr.annos.w.ref.anc[, scatter.i], xlab = "n.truth", ylab = scatter.i, pch = 20, 
         main = paste("Mercator", scatter.i, "distribution"))
    dev.off()
}



#' Save results:
save(mercator.annos, mercator.annos.f1, ipr.annos, ipr.annos.w.ref.anc, ipr.annos.no.ref.anc, 
    best.blast.tbl, bb.annos.no.ref.anc, bb.annos.w.ref.anc, file = file.path(input.args[[length(input.args)]], 
        "data", "predictionPerformances.RData"))


message("DONE")
