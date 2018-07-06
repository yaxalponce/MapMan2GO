require(MapMan2GO)
options(mc.cores = detectCores())

message("USAGE: Rscript path/2/MapMan2GO/exec/evalGoPredictionPerformances.R path/2/mercator_result.txt path/2/blast_result_table path/2/protein_identifiers_2_exclude path/2/preprocessed_uniprot_goa_table path/2/MapMan2GO")

input.args <- commandArgs(trailingOnly = TRUE)


#' Compute F1-Scores of competing GO predictors.
#' #############################################

#' First consider the universe of all possible GO term annotations 'as is',
#' that is without adding ancestral GO terms.

ipr.annos <- ukb.ref.goas[!is.na(ukb.ref.goas$V3), ]
ipr.annos.no.ref.anc <- interProScanF1Scores(ipr.annos, process.ref.annos.funk = identity)

best.blast.tbl <- extractBestBlastHits(input.args[[2]], input.args[[3]])
goa.tbl <- fread(input.args[[4]], sep = "\t", data.table = FALSE, stringsAsFactors = FALSE, 
    quote = "", na.strings = "", header = FALSE)
blast.hit.goa <- goa.tbl[which(goa.tbl$V3 %in% best.blast.tbl$hit.ukb.short.id), 
    ]
rm(goa.tbl)
bb.annos.no.ref.anc <- bestBlastF1Score(best.blast.tbl, blast.hit.goa, process.ref.annos.funk = identity)


#' Now consider the universe of all possible GO term annotations to include
#' ancestral terms. Also extend both reference GO term annotations as well as
#' predicted GO term annotations with ancestral terms.

options(MapMan2GO.performance.universe.annotations = ukb.ref.universe.gos.w.anc)

mercator.annos <- readMercatorResultTable(input.args[[1]], sanitize.accession = TRUE)
mercator.annos.f1 <- mercatorF1Score(mercator.annos)

ipr.annos.w.ref.anc <- interProScanF1Scores(ipr.annos)

bb.annos.w.ref.anc <- bestBlastF1Score(best.blast.tbl, blast.hit.goa)

#' Set the respective option back to default value:
options(MapMan2GO.performance.universe.annotations = NULL)


#' Plot histograms of F1-Scores:
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


#' For the above 'difficult to predict protein functions' (bb.queries) plot the
#' histograms of MCC scores into one graph:
pdf(file.path(input.args[[length(input.args)]], "inst", "difficultToPredictProtFunc_all_MCC_Hist.pdf"))
colors <- brewer.pal(3, "Dark2")
col.alpha <- addAlpha(colors)
bb.queries.all.mcc <- c(bb.annos.w.ref.anc[bb.i, "mcc"], mercator.annos.f1[which(mercator.annos.f1$gene %in% 
    tolower(bb.queries)), "mcc"], ipr.annos.w.ref.anc[which(ipr.annos.w.ref.anc$gene %in% 
    bb.queries), "mcc"])
x.min <- min(bb.queries.all.mcc, na.rm = TRUE)
x.max <- max(bb.queries.all.mcc, na.rm = TRUE)
bb.h <- hist(bb.annos.w.ref.anc[bb.i, "mcc"], plot = FALSE)
mer.h <- hist(mercator.annos.f1[which(mercator.annos.f1$gene %in% tolower(bb.queries)), "mcc"], plot = FALSE)
ipr.h <- hist(ipr.annos.w.ref.anc[which(ipr.annos.w.ref.anc$gene %in% bb.queries), 
    "mcc"], plot = FALSE)
y.max <- max(c(bb.h$counts, mer.h$counts, ipr.h$counts), na.rm = TRUE)
plot(bb.h, main = "Comparison of MCC for difficult to predict protein functions", 
    xlim = c(x.min, x.max), ylim = c(0, y.max), xlab = "Matthew's correlation coefficient (MCC)", 
    col = col.alpha[[1]], border = colors[[1]])
plot(mer.h, col = col.alpha[[2]], border = colors[[2]], add = TRUE)
plot(ipr.h, col = col.alpha[[3]], border = colors[[3]], add = TRUE)
dev.off()


#' Save results:
save(mercator.annos, mercator.annos.f1, ipr.annos, ipr.annos.w.ref.anc, ipr.annos.no.ref.anc, 
    best.blast.tbl, bb.annos.no.ref.anc, bb.annos.w.ref.anc, file = file.path(input.args[[length(input.args)]], 
        "data", "predictionPerformances.RData"))


message("DONE")
