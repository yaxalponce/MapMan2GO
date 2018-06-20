require(MapMan2GO)
options(mc.cores = detectCores())

message("USAGE: Rscript path/2/MapMan2GO/exec/evalGoPredictionPerformances.R path/2/mercator_result.txt path/2/MapMan2GO")

input.args <- commandArgs(trailingOnly = TRUE)


#' Compute F1-Scores of competing GO predictors:
mercator.annos <- readMercatorResultTable(input.args[[1]])
mercator.annos <- cbind(mercator.annos, mercatorF1Score(mercator.annos))

ipr.annos <- ukb.ref.goas[!is.na(ukb.ref.goas$V3), ]
ipr.annos.no.ref.anc <- interProScanF1Scores(ipr.annos, process.ref.annos.funk = identity)
ipr.annos.w.ref.anc <- interProScanF1Scores(ipr.annos)


#' Plot histograms of F1-Scores:
#' - for mercator
pdf(file.path(input.args[[2]], "inst", "mercatorF1ScoresHist.pdf"))
message("WARNING: Currently only plotting Mercator F1 Scores for Bins, we had AA sequences for.")
mm.bins.w.aa.seqs <- unique(names(read.fasta(file.path(path.package("MapMan2GO"), 
    "mapman4.fasta"), seqtype = "AA", strip.desc = TRUE, as.string = TRUE)))
hist(mercator.annos[which(mercator.annos$BINCODE %in% mm.bins.w.aa.seqs), "f1.score"], 
    xlab = "F1-Scores", main = "Mercator's Performance", col = "lightgrey")
dev.off()

#' - for InterProScan
pdf(file.path(input.args[[2]], "inst", "interProScanNoAncRefGoTermsF1ScoresHist.pdf"))
hist(ipr.annos.no.ref.anc$f1.score, xlab = "F1-Scores", main = "InterProScan's Performance", 
    sub = "Reference GO Terms as is", col = "lightgrey")
dev.off()
pdf(file.path(input.args[[2]], "inst", "interProScanWithAncRefGoTermsF1ScoresHist.pdf"))
hist(ipr.annos.w.ref.anc$f1.score, xlab = "F1-Scores", main = "InterProScan's Performance", 
    sub = "Reference GO Terms include ancestors", col = "lightgrey")
dev.off()


#' Save results:
save(mercator.annos, ipr.annos.w.ref.anc, ipr.annos.no.ref.anc, file = file.path(input.args[[length(input.args)]], 
    "data", "predictionPerformances.RData"))


message("DONE")
