require(MapMan2GO)
options(mc.cores = detectCores())

message("USAGE: Rscript path/2/MapMan2GO/exec/mapBinsToGOs.R path/2/MapMan2GO")

input.args <- commandArgs(trailingOnly = TRUE)

#' Map MapMan-Bins to compound Gene Ontology Term Annotations (GOA):
mm.2.go <- setNames(mclapply(mm.bins.vs.sprot$MapManBin, compoundGoAnnotationEntropy), 
    mm.bins.vs.sprot$MapManBin)
mm.2.go.df <- Reduce(rbind, mclapply(names(mm.2.go), function(x) {
    y <- mm.2.go[[x]]
    data.frame(MapManBin = x, MapManBin.GO = y[["MapManBin.GO"]], Shannon.Entropy = y[["Shannon.Entropy"]], 
        n.GO = y[["n.GO"]], n.genes = y[["n.genes"]], median.n.GO = y[["median.n.GO"]], 
        stringsAsFactors = FALSE)
}))


#' Some statistics:
#' - Histogram of Number of GO Terms per MapMan-Bin GOA
pdf(file.path(input.args[[1]], "inst", "NumberOfGoTermsPerMapManBinGoaHist.pdf"))
hist(mm.2.go.df$n.GO, breaks = 20, col = "lightgrey", xlab = "Number of GO Terms per MapMan-Bin GOA", 
    main = "")
dev.off()
#' - Histogram of Entropies
pdf(file.path(input.args[[1]], "inst", "MapManBinGoaEntropyHist.pdf"))
hist(mm.2.go.df$Shannon.Entropy, breaks = 20, xlab = "Shannon Entropy of compound GO Annotations per MapMan-Bin", 
    main = "", col = "lightgrey")
dev.off()
#' - Histogram of Sizes in terms of number of genes
pdf(file.path(input.args[[1]], "inst", "GenesPerMapManBinHist.pdf"))
hist(mm.2.go.df$n.genes, breaks = 20, xlab = "Number of genes per MapMan-Bin", main = "", 
    col = "lightgrey")
dev.off()
#' - Number of genes vs Number of GO Terms in the Bin-GOA:
pdf(file.path(input.args[[1]], "inst", "NumberOfGenesVsNumberOfGoTermsInMapManBinGOA.pdf"))
plot(mm.2.go.df$n.genes, mm.2.go.df$n.GO, xlab = "Number of genes per MapMan-Bin", 
    ylab = "Number of GO Terms in MapMan-Bin-GOA", pch = 20)
dev.off()
#' - Median of GO Terms per Gene-GOA vs Number of GO Terms in Bin GOA
pdf(file.path(input.args[[1]], "inst", "MedianOfGoTermsPerGeneGoaVsNumberOfGoTermsInMapManBinGOA.pdf"))
plot(mm.2.go.df$median.n.GO, mm.2.go.df$n.GO, xlab = "Median of GO Terms per Gene-GOA per MapMan-Bin", 
    ylab = "Number of GO Terms in MapMan-Bin-GOA", pch = 20)
dev.off()


#' Save results:
save(mm.2.go, mm.2.go.df, file = file.path(input.args[[1]], "data", "MapManBins2GO.RData"))


message("DONE")
