require(MapMan2GO)
options(mc.cores = detectCores())

message("USAGE: Rscript path/2/MapMan2GO/exec/mapBinsToGOs.R path/2/MapMan2GO")

input.args <- commandArgs(trailingOnly = TRUE)

#' Map MapMan-Bins to compound Gene Ontology Term Annotations (GOA):
mm.leaf.bins <- unique(mm.bins.vs.sprot$MapManBin)
mm.2.go <- setNames(mclapply(mm.leaf.bins, compoundGoAnnotationEntropy), mm.leaf.bins)
mm.2.go.df <- Reduce(rbind, mclapply(names(mm.2.go), function(x) {
    y <- mm.2.go[[x]]
    data.frame(MapManBin = x, MapManBin.GO = y[["MapManBin.GO"]], Shannon.Entropy = y[["Shannon.Entropy"]], 
        mutual.information = y[["mutual.information"]], n.GO = y[["n.GO"]], n.genes = y[["n.genes"]], 
        median.n.GO = y[["median.n.GO"]], stringsAsFactors = FALSE)
}))


#' Some statistics:
#' - Histogram of Number of GO Terms per MapMan-Bin GOA
pdf(file.path(input.args[[1]], "inst", "NumberOfGoTermsPerMapManBinGoaHist.pdf"))
plotDistAsHistAndBox(mm.2.go.df$n.GO, "Number of GO Terms per MapMan-Bin GOA")
dev.off()
#' - Histogram of Entropies
pdf(file.path(input.args[[1]], "inst", "MapManBinGoaEntropyHist.pdf"))
plotDistAsHistAndBox(mm.2.go.df$Shannon.Entropy, "Shannon Entropy of compound GO Annotations per MapMan-Bin")
dev.off()
#' - Histogram of Sizes in terms of number of genes
pdf(file.path(input.args[[1]], "inst", "GenesPerMapManBinHist.pdf"))
plotDistAsHistAndBox(mm.2.go.df$n.genes, "Number of genes per MapMan-Bin")
dev.off()
#' - Number of genes vs Number of GO Terms in the Bin-GOA:
pdf(file.path(input.args[[1]], "inst", "NumberOfGenesVsNumberOfGoTermsInMapManBinGOA.pdf"))
plot(mm.2.go.df$n.genes, mm.2.go.df$n.GO, xlab = "Number of genes per MapMan-Bin", 
    ylab = "Number of GO Terms in MapMan-Bin-GOA", pch = 20)
dev.off()
#' - Mutual information Histogram
pdf(file.path(input.args[[1]], "inst", "DistributionOfMutualInformationBetweenBinGoaAndReferenceGoas.pdf"))
plotDistAsHistAndBox(mm.2.go.df$mutual.information, "Mutual Information between Bin GOA and reference GOAs [bits]")
dev.off()
#' - Histogram of number of MapMan-Bins sharing identical GOAs
pdf(file.path(input.args[[1]], "inst", "NumberOfMapManBinsSharingIdentGOAsHist.pdf"))
x <- as.numeric(table(mm.2.go.df[which(mm.2.go.df$MapManBin.GO != ""), "MapManBin.GO"]))
hist(x, col = "lightgrey", xlab = "Number of Bins sharing identical GOAs", main = "Histogram of number of MapManBins sharing identical GOAs")
dev.off()


#' Save results:
save(mm.2.go, mm.2.go.df, file = file.path(input.args[[1]], "data", "MapManBins2GO.RData"))


message("DONE")
