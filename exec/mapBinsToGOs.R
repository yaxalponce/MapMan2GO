require(MapMan2GO)
options(mc.cores = detectCores())

message("USAGE: Rscript path/2/MapMan2GO/exec/mapBinsToGOs.R path/2/MercatorRawResults.tsf path/2/MapMan2GO")

input.args <- commandArgs(trailingOnly = TRUE)

#' Map MapMan-Bins to compound Gene Ontology Term Annotations (GOA):
mm.leaf.bins <- unique(mm.bins.vs.sprot$MapManBin)
mm.2.go <- setNames(mclapply(mm.leaf.bins, compoundGoAnnotationEntropy), mm.leaf.bins)  
mm.2.go.df <- Reduce(rbind, mclapply(names(mm.2.go), function(x) {
    y <- mm.2.go[[x]]
    data.frame(MapManBin = x, MapManBin.GO = y[["MapManBin.GO"]], Shannon.Entropy = y[["Shannon.Entropy"]], 
        n.GO = y[["n.GO"]], n.genes = y[["n.genes"]], 
        median.n.GO = y[["median.n.GO"]], stringsAsFactors = FALSE)
}))
#' Add a full description for each MapMan-Bin's GOA including GO-Term names:
go.terms.not.in.db <- c()
mm.2.full.desc <- Reduce(rbind, mclapply(names(mm.2.go), function(m.b) {
    m.b.gos <- Reduce(intersect, mm.2.go[[m.b]]$genes.goa)
    Reduce(rbind, mclapply(m.b.gos, function(g.id) {    #
        if (g.id %in% GO.OBO$id) {
            g.name <- GO.OBO$name[[g.id]]
            g.ancestors <- GO.OBO$ancestors[[g.id]]
            g.depth <- length(g.ancestors)
            g.ont <- if (g.depth > 1) {
                GO.OBO$name[[g.ancestors[[1]]]]
            } else {
                g.name
            }
            data.frame(MapManBin = m.b, GO.Term = g.id, GO.Name = g.name, GO.Ontolgy = g.ont, 
                GO.depth = g.depth, stringsAsFactors = FALSE)
        } else {
            go.terms.not.in.db <- c(go.terms.not.in.db, g.id)
            data.frame(MapManBin = m.b, GO.Term = g.id, GO.Name = NA, GO.Ontolgy = NA, 
                GO.depth = NA, stringsAsFactors = FALSE)
        }
    }))
}))
#' Warn about missing GO Terms:
if (length(go.terms.not.in.db) > 0) {
    message("WARNING: The following GO Terms were assigned to MapManBin(s), but had no entry in the installed 'GO.db' package:\n", 
        paste(sort(unique(go.terms.not.in.db)), collapse = ", "))
}


#' Add the MapMan Bins Descriptions to the above table:
    mr.full <- readMercatorResultTable(input.args[[1]], add.go.terms = FALSE)
    mm.desc.df <- unique(mr.full[, c("BINCODE", "NAME")])
    names(mm.desc.df) <- c("MapManBin", "Description")

    # mm.2.full.desc$Bin.Description <-
    # as.character(unlist(mclapply(mm.2.full.desc$MapManBin, function(x) {
    # mm.desc.df[which(mm.desc.df$MapManBin == x), 'Description'] })))
    
    value = integer(0)
    mm.2.full.desc$Bin.Description <- as.character(unlist(mclapply(mm.2.full.desc$MapManBin, 
        function(x) {
            if (identical(which(mm.desc.df$MapManBin == x), value)) {
                mm.2.full.desc$Bin.Description <- NA
            } else {
                mm.desc.df[which(mm.desc.df$MapManBin == x), "Description"]
            }
            
        })))


#' Some statistics:
#' - Histogram of Number of GO Terms per MapMan-Bin GOA
pdf(file.path(input.args[[length(input.args)]], "inst", "NumberOfGoTermsPerMapManBinGoaHist.pdf"))
plotDistAsHistAndBox(mm.2.go.df$n.GO, "Number of GO Terms per MapMan-Bin GOA")
dev.off()
#' - Histogram of Entropies
pdf(file.path(input.args[[length(input.args)]], "inst", "MapManBinGoaEntropyHist.pdf"))
plotDistAsHistAndBox(mm.2.go.df$Shannon.Entropy, "Shannon Entropy of compound GO Annotations per MapMan-Bin")
dev.off()
#' - Histogram of Sizes in terms of number of genes
pdf(file.path(input.args[[length(input.args)]], "inst", "GenesPerMapManBinHist.pdf"))
plotDistAsHistAndBox(mm.2.go.df$n.genes, "Number of genes per MapMan-Bin")
dev.off()
#' - Number of genes vs Number of GO Terms in the Bin-GOA:
pdf(file.path(input.args[[length(input.args)]], "inst", "NumberOfGenesVsNumberOfGoTermsInMapManBinGOA.pdf"))
plot(mm.2.go.df$n.genes, mm.2.go.df$n.GO, xlab = "Number of genes per MapMan-Bin", 
    ylab = "Number of GO Terms in MapMan-Bin-GOA", pch = 20)
dev.off()
#' - Histogram of number of MapMan-Bins sharing identical GOAs
pdf(file.path(input.args[[length(input.args)]], "inst", "NumberOfMapManBinsSharingIdentGOAsHist.pdf"))
x <- as.numeric(table(mm.2.go.df[which(mm.2.go.df$MapManBin.GO != ""), "MapManBin.GO"]))
hist(x, col = "lightgrey", xlab = "Number of Bins sharing identical GOAs", main = "Histogram of number of MapManBins sharing identical GOAs")
dev.off()


#' Save results:
write.table(mm.2.full.desc, file.path(input.args[[length(input.args)]], "inst", "MapManBins2GO.txt"), 
    sep = "\t", row.names = FALSE)
save(mm.2.go, mm.2.go.df, mm.2.full.desc, mm.desc.df, file = file.path(input.args[[length(input.args)]], 
    "data", "MapManBins2GO.RData"))


message("DONE")
