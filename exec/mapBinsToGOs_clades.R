require(MapMan2GO)
options(mc.cores = detectCores())

message("USAGE: Rscript path/2/MapMan2GO/exec/mapBinsToGOs.R path/2/MapMan2GO [/path/2/SwissProtIDs] [SwissProtIDs_subset_name]")

input.args <- commandArgs(trailingOnly = TRUE)

#' When a SwissProtIDs file is provided, perform the analysis for both, the whole SwissProt and the subset

#' If subset is given,
#' generate two mm.2.go, one called 'all' and one called 'subset' or input.arg[[3]] '... [SwissProtIDs subset_name]'
ref.prot.sets <- list(all = mm.bins.vs.sprot)
if (length(input.args) > 1) {
    subset.name <- input.args[[3]]
    clade.ID <- unlist((fread(input.args[[2]], header = FALSE, data.table = FALSE, 
        stringsAsFactors = FALSE)), use.names = FALSE)
    ref.prot.sets[subset.name] <- list(mm.bins.vs.sprot[which(mm.bins.vs.sprot$Swissprot.Short.ID %in% 
        clade.ID), ])
}

#' Run the algorithm mapping MapManBins to Gene Ontology Term Annotations
#' (GOAs) for each reference protein set:
for (ref.set.name in names(ref.prot.sets)) {
    mmvs.work <- ref.prot.sets[[ref.set.name]]
    
    #' Map MapMan-Bins to compound Gene Ontology Term Annotations (GOA):
    mm.leaf.bins <- unique(mmvs.work$MapManBin)
    mm.2.go <- setNames(mclapply(mm.leaf.bins, compoundGoAnnotationEntropy), mm.leaf.bins)
    mm.2.go.df <- Reduce(rbind, mclapply(names(mm.2.go), function(x) {
        y <- mm.2.go[[x]]
        data.frame(MapManBin = x, MapManBin.GO = y[["MapManBin.GO"]], Shannon.Entropy = y[["Shannon.Entropy"]], 
            n.GO = y[["n.GO"]], n.genes = y[["n.genes"]], median.n.GO = y[["median.n.GO"]], 
            stringsAsFactors = FALSE)
    }))
    #' Add a full description for each MapMan-Bin's GOA including GO-Term names:
    go.terms.not.in.db <- c()
    mm.2.full.desc <- Reduce(rbind, mclapply(names(mm.2.go), function(m.b) {
        m.b.gos <- Reduce(intersect, mm.2.go[[m.b]]$genes.goa)
        Reduce(rbind, mclapply(m.b.gos, function(g.id) {
            # 
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
    mm.fst <- read.fasta(file.path(path.package("MapMan2GO"), "mapman4.fasta"), 
        seqtype = "AA", strip.desc = TRUE, as.string = TRUE)
    mm.desc.df <- unique(Reduce(rbind, mclapply(mm.fst, function(x) {
        x.data <- strsplit(attr(x, "Annot"), " \\| ")[[1]]
        data.frame(MapManBin = x.data[[1]], Description = x.data[[2]], stringsAsFactors = FALSE)
    })))
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
    fl.name <- paste(ref.set.name, "NumberOfMapManBinsSharingIdentGOAsHist.pdf", 
        sep = "_")
    pdf(file.path(input.args[[1]], "inst", fl.name))
    plotDistAsHistAndBox(mm.2.go.df$n.GO, "Number of GO Terms per MapMan-Bin GOA")
    dev.off()
    #' - Histogram of Entropies
    fl.name <- paste(ref.set.name, "MapManBinGoaEntropyHist.pdf", 
        sep = "_")
    pdf(file.path(input.args[[1]], "inst", fl.name))
    plotDistAsHistAndBox(mm.2.go.df$Shannon.Entropy, "Shannon Entropy of compound GO Annotations per MapMan-Bin")
    dev.off()
    #' - Histogram of Sizes in terms of number of genes
    fl.name <- paste(ref.set.name, "GenesPerMapManBinHist.pdf", 
        sep = "_")
    pdf(file.path(input.args[[1]], "inst", fl.name))
    plotDistAsHistAndBox(mm.2.go.df$n.genes, "Number of genes per MapMan-Bin")
    dev.off()
    #' - Number of genes vs Number of GO Terms in the Bin-GOA:
    fl.name <- paste(ref.set.name, "GenesPerMapManBinHist.pdf", 
        sep = "_")
    pdf(file.path(input.args[[1]], "inst", fl.name))
    plot(mm.2.go.df$n.genes, mm.2.go.df$n.GO, xlab = "Number of genes per MapMan-Bin", 
        ylab = "Number of GO Terms in MapMan-Bin-GOA", pch = 20)
    dev.off()
    #' - Mutual information Histogram
#    fl.name <- paste(ref.set.name, "DistributionOfMutualInformationBetweenBinGoaAndReferenceGoas.pdf", 
#        sep = "_")
#    pdf(file.path(input.args[[1]], "inst", fl.name))
#    plotDistAsHistAndBox(mm.2.go.df$mutual.information, "Mutual Information between Bin GOA and reference GOAs [bits]")
#    dev.off()
    #' - Histogram of number of MapMan-Bins sharing identical GOAs
    #' @Yaxal: Example to be used in all cases of writing files (PDF, tables, etc):
    fl.name <- paste(subset.name, "NumberOfMapManBinsSharingIdentGOAsHist.pdf", 
        sep = "_")
    pdf(file.path(input.args[[1]], "inst", fl.name))
    x <- as.numeric(table(mm.2.go.df[which(mm.2.go.df$MapManBin.GO != ""), "MapManBin.GO"]))
    hist(x, col = "lightgrey", xlab = "Number of Bins sharing identical GOAs", 
        main = "Histogram of number of MapManBins sharing identical GOAs")
    dev.off()
    
    
    #' Save results:
    fl.name <- paste(ref.set.name, "MapManBins2GO.txt", 
        sep = "_")
    write.table(mm.2.full.desc, file.path(input.args[[1]], "inst", fl.name), 
        sep = "\t", row.names = FALSE)
    
   #' Save binary results!
    #' If we deal with more than a single reference protein set, and we are in
    #' the iteration that is not doing 'all', rename the results appending the
    #' ref.set.name. This preserves backward compatibility, and enables storing
    #' many mappings in a single RData.
    rdata.file <- file = file.path(input.args[[1]], "data", "MapManBins2GO.RData")
    if (ref.set.name != "all") {
        ref.set.results <- c(paste("mm.2.go", ref.set.name, sep = "."), paste("mm.2.go.df", 
            ref.set.name, sep = "."), paste("mm.2.full.desc", ref.set.name, sep = "."), 
            paste("mm.desc.df ", ref.set.name, sep = "."))
        assign(ref.set.results[[1]], mm.2.go)
        assign(ref.set.results[[2]], mm.2.go.df)
        assign(ref.set.results[[3]], mm.2.full.desc)
        assign(ref.set.results[[4]], mm.desc.df)
        appendToRData(list = ref.set.results, file = rdata.file)
    } else {
        save(mm.2.go, mm.2.go.df, mm.2.full.desc, mm.desc.df, file = rdata.file)
    }
}


message("DONE")

