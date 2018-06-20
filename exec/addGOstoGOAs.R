require(MapMan2GO)

message("USAGE: Rscript path/2/MapMan2GO/addGOstoGOAs.R path/2/MapMan2GO")

input.args <- commandArgs(trailingOnly = TRUE)

words.to.filter <- fread(file.path(input.args[[1]], "inst", "words2filter"), header = FALSE)

#' Analysis on shared words for all the MapMan-Bins

#' Add GOs to MapMan-Bin GOAs that add new information
GOs.to.add.list <- c()

info.bins.words.all <- analyzeSharedWords(mm.2.go.df$MapManBin.GO, mm.2.go.df$MapManBin)

for (i in 1:length(mm.2.go.df$MapManBin)) {
    
    #' Union of protein GOAs
    a <- paste("mm.2.go$`", paste(mm.2.go.df$MapManBin[i], "`$genes.goa", sep = ""), 
        sep = "")
    GO.joined <- sort(unique(unlist(eval(parse(text = a)))))
    n <- which(info.bins.words.all$mm.2.go.df.MapManBin == mm.2.go.df$MapManBin[i])
    GOA.2.compare <- strsplit(mm.2.go.df$MapManBin.GO[n], ",")
    
    #' Include just words no shared before
    bin.description.2.compare <- strsplit(as.character(info.bins.words.all$bins.gos.no.shared.words[n]), 
        split = ",")
    
    #' Identify GOs not present in original GOA
    GOs.to.eval <- GO.joined[which(GO.joined %in% GOA.2.compare[[1]] == FALSE)]
    
    #' Split GOs name, search them in MapManBin description, if there's any, select the GO
    GO.names.extended <- c()
    GO.description.splited.extended <- c()
    GOs.to.add <- c()
    for (j in 1:length(GOs.to.eval)) {
        GO.names.extended[j] <- GO.OBO$name[GOs.to.eval[j]]
        GO.description.splited.extended <- unique(unlist(strsplit(GO.names.extended[j], 
            split = "[., ,_,:]")))
        if (any(sapply(GO.description.splited.extended, function(x) x %in% bin.description.2.compare[[1]]) == 
            TRUE) == TRUE) {
            GOs.to.add[j] <- GOs.to.eval[j]
        } else {
            GOs.to.add[j] <- NA
        }
    }
    if (!is.null(GOs.to.add) || !is.na(GOs.to.add)) {
        GOs.to.add.list[i] <- paste(sort(unique(unlist(GOs.to.add))), collapse = ",")
    } else {
        GOs.to.add.list[i] <- NA
    }
}

#' Retrieve the Original GOAs to add the new GOs
bins.GOs.to.add <- c()
bins.GOs.to.add <- data.frame(mm.2.go.df$MapManBin, GOs.to.add.list)
bins.to.extend <- bins.GOs.to.add$mm.2.go.df.MapManBin[which(bins.GOs.to.add$GOs.to.add.list != 
    "")]

new.goa <- c()
n.new.GOs <- c()
for (i in 1:length(bins.to.extend)) {
    goao <- mm.2.go.df$MapManBin.GO[which(mm.2.go.df$MapManBin == bins.to.extend[i])]
    goe <- bins.GOs.to.add$GOs.to.add.list[which(mm.2.go.df$MapManBin == bins.to.extend[i])]
    ng <- union(goao, goe)
    ng <- paste(ng, collapse = ",")
    ng <- sort(unlist(strsplit(ng, split = ",")))
    n.new.GOs[i] <- length(ng)
    new.goa[i] <- paste(ng, collapse = ",")
}

#' Substitute original GOAs and n.GOs with extended GOAs in data frame

for (i in 1:length(bins.to.extend)) {
    coords <- which(mm.2.go.df$MapManBin == bins.to.extend[i])
    mm.2.go.df$MapManBin.GO[coords] <- new.goa[i]
    mm.2.go.df$n.GO[coords] <- n.new.GOs[i]
}

#' Identify GOAs duplicated
coords.duplicated.GOAs <- which(duplicated(mm.2.go.df$MapManBin.GO) | duplicated(mm.2.go.df$MapManBin.GO, 
    fromLast = TRUE))
duplicated.GOAs <- mm.2.go.df$MapManBin.GO[coords.duplicated.GOAs]
duplicated.bins <- mm.2.go.df$MapManBin[coords.duplicated.GOAs]

#' Level at which Bins are sharing GOAs
#' Apply just for duplicated GOAs in different MapMan-Bins

duplicated.bins.combn <- combn(duplicated.bins, 2)
names(duplicated.GOAs) <- duplicated.bins

n <- length(duplicated.bins.combn)/2

MRCA.D <- c()
MRCA.Dist <- c()
for (j in 1:n) {
    if ((duplicated.GOAs[names(duplicated.GOAs) == duplicated.bins.combn[1, j]]) == 
        (duplicated.GOAs[names(duplicated.GOAs) == duplicated.bins.combn[2, j]])) {
        a <- unlist(strsplit(duplicated.bins.combn[1, j], ".", fixed = TRUE))
        b <- unlist(strsplit(duplicated.bins.combn[2, j], ".", fixed = TRUE))
        
        COUNT = 1
        e <- min(length(a), length(b))
        while (a[COUNT] == b[COUNT] && COUNT <= e) {
            COUNT = COUNT + 1
        }
        d <- max(length(a), length(b))
        MRCA.D[j] <- COUNT - 1
        MRCA.Dist[j] <- d - (COUNT - 1)
    }
}

#' Prepare data to graph number of GO sharing words with Bins
GO.description.splited <- c()
n.GO.sharing <- c()
for (j in 1:length(mm.2.go.df$MapManBin.GO)) {
    splited.st <- strsplit(mm.2.go.df$MapManBin.GO[j], ",")
    bin.description.splited <- strsplit(gsub("\\(", "", gsub("\\)", "", gsub(",", 
        "", mm.desc.df$Description[which(mm.desc.df$MapManBin == mm.2.go.df$MapManBin[j])]))), 
        split = "[., ,_,:]")

    bin.description.splited <- if (!length(bin.description.splited)) {
		    c("na")
    } else tolower(bin.description.splited[[1]])
    
    # bin.description.splited <- tolower(bin.description.splited[[1]])
    bin.description.splited <- bin.description.splited[which(bin.description.splited %in% 
        unlist(words.to.filter) == "FALSE")]

    GO.names <- c()
    COUNT = 0
    for (i in 1:length(splited.st[[1]])) {
        GO.names[i] <- GO.OBO$name[splited.st[[1]][i]]
	GO.description.splited <- unique(unlist(strsplit(GO.names, split = "[., ,_,:]")))
	words.grepl <- sapply(bin.description.splited, function(x) any(grepl(x, 
            GO.description.splited[i], ignore.case = TRUE)))
	if(any(words.grepl == TRUE)) {
            COUNT = COUNT + 1
        }
    n.GO.sharing[j] <- COUNT
    }
}
n.GO.inters <- mm.2.go.df$n.GO

n.GO.percent.sharing <- c()
for (i in 1:length(n.GO.inters)) {
      n.GO.percent.sharing[i] <- n.GO.sharing[i]/n.GO.inters[i]
}

##' - Histogram of Entropies
pdf(file.path(input.args[[1]], "inst", "MapManBinGoaEntropyHist.pdf"))
plotDistAsHistAndBox(mm.2.go.df$Shannon.Entropy, "Shannon Entropy of compound GO Annotations per MapMan-Bin")
dev.off()

##' - Histogram of Sizes in terms of number of genes
pdf(file.path(input.args[[1]], "inst", "GenesPerMapManBinHist.pdf"))
plotDistAsHistAndBox(mm.2.go.df$n.genes, "Number of genes per MapMan-Bin")
dev.off()

##' - Number of genes vs Number of GO Terms in the Bin-GOA:
pdf(file.path(input.args[[1]], "inst", "NumberOfGenesVsNumberOfGoTermsInMapManBinGOA.pdf"))
plot(mm.2.go.df$n.genes, mm.2.go.df$n.GO, xlab = "Number of genes per MapMan-Bin",
    ylab = "Number of GO Terms in MapMan-Bin-GOA", pch = 20)
dev.off()

##' - Mutual information Histogram
pdf(file.path(input.args[[1]], "inst", "DistributionOfMutualInformationBetweenBinGoaAndReferenceGoas.pdf"))
plotDistAsHistAndBox(mm.2.go.df$mutual.information, "Mutual Information between Bin GOA and reference GOAs [bits]")
dev.off()

#' Plot the percent of shared words among the MapMan-Bins and GOAs descriptions
pdf(file.path(input.args[[1]], "inst", "MapManBin.GOAs.Shared.Words.pdf"))
plotDistAsHistAndBox(info.bins.words.all$percent.shared.words, "Percent of shared words among MapMan-Bins and GOAs")
dev.off()

#' - Number of genes vs Number of GO Terms in the Bin-GOA:
pdf(file.path(input.args[[1]], "inst", "NumberOfGenesVsNumberOfGoTermsInMapManBinGOA.pdf"))
plot(mm.2.go.df$n.genes, mm.2.go.df$n.GO, xlab = "Number of genes per MapMan-Bin", 
    ylab = "Number of GO Terms in MapMan-Bin-GOA", pch = 20)
dev.off()

#' - Histogram of Number of GO Terms per MapMan-Bin GOA
pdf(file.path(input.args[[1]], "inst", "NumberOfGoTermsPerMapManBinGoaHist.pdf"))
plotDistAsHistAndBox(mm.2.go.df$n.GO, "Number of GO Terms per MapMan-Bin GOA")
dev.off()

#' - Histogram of number of MapMan-Bins sharing identical GOAs
pdf(file.path(input.args[[1]], "inst", "NumberOfMapManBinsSharingIdentGOAsHist.pdf"))
x <- as.numeric(table(mm.2.go.df[which(mm.2.go.df$MapManBin.GO != ""), "MapManBin.GO"]))
hist(x, col = "lightgrey", xlab = "Number of Bins sharing identical GOAs", main = "Histogram of number of MapManBins sharing identical GOAs")
dev.off()

#' - Histogram of most recent common ancestor of depth (MRCA-D)
pdf(file.path(input.args[[1]], "inst", "MRCA-D.Hist.pdf"))
hist(MRCA.D, col = "lightgrey", xlab = "MRCA-D", main = "Histogram of most recent common ancestor of depth")
dev.off()

#' - Histogram of the distance to the most recent common ancestor  (MRCA-Dist)
pdf(file.path(input.args[[1]], "inst", "MRCA-Dist.Hist.pdf"))
hist(MRCA.Dist, col = "lightgrey", xlab = "MRCA-Dist", main = "Histogram of the distance to the most recent common ancestor")
dev.off()

#' - Scatterplot of GO terms sharing words with MapMan-Bin description
pdf(file.path(input.args[[1]], "inst", "GO.Terms.Sharing.MapMan-Bin.Description.pdf"))
  plot(n.GO.sharing, n.GO.inters, main="GO Terms Sharing MapMan-Bin Description", xlab="Number of GO-Terms sharing words with MapMan Description", ylab="Number of GO Terms", pch = 20)
dev.off()

#' - Histogram of percent of GO terms that share words with the MapMan-Bin description
pdf(file.path(input.args[[1]], "inst", "GO.Terms.Sharing.Words.Percent.pdf"))
hist(n.GO.percent.sharing, col = "lightgrey", xlab= "Percent of GO terms sharing at least one word", main = "GO terms sharing words with the MapMan-Bin description")
dev.off()


#' Save results:
save(mm.2.go, mm.2.go.df, mm.2.full.desc, mm.desc.df, file = file.path(input.args[[1]], "data", "MapManBins2GO.RData"))
