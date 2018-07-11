require(MapMan2GO)
options(mc.cores = detectCores())      

message("USAGE: Rscript /path/2/MapMan2GO/mapBinsRevision.R /path/2/MapMan2GO ")

input.args <- commandArgs(trailOnly = TRUE)

#' Retrieve the GOA and classify the GO Terms as 'Used' or 'Not Used'
sh.entropy <- c()
trusted.used <- c()
untrusted.used <- c()
trusted.no.used <- c()
untrusted.no.used <- c()

for (i in 1:length(mm.2.go.df$MapManBin)) {
    
    #' used as GOs in the GOA before the extension
    a <- paste("mm.2.go$`", paste(mm.2.go.df$MapManBin[i], "`$MapManBin.GO", sep = ""), 
        sep = "")
    o.go <- sort(unique(unlist(eval(parse(text = a)))))
    used.gos <- unlist(strsplit(o.go, split = ","))
    
    #' used as GOs included in the GOA
    b <- paste("mm.2.go$`", paste(mm.2.go.df$MapManBin[i], "`$genes.goa", sep = ""), 
        sep = "")
    GO.joined <- sort(unique(unlist(eval(parse(text = b)))))
    no.used.gos <- setdiff(GO.joined, used.gos)
    sh.entropy[i] <- entropy(table(as.character(no.used.gos)))

#' Analyze Evidence Codes
    eco.used.gos.all <- unique(unlist(mclapply(used.gos, function(x) { 
        ukb.goa.hits$ECO[which(ukb.goa.hits$GO == x)]
    } ) ) )
    eco.used.gos <- eco.used.gos.all[which(eco.used.gos.all %in% ECO.OBO$id)]
    li.used <- table(unlist(mclapply(eco.used.gos, function(x) { evidenceCodeBins(x) } ) ) )
    
    eco.no.used.gos.all <- unique(unlist(mclapply(no.used.gos, function(x) { 
        ukb.goa.hits$ECO[which(ukb.goa.hits$GO == x)]
    } ) ) )
    eco.no.used.gos <- eco.no.used.gos.all[which(eco.no.used.gos.all %in% ECO.OBO$id)]
    li.no.used <- table(unlist(mclapply(eco.no.used.gos, function(x) { evidenceCodeBins(x) } ) ) )

    trusted.used[i] <- li.used[1]/(length(eco.no.used.gos)+length(eco.used.gos))
    untrusted.used[i] <- li.used[2]/(length(eco.no.used.gos)+length(eco.used.gos))
    trusted.no.used[i] <- li.no.used[1]/(length(eco.no.used.gos)+length(eco.used.gos))
    untrusted.no.used[i] <- li.no.used[2]/(length(eco.no.used.gos)+length(eco.used.gos))
    
}

#' Mean Relative Reference Protein Abundance (MRRPA)

#' Identify Bins sharing GOAs
z <- mm.2.go.df[which(mm.2.go.df$MapManBin.GO != ""), "MapManBin.GO"]
coords.duplicated.goas <- which(duplicated(z) | duplicated(z, fromLast = TRUE))
duplicated.goas <- mm.2.go.df$MapManBin.GO[coords.duplicated.goas]
# duplicated.bins <- mm.2.go.df$MapManBin[coords.duplicated.goas]
means.bins.rel.ref.abund <- c()
for (i in 1:length(duplicated.goas)) { 

  #' Identify Bins
    bins.1.goa <- mm.2.go.df$MapManBin[which(mm.2.go.df$MapManBin.GO == duplicated.goas[i])]

    for (j in 1:length(bins.1.goa)) {
      a <- paste("mm.2.go$`", paste(bins.1.goa[j], "`$MapManBin.GO", sep = ""), sep = "")
      o.go <- sort(unique(unlist(eval(parse(text = a)))))
      used.gos[j] <- strsplit(o.go, split = ",")
    }
    means.bins.rel.ref.abund[i] <- mean(table(unlist(used.gos)) / length(used.gos))
}

#' - Histogram of Entropies
pdf(file.path(input.args[[1]], "inst", "GoTermsNotUsedEntropyHist.pdf"))
plotDistAsHistAndBox(sh.entropy, "Shannon Entropy of not used GO Terms in GOAs per MapMan-Bin")
dev.off()

#' - Histogram of Mean Relative Reference Protein Abundance
pdf(file.path(input.args[[1]], "inst", "MeanRelativeReferenceProteinAbundanceHist.pdf"))
plotDistAsHistAndBox(means.bins.rel.ref.abund, "Mean Relative Reference Protein Abundance")
dev.off()

#' - Boxplot of GO Term Annotations trustworthiness
pdf(file.path(input.args[[1]], "inst", "GOAsTrustworthiness.pdf"))
boxplot(trusted.used, untrusted.used, trusted.no.used, untrusted.no.used, col=c("lightgrey", "lightgrey", "white", "white"), 
        main= "Trustworthiness of GO Terms Annotations", pch="-", names=c("Trusted", "Untrusted", "Trusted", "Untrusted"))
legend("topright", inset=.05, c("Used","No used"), fill=c("lightgrey", "white"))
dev.off()

save(sh,entropy, means.bins.rel.ref.abund, trusted.used, untrusted.used, trusted.no.used, untrusted.no.used, file = file.path(input.args[[1]], 
    "data", "MapManBinsRevision.RData"))

