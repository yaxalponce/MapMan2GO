require(MapMan2GO)

message("USAGE: Rscript path/2/MapMan2GO/exec/loadSeqSimResults.R path/2/preProcessedSeqSimSearchResults.txt path/2/UniProtKB_GOA_preprocessed.txt path/2/MapMan2GO")

#' Input command line arguments:
input.args <- commandArgs(trailingOnly = TRUE)

#' Sequence Similarity Search Results:
mm.bins.vs.sprot <- fread(input.args[[1]], data.table = FALSE)
colnames(mm.bins.vs.sprot) <- c("MapManBin", "Swissprot.Hit")
mm.bins.vs.sprot$Swissprot.Short.ID <- sanitizeAccession(mm.bins.vs.sprot$Swissprot.Hit)


#' UniProtKB Gene Identifier to Gene Ontology Term Annotations (GOA):
ukb.goa <- fread(input.args[[2]], header = FALSE, sep = "\t", colClasses = rep("character", 
    3))
ukb.goa.hits <- as.data.frame(ukb.goa[which(ukb.goa$V3 %in% mm.bins.vs.sprot$Swissprot.Short.ID), 
    ], stringsAsFactors = FALSE)
colnames(ukb.goa.hits) <- c("ECO", "GO", "Swissprot.Hit")


#' Save results:
save(mm.bins.vs.sprot, ukb.goa.hits, file = file.path(input.args[[3]], "data", "MapManBinsVsSwissprotAndGOA.RData"))


message("DONE")
