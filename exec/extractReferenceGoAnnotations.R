require(MapMan2GO)

message("Usage: Rscript ./exec/extractReferenceGoAnnotations.R path/2/UniProtKB_GOA_preprocessed_including_InterPro.txt path/2/list_of_gene_IDs.txt path/2/MapMan2GO")

input.args <- commandArgs(trailingOnly = TRUE)


ukb.goa <- fread(input.args[[1]], header = FALSE, sep = "\t", colClasses = rep("character", 
    4))
ref.gene.ids <- readLines(input.args[[2]])
ukb.ref.goas <- as.data.frame(ukb.goa[which(ukb.goa$V4 %in% ref.gene.ids), ])


#' Save results:
save(ukb.ref.goas, file = file.path(input.args[[3]], "data", "UKB_Reference_GOA_InterPro.RData"))


message("DONE")
