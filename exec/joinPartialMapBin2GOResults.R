require(MapMan2GO)

message("USAGE: Rscript path/2/joinPartialMapBins2GOResults.R path/2/partial_result_files_dir/ path/2/MapMan2GO")

input.args <- commandArgs(trailingOnly = TRUE)

#' Define persistent objects not affected by the subsequent "loads":
#' mm.2.go.p, mm.2.go.df.p, mm.2.full.desc.p, mm.desc.df.p
mm.2.go.p <- c()
mm.2.go.df.p <- c()
mm.2.full.desc.p <- c()
# mm.desc.df.p <- c()

#' For each file in partial_result_files_dir
tmp.results <- list.files(path=input.args[[1]], pattern= ".RData", full.names=T, recursive=FALSE)

#' load(partial_result_file)
for(i in 1:length(tmp.results)) {
  load(tmp.results[i])
  mm.2.go.p <- append (mm.2.go.p, mm.2.go)
  mm.2.go.df.p <- rbind(mm.2.go.df.p, mm.2.go.df)
  mm.2.full.desc.p <- append (mm.2.full.desc.p, mm.2.full.desc)
  # mm.desc.df.p <- rbind (mm.desc.df.p, mm.desc.df)
} 

write.table(mm.2.full.desc.p, file.path(input.args[[1]], "MapManBins2GO.txt"), sep = "\t", row.names = FALSE)
save(mm.2.go.p, mm.2.go.df.p, mm.2.full.desc.p, mm.desc.df.p, file = file.path(input.args[[1]], "MapManBins2GO.Rdata"))
