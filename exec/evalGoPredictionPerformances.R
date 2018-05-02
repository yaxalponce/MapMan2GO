require(MapMan2GO)
options( mc.cores=detectCores() )

message("USAGE: Rscript path/2/MapMan2GO/exec/evalGoPredictionPerformances.R path/2/mercator_result.txt path/2/MapMan2GO")

input.args <- commandArgs(trailingOnly = TRUE)


mercator.annos <- readMercatorResultTable(input.args[[1]])
mercator.annos$F1.score <- mercatorF1Score(mercator.annos)


message("DONE")
