message("USAGE: Rscript path/2/MapMan2GO/exec/plotDepths.R path/2/MapManBins2GO.RData path/2/MapMan2GO")

#' Input command line arguments:
input.args <- commandArgs(trailingOnly = TRUE)

load(input.args[[1]])

#' Calculate median, max and mean
meanBins.GO.depth <- vector()
for (i in 1:length(mm.2.go.df$MapManBin)){
  meanBins.GO.depth[i] <- mean(mm.2.full.desc$GO.depth[grep(pattern = mm.2.go.df$MapManBin[i], mm.2.full.desc$MapManBin, fixed =T)])
}
mm.2.go.df$mean.GO.depth <- meanBins.GO.depth

medianBins.GO.depth <- vector()
for (i in 1:length(mm.2.go.df$MapManBin)){
  medianBins.GO.depth[i] <- median(mm.2.full.desc$GO.depth[grep(pattern = mm.2.go.df$MapManBin[i], mm.2.full.desc$MapManBin, fixed =T)])
}
mm.2.go.df$median.GO.depth <- medianBins.GO.depth

maxBins.GO.depth <- vector()
for (i in 1:length(mm.2.go.df$MapManBin)){
  maxBins.GO.depth[i] <- max(mm.2.full.desc$GO.depth[grep(pattern = mm.2.go.df$MapManBin[i], mm.2.full.desc$MapManBin, fixed =T)])
}
mm.2.go.df$max.GO.depth <- maxBins.GO.depth

#' - Boxplot of GO.depths statistics
pdf(file.path(input.args[[2]], "inst", "GO.depthMapMan-BinsDistribution.pdf"))
  boxplot(mm.2.go.df$median.GO.depth, mm.2.go.df$max.GO.depth, mm.2.go.df$mean.GO.depth, ylab = "GO.depth", names=c("Median", "Maximum", "Mean"), col = "lightgrey", pch = '-')
dev.off()

pdf(file.path(input.args[[2]], "inst", "GO.depthGOtermsDistribution.pdf"))
#  plotDistAsHistAndBox(mm.2.full.desc$GO.depth, "GO.depth GO terms distribution")
  def.mar <- par('mar')
  m.1 <- def.mar; m.1[[1]] <- 0
  op <- par(mfcol = 2:1, mar=m.1)
  hist(mm.2.full.desc$GO.depth, col = "lightgrey", main = "GO.depth GO terms distribution", xlab = NULL)
  m.2 <- def.mar; m.2[[3]] <- 0
  par( mar=m.2 )
  boxplot(mm.2.full.desc$GO.depth, col = "lightgrey", horizontal = TRUE, frame = FALSE, pch = '|')
  par(op)
dev.off()

#' Scatter plot Number of genes vs GO depth
pdf(file.path(input.args[[2]], "inst", "NumberOfGenesVsGODepthMean.pdf"))
  plot(mm.2.go.df$n.genes, mm.2.go.df$mean.GO.depth, xlab="Number of genes per MapMan-Bin", ylab="GO.depth", pch=19)
#  abline(lm(mm.2.go.df$mean.GO.depth~mm.2.go.df$n.genes), col="red") # regression line (y~x)
dev.off()

pdf(file.path(input.args[[2]], "inst", "NumberOfGenesVsGODepthMedian.pdf"))
  plot(mm.2.go.df$n.genes, mm.2.go.df$median.GO.depth, xlab="Number of genes per MapMan-Bin", ylab="GO.depth", pch=19)
#  abline(lm(mm.2.go.df$median.GO.depth~mm.2.go.df$n.genes), col="red") # regression line (y~x)
dev.off()

pdf(file.path(input.args[[2]], "inst", "NumberOfGenesVsGODepthMaximum.pdf"))
  plot(mm.2.go.df$n.genes, mm.2.go.df$max.GO.depth, xlab="Number of genes per MapMan-Bin", ylab="GO.depth", pch=19)
#  abline(lm(mm.2.go.df$max.GO.depth~mm.2.go.df$n.genes), col="red") # regression line (y~x)
dev.off()
