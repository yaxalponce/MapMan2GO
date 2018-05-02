message("USAGE: Rscript path/2/MapMan2GO/exec/mapBinsToGOs.R path/2/MapMan2GO")

require(MapMan2GO)
input.args <- commandArgs(trailingOnly = TRUE)

#' Identify duplicated GOAs and their MapManBins 
duplicated.GOAs <- mm.2.go.df$MapManBin.GO[ which(duplicated(mm.2.go.df$MapManBin.GO))]
coords.duplicated.GOAs <- which(duplicated(mm.2.go.df$MapManBin.GO))
duplicated.bins <- mm.2.go.df$MapManBin[coords.duplicated.GOAs]

#' Level at which Bins are sharing GOAs
#' Apply just for duplicated GOAs in different MapMan-Bins

duplicated.bins.combn <- combn(duplicated.bins, 2)
names(duplicated.GOAs) <- duplicated.bins

n <- length(duplicated.bins.combn)/2

MRCA.D <- c()
MRCA.Dist <- c()
for (j in 1:n) {
  if ( (duplicated.GOAs[names(duplicated.GOAs) == duplicated.bins.combn[1,j]]) == (duplicated.GOAs[names(duplicated.GOAs) == duplicated.bins.combn[2,j]]) ) {
      a <- unlist(strsplit(duplicated.bins.combn[1,j], ".", fixed=TRUE))
      b <- unlist(strsplit(duplicated.bins.combn[2,j], ".", fixed=TRUE))

      COUNT=1
      e <- min( length(a), length(b) )
      while(a[COUNT] == b[COUNT] && COUNT <= e) {
        COUNT = COUNT+1
      } 
      d <- max(length(a), length(b))
      MRCA.D[j] <- COUNT-1
      MRCA.Dist[j] <- d-(COUNT-1)
    } 
  }
# duplicated.bins.combn[,which(MRCA.D != "")]

#' MapManBins that share GOAs
#' Analysis on shared words

n.words.shared <- c()
bins.gos.shared.words <- c()
for (j in 1:length(duplicated.GOAs)) {
 splited.st <- strsplit(duplicated.GOAs[j], ",")
 duplicated.bin.description <- unique( mm.2.full.desc$Bin.Description[which(mm.2.full.desc$MapManBin == duplicated.bins[j])])

  GO.names <- c()
  for (i in 1:length(splited.st[[1]])) {
    GO.names[i] <-GO.OBO$name[splited.st[[1]][i]]
  }

  GO.description.splited <- unique(unlist(strsplit(GO.names, split="[., ,_,:]")))
  bin.description.splited <- strsplit (sub("\\(", "", sub("\\)", "", sub(",", "", duplicated.bin.description))), split="[., ,_,:]")
  
  words.grepl <- sapply(bin.description.splited[[1]], function(x) any(grepl(x, GO.description.splited, ignore.case=TRUE)))
  n.words.shared[j] <- length(unique(tolower(names(which(words.grepl == "TRUE")))))
  bins.gos.shared.words[j] <- paste(unique(tolower(names(which(words.grepl == "TRUE")))), collapse = ",")

}

info.bins.words <- data.frame(duplicated.bins, n.words.shared, bins.gos.shared.words)

#' Analysis on shared words for all MapManBins

n.words.shared.all <- c()
bins.gos.shared.words.all <- c()
bins.gos.no.shared.words.all <- c()
percent.shared.words <- c()

info.bins.words.all <- analyzeSharedWords(mm.2.go.df$MapManBin.GO, mm.2.go.df$MapManBin)

for (j in 1:length(mm.2.go.df$MapManBin.GO)) {
 splited.st.all <- strsplit(mm.2.go.df$MapManBin.GO[j], ",")

  GO.names.all <- c()
  for (i in 1:length(splited.st.all[[1]])) {
    GO.names.all[i] <-GO.OBO$name[splited.st.all[[1]][i]]
  }

  GO.description.splited.all <- unique(unlist(strsplit(GO.names.all, split="[., ,_,:]")))
  bin.description.splited.all <- strsplit (sub("\\(", "", sub("\\)", "", sub(",", "", mm.desc.df$Description[ which(mm.desc.df$MapManBin == mm.2.go.df$MapManBin[j])] ))), split="[., ,_,:]")
bin.description.splited.all <- tolower(bin.description.splited.all[[1]])

  words.grepl.all <- sapply(bin.description.splited.all, function(x) any(grepl(x, GO.description.splited.all, ignore.case=TRUE)))
  n.words.shared.all[j] <- length(unique(tolower(names(which(words.grepl.all == "TRUE")))))
  bins.gos.shared.words.all[j] <- paste(unique(tolower(names(which(words.grepl.all == "TRUE")))), collapse = ",")
  bins.gos.no.shared.words.all[j] <- paste(unique(tolower(names(which(words.grepl.all == "FALSE")))), collapse = ",")
  percent.shared.words[j] <- n.words.shared.all[j]/(length(bin.description.splited.all)+length(GO.description.splited.all))

}
info.bins.words.all <- data.frame(mm.2.go.df$MapManBin, n.words.shared.all, percent.shared.words, bins.gos.shared.words.all, bins.gos.no.shared.words.all)

#pdf("/LUSTRE/Genetica/ahallab/projects/MapMan2GO/MapMan2GO/inst/MapManBin.GOAs.Shared.Words.pdf")
pdf(file.path(input.args[[1]], "inst", "MapManBin.GOAs.Shared.Words.pdf"))
  plotDistAsHistAndBox(info.bins.words.all$percent.shared.words, "Percent of shared words among MapMan-Bins and GOAs")
dev.off()

######################

#' Add GOs to MapMan-Bin GOAs that have identical GOAs in it add new information

GOs.to.add.list <- c()
bins.GOs.to.add <- data.frame()

for (i in 1:length(duplicated.bins)) {

#' Union of protein GOAs
  a <- paste("mm.2.go$`", paste(duplicated.bins[i], "`$genes.goa", sep=""), sep="")

  GO.joined <- sort( unique( unlist( eval( parse( text=a)))))
  n <- which(info.bins.words.all$mm.2.go.df.MapManBin == duplicated.bins[i])
  GOA.2.compare <- strsplit(mm.2.go.df$MapManBin.GO[n], ",")

##' Working to include just words no shared before
  bin.description.2.compare <- strsplit( as.character(info.bins.words.all$bins.gos.no.shared.words[n]), split=",")

#' Identify GOs not present in original GOA
  GOs.to.eval <- GO.joined[ which(GO.joined %in% GOA.2.compare[[1]] == FALSE) ]

#' Split GOs name, search them in MapManBin description, if there's any, select the GO
  GO.names.extended <- c()
  GO.description.splited.extended <- c()
  GOs.to.add <- c()
  for (j in 1:length(GOs.to.eval)) {
    GO.names.extended[j] <- GO.OBO$name[GOs.to.eval[j]]
    GO.description.splited.extended <- unique( unlist( strsplit( GO.names.extended[j], split="[., ,_,:]")))
    if ( any(sapply(GO.description.splited.extended, function(x) x %in% bin.description.2.compare[[1]]) == TRUE) == TRUE) {
      GOs.to.add[j] <- GOs.to.eval[j]
     } else {
       GOs.to.add[j] <- NA
    }
  }
  if (!is.null(GOs.to.add) || !is.na(GOs.to.add)) { 
    GOs.to.add.list[i] <- paste( sort( unique( unlist(GOs.to.add))), collapse=",")
  } else {
    GOs.to.add.list[i] <- NA
  }
}
bins.GOs.to.add <- c()
bins.GOs.to.add <- data.frame(duplicated.bins, GOs.to.add.list)

#' Retrieve the Original GOAs to add the new GOs

bins.to.extend <- bins.GOs.to.add$duplicated.bins[which(bins.GOs.to.add$GOs.to.add.list != "")]

new.goa <- c()
n.new.GOs <- c()
for (i in 1:length(bins.to.extend)){
  goao <- mm.2.go.df$MapManBin.GO[which(mm.2.go.df$MapManBin == bins.to.extend[i])]
  goe <-  bins.GOs.to.add$GOs.to.add.list[which(duplicated.bins == bins.to.extend[i])]
  ng <- union(goao, goe)
  ng <- paste(ng, collapse=",")
  ng <- sort( unlist( strsplit(ng, split=",")))
  n.new.GOs[i] <- length(ng)
  new.goa[i] <- paste(ng, collapse=",")
}

#' Substitute original GOAs  and n.GOs with extended GOAs in data frame
# mm.2.go.df.copy <- cbind(mm.2.go.df)

for (i in 1:length(bins.to.extend)) {
  coords <- which(mm.2.go.df.copy$MapManBin == bins.to.extend[i])
  mm.2.go.df.copy$MapManBin.GO[coords] <- new.goa[i]
  mm.2.go.df.copy$n.GO[coords] <- n.new.GOs[i]
}

######################

#' Union of protein GOAs
shared.GOAs.extended <- c()
for (i in 1:length(duplicated.bins)) {
  a <- paste("mm.2.go$`", paste(duplicated.bins[i], "`$genes.goa", sep=""), sep="")
shared.GOAs.extended[i] <- paste( sort( unique( unlist (eval( parse( text=a))))), collapse=",")
}
names(shared.GOAs.extended) <- duplicated.bins
# duplicated.GOAs.extended <- shared.goas.extended[which(duplicated(shared.goas.extended))]
# duplicated.bins.extended <- names(duplicated.GOAs.extended)

#' Bins extended
n.words.shared.extended <- c()
bins.gos.shared.words.extended <- c()
percent.shared.words <- c()
for (j in 1:length(shared.GOAs.extended)) {
 splited.st.extended <- strsplit(shared.GOAs.extended[j], ",")

  GO.names.extended <- c()
  for (i in 1:length(splited.st.extended[[1]])) {
    GO.names.extended[i] <-GO.OBO$name[splited.st.extended[[1]][i]]
  }

  GO.description.splited.extended <- unique(unlist(strsplit(GO.names.extended, split="[., ,_,:]")))
  bin.description.splited.extended <- strsplit (sub("\\(", "", sub("\\)", "", sub(",", "", mm.desc.df$Description[ which(mm.desc.df$MapManBin == mm.2.go.df$MapManBin[j])] ))), split="[., ,_,:]")
bin.description.splited.extended <- tolower(bin.description.splited.extended[[1]])

  words.grepl.extended <- sapply(bin.description.splited.extended, function(x) any(grepl(x, GO.description.splited.extended, ignore.case=TRUE)))
  n.words.shared.extended[j] <- length(unique(tolower(names(which(words.grepl.extended == "TRUE")))))
  bins.gos.shared.words.extended[j] <- paste(unique(tolower(names(which(words.grepl.extended == "TRUE")))), collapse = ",")
  percent.shared.words[j] <- n.words.shared.extended[j]/(length(bin.description.splited.extended)+length(GO.description.splited.extended))

}
info.bins.words.extended <- data.frame(duplicated.bins, n.words.shared.extended, percent.shared.words, bins.gos.shared.words.extended)

b1 <- strsplit(mm.2.go.df$MapManBin.GO[which(mm.2.go.df$MapManBin == duplicated.bins[1])], ",")
b2 <- strsplit(shared.GOAs.extended[1], ",")
different.GOs <- setdiff(unlist(b1), unlist(b2))
different.GOs <- append(different.GOs,setdiff(unlist(b2), unlist(b1)))

GO.names.diff <- list()
# GO.description.splited.diff
for (i in length(different.GOs)) {
  GO.names.diff[i] <-GO.OBO$name[different.GOs[i]]
  while (GO.names.diff[i] == "") { next }
    GO.names.diff[i] <- unique(unlist(strsplit(GO.names.diff[i], split="[., ,_,:]")))
#info.bins.words.all[which(mm.2.go.df$MapManBin == duplicated.bins[1])]

}




words.diff %in% GO.extended
