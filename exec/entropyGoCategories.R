require(MapMan2GO)
options(mc.cores = detectCores())

message("USAGE: Rscript /path/2/MapMan2GO/entropyGoCategories.R /path/2/MapMan2GO ")

input.args <- commandArgs(trailOnly = TRUE)

#' Get all GO terms from the three GO categories

mf.i <- which(grepl("molecular_function", GO.OBO$name, ignore.case = T))
bp.i <- which(grepl("biological_process", GO.OBO$name, ignore.case = T))
cc.i <- which(grepl("cellular_component", GO.OBO$name, ignore.case = T))

mf.des <- get_descendants(GO.OBO, GO.OBO$id[[mf.i]])
bp.des <- get_descendants(GO.OBO, GO.OBO$id[[bp.i]])
cc.des <- get_descendants(GO.OBO, GO.OBO$id[[cc.i]])


#' Normalized Shannon Entropy for used GO terms
mm.2.go.categories <- Reduce(rbind, mclapply(mm.2.go.df$MapManBin.GO, function(x) {
    y <- unlist(strsplit(x, split = ","))
    mf.int <- intersect(mf.des, y)
    bp.int <- intersect(bp.des, y)
    cc.int <- intersect(cc.des, y)
    sh.entropy.mf <- shannonEntropy(table(as.character(mf.int)))
    sh.entropy.bp <- shannonEntropy(table(as.character(bp.int)))
    sh.entropy.cc <- shannonEntropy(table(as.character(cc.int)))
    mf.int <- paste(mf.int, collapse = ",")
    bp.int <- paste(bp.int, collapse = ",")
    cc.int <- paste(cc.int, collapse = ",")
    data.frame(Mol.Fun.GO = mf.int, Biol.Proc.GO = bp.int, Cell.Comp.GO = cc.int, 
        Sh.Ent.MF = sh.entropy.mf, Sh.Ent.BP = sh.entropy.bp, Sh.Ent.CC = sh.entropy.cc)
}))
MapManBin <- mm.2.go.df$MapManBin
mm.2.go.categories <- cbind(mm.2.go.categories, MapManBin)

#' Normalized Shannon Entropy for compound GO annotations
compoundShannonEntropy.categories <- Reduce(rbind, mclapply(mm.2.go.df$MapManBin, 
    function(x) {
        gene.ids <- mm.bins.vs.sprot$Swissprot.Short.ID[which(mm.bins.vs.sprot$MapManBin == 
            x)]
        genes.goa <- setNames(lapply(gene.ids, function(g.id) {
            compoundGoAnnotation(g.id, ukb.goa.hits, 3, 2)
        }), gene.ids)
        mf.int <- intersect(mf.des, unlist(genes.goa))
        bp.int <- intersect(bp.des, unlist(genes.goa))
        cc.int <- intersect(cc.des, unlist(genes.goa))
        
        s.e.mf <- shannonEntropy(table(as.character(mf.int)))
        s.e.bp <- shannonEntropy(table(as.character(bp.int)))
        s.e.cc <- shannonEntropy(table(as.character(cc.int)))
        
        mf.int <- paste(mf.int, collapse = ",")
        bp.int <- paste(bp.int, collapse = ",")
        cc.int <- paste(cc.int, collapse = ",")
        
        data.frame(Mol.Fun.GO = mf.int, Biol.Proc.GO = bp.int, Cell.Comp.GO = cc.int, 
            Sh.Ent.MF = s.e.mf, Sh.Ent.BP = s.e.bp, Sh.Ent.CC = s.e.cc)
    }))
compoundShannonEntropy.categories  <- cbind(compoundShannonEntropy.categories, MapManBin)

#' Identify not used GO terms and calculate their normzalized Shannon Entropy
mm.2.go.categories.not.used.normalized <- Reduce(rbind, mclapply(mm.2.go.df$MapManBin, 
    function(x) {
        a <- paste("mm.2.go$`", paste(x, "`$MapManBin.GO", sep = ""), sep = "")
        o.go <- sort(unique(unlist(eval(parse(text = a)))))
        used.gos <- unlist(strsplit(o.go, split = ","))
        
        b <- paste("mm.2.go$`", paste(x, "`$genes.goa", sep = ""), sep = "")
        GO.joined <- sort(unique(unlist(eval(parse(text = b)))))
        not.used.gos <- setdiff(GO.joined, used.gos)
        
        mf.int.not.used <- intersect(mf.des, not.used.gos)
        bp.int.not.used <- intersect(bp.des, not.used.gos)
        cc.int.not.used <- intersect(cc.des, not.used.gos)
        sh.entropy.mf.not.used <- shannonEntropy(table(as.character(mf.int.not.used)))
        sh.entropy.bp.not.used <- shannonEntropy(table(as.character(bp.int.not.used)))
        sh.entropy.cc.not.used <- shannonEntropy(table(as.character(cc.int.not.used)))
        mf.int.not.used <- paste(mf.int.not.used, collapse = ",")
        bp.int.not.used <- paste(bp.int.not.used, collapse = ",")
        cc.int.not.used <- paste(cc.int.not.used, collapse = ",")
        data.frame(Mol.Fun.GO = mf.int.not.used, Biol.Proc.GO = bp.int.not.used, 
            Cel.Comp.GO = cc.int.not.used, Sh.Ent.MF = sh.entropy.mf.not.used, 
            Sh.Ent.BP = sh.entropy.bp.not.used, Sh.Ent.CC = sh.entropy.cc.not.used)
    }))
mm.2.go.categories.not.used.normalized <- cbind(mm.2.go.categories.not.used.normalized, 
    MapManBin)

#' Plot histograms of MF, BP and CC Gene Ontology categories
entropies <- c("Sh.Ent.MF", "Sh.Ent.BP", "Sh.Ent.CC")

#' - Histogram of Entropies from used GO terms
for (entropies.i in entropies) {
    pdf(file.path(input.args[[length(input.args)]], "inst", paste("used_", entropies.i, 
        "_Hist.pdf", sep = "")))
    plotDistAsHistAndBox(mm.2.go.categories[, entropies.i], main = paste("Used GO Terms: ", 
        entropies.i), summary.as.title = TRUE)
    dev.off()
}

#' - Histogram of Entropies for compound GO annotations
for (entropies.i in entropies) {
    pdf(file.path(input.args[[length(input.args)]], "inst", paste("used_", entropies.i, 
        "_CompoundGoaHist.pdf", sep = "")))
    plotDistAsHistAndBox(compoundShannonEntropy.categories[, entropies.i], main = paste("Compound GO annotation: ", 
        entropies.i), summary.as.title = TRUE)
    dev.off()
}

#' - Histogram of Entropies from not used GO terms
for (entropies.i in entropies) {
    pdf(file.path(input.args[[length(input.args)]], "inst", paste("notUsed_", entropies.i, 
        "_Hist.pdf", sep = "")))
    plotDistAsHistAndBox(mm.2.go.categories.not.used.normalized[, entropies.i], main = paste("Not used GO Terms: ", 
        entropies.i), summary.as.title = TRUE)
    dev.off()
}

#' Save results
save(mm.2.go.categories, mm.2.go.categories.not.used.normalized, compoundShannonEntropy.categories, file = file.path(input.args[[length(input.args)]], 
    "data", "MapManBinGoCategories.RData"))


