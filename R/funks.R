#' Extracts short UniProtKB gene accessions from the long version identifiers,
#' e.g. 'sp|Q8GWW7|AGUA_ARATH' will become 'Q8GWW7'.
#'
#' @param ukb.accs A character vector of long UniProtKB gene accessions
#'
#' @export
#' @return The short gene accessions extracted from the argument
#' \code{ukb.accs}.
sanitizeAccession <- function(ukb.accs) {
    sub("^[^|]+\\|", "", sub("\\|[^|]+$", "", ukb.accs))
}

#' Extracts the unique and sorted Gene Ontology Term Annotations for a given
#' gene.
#'
#' @param gene.id The unique identifier of the gene to lookup GOA for
#' @param goa.tbl An instance of \code{data.frame} holding GOAs. Default is
#' \code{getOption('MapMan2GO.goa.tbl', ukb.goa.hits)}
#' @param gene.col The column of \code{goa.tbl} number or name in which to
#' lookup the genes. Default is \code{getOption('MapMan2GO.goa.tbl.gene.col',
#' 3)}.
#' @param go.col The column of \code{goa.tbl} number or name in which to lookup
#' the GO terms. Default is \code{getOption('MapMan2GO.goa.tbl.go.col',2)}
#' @param extend.goas.with.ancestors boolean indicating whether to extend each
#' proteins' GOA with the ancestors of the respective GO Terms. Default is
#' \code{getOption('MapMan2GO.extend.goa.with.ancestors', TRUE)}.
#'
#' @export
#' @return A character holding the GO terms for \code{gene.id}
compoundGoAnnotation <- function(gene.id, goa.tbl = getOption("MapMan2GO.goa.tbl", 
    ukb.goa.hits), gene.col = getOption("MapMan2GO.goa.tbl.gene.col", 3), go.col = getOption("MapMan2GO.goa.tbl.go.col", 
    2), extend.goas.with.ancestors = getOption("MapMan2GO.extend.goa.with.ancestors", 
    TRUE)) {
    res.goa <- unique(sort(goa.tbl[which(goa.tbl[, gene.col] == gene.id), go.col]))
    if (extend.goas.with.ancestors && length(res.goa) > 0) {
        addAncestors(res.goa)
    } else {
        res.goa
    }
}

#' Adds all GO Terms that are ancestral to the argument \code{go.terms}.
#'
#' @param go.terms a character vector of GO identifier, e.g.
#' \code{'GO:006969'}.
#' @param go.obo The result of reading the Gene Ontology OBO file with
#' \code{ontologyIndex::get_ontology('go.obo')}. Default is
#' \code{getOption('MapMan2GO.go.obo', GO.OBO)}. See \code{MapMan2GO} data for
#' more details.
#'
#' @return A character vector including the argument \code{go.terms} and all
#' found ancestors. Returns \code{character(0)} if the argument \code{go.terms}
#' is NULL, all of its entries are \code{NA}, or if its length equals zero.
#' @export
addAncestors <- function(go.terms, go.obo = getOption("MapMan2GO.go.obo", GO.OBO)) {
    if (is.null(go.terms) || all(is.na(go.terms)) || length(go.terms) == 0) {
        character(0)
    } else {
        sort(unique(unlist(lapply(go.terms, function(g.id) {
            if (g.id %in% go.obo$id) {
                go.obo$ancestors[[g.id]]
            } else {
                g.id
            }
        }))))
    }
}

#' Infers the compound GO annotations found for the genes related to the
#' argument MapMan Bin and measures the resulting Shannon Entropy of these
#' compound GO annotations.
#'
#' @param map.man.bin The identifier of the MapMan Bin to assign GO terms to.
#' @param mm.bins.vs.genes An instance of \code{data.frame} holding
#' MapManBin-Gene-Relations. Default is
#' \code{getOption('MapMan2GO.seq.sim.tbl', mm.bins.vs.sprot)}.
#' @param mm.bin.col The column of \code{mm.bins.vs.genes} in which to lookup
#' the MapMan-Bins. Default is \code{getOption('MapMan2GO.seq.sim.tbl.bin.col',
#' 'MapManBin')}.
#' @param mm.gene.col The column of \code{mm.bins.vs.genes} in which to lookup
#' the gene identifiers. Default is
#' \code{getOption('MapMan2GO.seq.sim.tbl.gene.col', 'Swissprot.Short.ID')}.
#' @param mm.bins.vs.genes An instance of \code{data.frame} with at least two
#' columns. It must hold mappings of \code{map.man.bin} to genes at least
#' partially found in \code{goa.tbl}.
#' @param goa.tbl An instance of \code{data.frame} holding GOAs. Default is
#' \code{getOption('MapMan2GO.goa.tbl', ukb.goa.hits)}
#' @param gene.col The column of \code{goa.tbl} number or name in which to
#' lookup the genes. Default is \code{getOption('MapMan2GO.goa.tbl.gene.col',
#' 3)}.
#' @param go.col The column of \code{goa.tbl} number or name in which to lookup
#' the GO terms. Default is \code{getOption('MapMan2GO.goa.tbl.go.col',2)}
#' @param extend.goas.with.ancestors boolean indicating whether to extend each
#' proteins' GOA with the ancestors of the respective GO Terms. Default is
#' \code{getOption('MapMan2GO.extend.goa.with.ancestors', TRUE)}.
#'
#' @export
#' @return An instance of \code{list} with the following named entries:
#' 'Shannon.Entropy' is the measured entropy of compound GO annotations
#' retrieved for the genes related to \code{map.man.bin}. The second entry
#' 'genes.goa' is a list with each genes' compound GO Annotations,
#' 'mutual.information' holds the mutual information between the Bin's GOA and
#' those found for the Bin's reference proteins, 'MapManBin.GO' is the
#' intersection of the genes' compound GO Annotations to be used as the
#' MapMan-Bin's compound GO Annotation, n.GO is the number of GO Terms in the
#' MapMan-Bin's compound GO Annotation, median.n.GO is the median of the number
#' of GO Terms in the genes' GOAs, and n.genes is the number of genes related
#' to \code{map.man.bin}.
compoundGoAnnotationEntropy <- function(map.man.bin, mm.bins.vs.genes = getOption("MapMan2GO.seq.sim.tbl", 
    mm.bins.vs.sprot), mm.bin.col = getOption("MapMan2GO.seq.sim.tbl.bin.col", 
    "MapManBin"), mm.gene.col = getOption("MapMan2GO.seq.sim.tbl.gene.col", "Swissprot.Short.ID"), 
    goa.tbl = getOption("MapMan2GO.goa.tbl", ukb.goa.hits), gene.col = getOption("MapMan2GO.goa.tbl.gene.col", 
        3), go.col = getOption("MapMan2GO.goa.tbl.go.col", 2), extend.goas.with.ancestors = getOption("MapMan2GO.extend.goa.with.ancestors", 
        TRUE)) {
    tryCatch({
        gene.ids <- mm.bins.vs.genes[which(mm.bins.vs.genes[, mm.bin.col] == map.man.bin), 
            mm.gene.col]
        genes.goa <- setNames(lapply(gene.ids, function(g.id) {
            compoundGoAnnotation(g.id, goa.tbl, gene.col, go.col)
        }), gene.ids)
        s.e <- shannonEntropy(table(as.character(unlist(lapply(genes.goa, paste, collapse = ",")))))
        m.i <- mutualInformationBinGoaGenesGoas(genes.goa)
        bin.goa <- Reduce(intersect, genes.goa[which(as.logical(lapply(genes.goa, 
            function(x) length(x) > 0 && !is.na(x) && !is.null(x))))])
        if (extend.goas.with.ancestors) {
            bin.goa <- addAncestors(bin.goa)
        }
        bin.goa <- sort(bin.goa)
        list(Shannon.Entropy = s.e, genes.goa = genes.goa, mutual.information = m.i, 
            MapManBin.GO = paste(bin.goa, collapse = ","), n.GO = length(bin.goa), 
            median.n.GO = median(unlist(lapply(genes.goa, length)), na.rm = TRUE), 
            n.genes = length(gene.ids))
    }, error = function(e) {
        message("MapMan-Bin '", e, "' caused an error:\n", e)
    })
}

#' Computes the empirical counts of events in a sample. Events are GO Term
#' annotations.
#'
#' @param go.sample.space A vector representing the sample space.
#' @param go.annos A vector of samples
#'
#' @export
#' @return A named numeric of empirical counts, names are the events in
#' \code{go.sample.space} and values are the empirical counts as found in
#' \code{go.annos}.
goCounts <- function(go.sample.space, go.annos) {
    tryCatch({
        setNames(as.numeric(lapply(go.sample.space, function(go.t) length(which(go.annos == 
            go.t)))), go.sample.space)
    }, error = function(e) {
        browser()
    })
}

#' Compute the mutual information between a MapMan Bin's GOA and the GOAs found
#' for it's reference genes. In this, regard each single GO Term as a event in
#' the sample space. Because each GO Term can only have binary count values in
#' the Bin's GOA the Bin's counts are normalized by multiplying with them with
#' the number of reference genes.
#'
#' @param genes.goa A list of character vectors. Names are reference gene
#' identifiers and values, the character vectors, represent each gene's GO Term
#' Annotations. See \code{compoundGoAnnotationEntropy} for more details.
#'
#' @export
#' @return A numeric value, the computed mutual information in bits (log2).
mutualInformationBinGoaGenesGoas <- function(genes.goa) {
    genes.gos <- Reduce(c, genes.goa)
    bin.gos <- Reduce(intersect, genes.goa)
    go.sample.space <- unique(genes.gos)
    genes.go.counts <- goCounts(go.sample.space, genes.gos)
    bin.go.counts.norm <- length(genes.goa) * goCounts(go.sample.space, bin.gos)
    mi.plugin(rbind(genes.go.counts, bin.go.counts.norm), unit = "log2")
}

#' Generates a two row plot with the first one being a Histogram and the second
#' row a horizontal Boxplot.
#'
#' @param x The values passed into \code{hist} and \code{boxplot}
#' @param The main title of the resulting plot
#' @param summary.as.title boolean indicating whether to add the output of
#' \code{base::summary(x)} to the title of the plot. Default is
#' \code{getOption('MapMan2GO.plot.dist.summary.as.title', FALSE)}.
#'
#' @export
#' @return TRUE if and only if no error has occurred
plotDistAsHistAndBox <- function(x, main, summary.as.title = getOption("MapMan2GO.plot.dist.summary.as.title", 
    FALSE)) {
    def.mar <- par("mar")
    m.1 <- def.mar
    m.1[[1]] <- 0
    op <- par(mfcol = 2:1, mar = m.1)
    if (summary.as.title) {
        main <- paste(c(main, capture.output(summary(x))), collapse = "\n")
    }
    hist(x, col = "lightgrey", main = main, xlab = NULL)
    m.2 <- def.mar
    m.2[[3]] <- 0
    par(mar = m.2)
    boxplot(x, col = "lightgrey", horizontal = TRUE, frame = FALSE, pch = "|")
    par(op)
    TRUE
}

#' Analyze the number of shared words among the MapMan-Bins descriptions and the
#' GOs that conform their GOA.
#'
#' @param mm.2.go.df.MapManBin.GO a list of the GOAs to analyze that conform the
#' MapMan-Bin
#' @param mm.2.go.df.MapManBin a list of the MapMan-Bins to analyze. Must have
#' the same length as \code{mm.2.go.df.MapManBin.GO}
#'
#' @return An instance of \code{data frame} with the following entries:
#' 'mm.2.go.df.MapManBin' is the MapMan-Bin analyzed.
#' 'n.words.shared' are the number of words shared between the MapMan-Bins
#' and their GOAs.
#' 'percent.shared.words' is the percentage of shared words getween the
#' MapMan-Bin and its corresponding GOA.
#' 'bins.gos.shared.words' is a list of words shared.
#' 'bins.gos.no.shared.words' is a list of the words from the MapMan-Bin
#' description that are not shared with the GOA.
#' @export
analyzeSharedWords <- function(mm.2.go.df.MapManBin.GO, mm.2.go.df.MapManBin) {
    n.words.shared <- c()
    bins.gos.shared.words <- c()
    bins.gos.no.shared.words <- c()
    percent.shared.words <- c()
    
    for (j in 1:length(mm.2.go.df$MapManBin.GO)) {
        splited.st <- strsplit(mm.2.go.df$MapManBin.GO[j], ",")
        
        GO.names <- c()
        for (i in 1:length(splited.st[[1]])) {
            GO.names[i] <- GO.OBO$name[splited.st[[1]][i]]
        }
        
        GO.description.splited <- unique(unlist(strsplit(GO.names, split = "[., ,_,:]")))
        bin.description.splited <- strsplit(gsub("\\(", "", gsub("\\)", "", gsub(",", 
            "", mm.desc.df$Description[which(mm.desc.df$MapManBin == mm.2.go.df$MapManBin[j])]))), 
            split = "[., ,_,:]")
        bin.description.splited <- if (!length(bin.description.splited)) {
            c("na")
        } else tolower(bin.description.splited[[1]])
        bin.description.splited <- bin.description.splited[which(bin.description.splited %in% 
            unlist(words.to.filter) == "FALSE")]
        
        words.grepl <- sapply(bin.description.splited, function(x) any(grepl(x, 
            GO.description.splited, ignore.case = TRUE)))
        n.words.shared[j] <- length(unique(tolower(names(which(words.grepl == "TRUE")))))
        bins.gos.shared.words[j] <- paste(unique(tolower(names(which(words.grepl == 
            "TRUE")))), collapse = ",")
        bins.gos.no.shared.words[j] <- paste(unique(tolower(names(which(words.grepl == 
            "FALSE")))), collapse = ",")
        percent.shared.words[j] <- n.words.shared[j]/(length(bin.description.splited) + 
            length(GO.description.splited))
        
    }
    info.bins.words <- data.frame(mm.2.go.df.MapManBin, n.words.shared, percent.shared.words, 
        bins.gos.shared.words, bins.gos.no.shared.words)
    return(info.bins.words)
}

#' Wraps \code{base::read.table} to parse the result produced by 'Mercator'.
#'
#' @param path.2.mercator.result.tbl The valid file path to the text file
#' holding the table produced as result by Mercator.
#' @param add.go.terms If set to \code{TRUE} the MapMan-Bin to Gene Ontology Term
#' mappings are looked up and automatically added as an additional column.
#' Default is \code{getOption('MapMan2GO.add.GO.terms', TRUE)}
#' @param map.man.bins.2.go An instance of \code{base::data.frame} holding the
#' MapMan-Bin to Gene Ontology Term mappings. It must have at least two
#' columns: One named 'MapManBin' holding the BINCODEs and another named
#' 'MapManBin.GO' holding the compound GO Term annotations separated by comma.
#' Default is \code{getOption('MapMan2GO.map.man.bins.2.go', mm.2.go.df)}.
#' @param sanitize.accession boolean if TRUE the mercator table column
#' 'IDENTIFIER' will be passed through \code{MapMan2GO::sanitizeAccession} and
#' an additional column 'IDENTIFIER.LONG' will hold the long Uniprot protein
#' identifiers. Default is \code{getOption(
#' 'MapMan2GO.read.mercator.sanitize.accession', FALSE)}.
#'
#' @return An instance of \code{base::data.frame} holding the results of the
#' mercator annotation pipeline.
#' @export
readMercatorResultTable <- function(path.2.mercator.result.tbl, add.go.terms = getOption("MapMan2GO.add.GO.terms", 
    TRUE), map.man.bins.2.go = getOption("MapMan2GO.map.man.bins.2.go", mm.2.go.df), 
    sanitize.accession = getOption("MapMan2GO.read.mercator.sanitize.accession", 
        FALSE)) {
    m.df <- read.table(path.2.mercator.result.tbl, header = TRUE, sep = "\t", colClasses = c(rep("character", 
        4), "logical"), quote = "", na.strings = "", comment.char = "", stringsAsFactors = FALSE)
    for (i.col in colnames(m.df)) {
        if (class(m.df[, i.col]) == "character") {
            m.df[, i.col] <- sub("^\\s*'", "", sub("'\\s*$", "", m.df[, i.col]))
        }
    }
    if (add.go.terms) {
        m.df$MapManBin.GO <- unlist(lapply(m.df$BINCODE, function(mm.bin) {
            i <- which(map.man.bins.2.go$MapManBin == mm.bin)
            if (length(i) > 0) {
                map.man.bins.2.go[[i, "MapManBin.GO"]]
            } else NA
        }))
    }
    if (sanitize.accession) {
        m.df$IDENTIFIER.LONG <- m.df$IDENTIFIER
        m.df$IDENTIFIER <- sanitizeAccession(m.df$IDENTIFIER)
    }
    
    m.df
}

#' Uses the MapMan-Bin to Gene Ontology mappings generated within the
#' \code{MapMan2GO} package to assign a set of Gene Ontology terms to the
#' argument MapMan-Bin \code{mm.bincode}. Optional recursion looks up GO Terms
#' assigned to a parent MapMan-Bin, if no GO Terms are found for the argument
#' MapMan-Bin.
#'
#' @param mm.bincode A string identifying the MapMan-Bin, eg. \code{1.1.1.2}.
#' @param map.man.2.go An instance of \code{base::data.frame} holding the
#' MapMan-Bin to Gene Ontology mappings. The mapped GO Terms are expected to be
#' concatenated into a comma separated string. The table must have at least the
#' following two columns: \code{MapManBin} with the MapMan-Bin's identifiers
#' (BINCODE), and \code{MapManBin.GO} with the comma separated and concatenated
#' GO Terms assigned to the respective Bin. Default is
#' \code{getOption('MapMan2GO.map.man.2.go', mm.2.go.df)}.
#' @param use.parental.bins A BOOLEAN indicating wether to perform recursive
#' lookup of GO Terms assigned to parent Bins of \code{mm.bincode}, if
#' \code{mm.bincode} has no GO Terms assigned in \code{map.man.2.go}. Default
#' is \code{getOption('MapMan2GO.use.parental.bins', FALSE)}.
#' @param parental.regex A string representing the regular expresion to be used
#' to delete the matching trailing substring yielding the parental BINCODE of
#' \code{mm.bincode}. Default is \code{getOption('MapMan2GO.parental.regex',
#' '\\.[0-9]+$')}.
#'
#' @return A string of comma separated concatenated GO Terms assigned to the
#' argument \code{mm.bincode} or one of its ancestors (see argument
#' \code{use.parental.bins}). \code{NA} is returned, if no GO Terms can be
#' found.
#' @export
goaForMapManBin <- function(mm.bincode, map.man.2.go = getOption("MapMan2GO.map.man.2.go", 
    mm.2.go.df), use.parental.bins = getOption("MapMan2GO.use.parental.bins", FALSE), 
    parental.regex = getOption("MapMan2GO.parental.regex", "\\.[0-9]+$")) {
    mm.goa <- NA
    if (!is.null(mm.bincode) && !is.na(mm.bincode) && mm.bincode != "") {
        if (mm.bincode %in% map.man.2.go$MapManBin) {
            mm.goa <- map.man.2.go[which(map.man.2.go$MapManBin == mm.bincode), 
                "MapManBin.GO"]
        } else if (use.parental.bins && grepl(parental.regex, mm.bincode)) {
            parent.bincode <- sub(parental.regex, "", mm.bincode)
            mm.goa <- goaForMapManBin(parent.bincode, map.man.2.go, use.parental.bins)
        }
    }
    mm.goa
}

#' Tests function \code{MapMan2GO::goaForMapManBin}.
#'
#' @return \code{TRUE} if and only if all tests pass successfully.
#' @export
testGoaForMapManBin <- function() {
    m2go <- data.frame(MapManBin = c("1", "1.1", "2"), MapManBin.GO = LETTERS[1:3])
    t.1 <- goaForMapManBin("1", map.man.2.go = m2go) == LETTERS[1]
    t.2 <- goaForMapManBin("1.1", map.man.2.go = m2go) == LETTERS[2]
    t.3 <- goaForMapManBin("1.1.1", map.man.2.go = m2go) == LETTERS[2]
    t.4 <- is.na(goaForMapManBin("1.1.1", map.man.2.go = m2go, use.parental.bins = FALSE))
    t.5 <- is.na(goaForMapManBin(NULL, map.man.2.go = m2go))
    t.6 <- is.na(goaForMapManBin(NA, map.man.2.go = m2go))
    t.7 <- is.na(goaForMapManBin("", map.man.2.go = m2go))
    t.8 <- is.na(goaForMapManBin(mm.bincode = NULL))
    t.9 <- goaForMapManBin("1.1.1.1.1.1.1", map.man.2.go = m2go) == LETTERS[2]
    all(c(t.1, t.2, t.3, t.4, t.5, t.6, t.7, t.8, t.9))
}

#' Identifies the subset of evidence codes ancestral to \code{eco.id} that
#' intersect with argument \code{ancestral.eco.ids}.
#'
#' @param eco.id The evidence code's identifier for which to find intersecting
#' ancestors.
#' @param ontology An instance of class \code{ontologyIndex::ontology_index}
#' generated by \code{ontologyIndex::get_ontolgy('eco.obo'}. Default is
#' \code{getOption('MapMan2GO.eco.ontology', ECO.OBO)}.
#' @param ancestral.eco.ids A character vector of evidence code IDs to be
#' intersected with the ancestors of argument \code{eco.id}. Default is
#' \code{getOption('MapMan2GO.ancestral.eco.ids', eco.first.level)}.
#'
#' @return A character vector holding those terms ancestral to \code{eco.id}
#' found also in \code{ancestral.eco.ids}. \code{character(0)} is returned if
#' none are found.
#' @export
intersectAncestralEvidenceCodes <- function(eco.id, ontology = getOption("MapMan2GO.eco.ontology", 
    ECO.OBO), ancestral.eco.ids = getOption("MapMan2GO.ancestral.eco.ids", eco.first.level)) {
    intersect(ancestral.eco.ids, ontology$ancestors[[eco.id]])
}

#' Classifies any given evidence code into the categories defined in argument
#' named list \code{ancestral.evidence.code.bins}. The names function of the
#' categories.
#'
#' @param eco.id The evidence code identifier to be classified.
#' @param ancestral.evidence.code.bins A named list where the names are the
#' categories and the values are ancestral evidence codes helping to classify
#' all their children into categories not present in the ontology. Default is
#' \code{getOption('MapMan2GO.ancestral.evidence.code.bins', list(TRUSTED =
#' c('ECO:0000006', 'ECO:0000088', 'ECO:0000212', 'ECO:0000352'), 
#' UNTRUSTED = c('ECO:0000041', 'ECO:0000177', 'ECO:0000204', 'ECO:0000311', 'ECO:0000361',
#' 'ECO:0000501', 'ECO:0006055')))}.
#' @param there.can.be.only.one boolean indicating whether to return only a
#' single group, even if more than one are assigned. In this case define group
#' priority by the order of groups appearing in argument
#' \code{ancestral.evidence.code.bins}, the first matching will be returned.
#' Default is \code{getOption( 'MapMan2GO.there.can.be.only.one', TRUE )}.  
#'
#' @export
#' @return A character vector the names of the categries into which the
#' argument \code{eco.id} could be put.
evidenceCodeBins <- function(eco.id, ancestral.evidence.code.bins = getOption("MapMan2GO.ancestral.evidence.code.bins", 
    list(TRUSTED = c("ECO:0000006", "ECO:0000088", "ECO:0000212", "ECO:0000352"), 
        UNTRUSTED = c("ECO:0000041", "ECO:0000177", "ECO:0000204", "ECO:0000311", 
            "ECO:0000361", "ECO:0000501", "ECO:0006055"))), there.can.be.only.one = getOption("MapMan2GO.there.can.be.only.one", 
    TRUE)) {
    eco.classes <- intersectAncestralEvidenceCodes(eco.id)
    bin.i <- as.logical(lapply(ancestral.evidence.code.bins, function(eco.bin.ids) {
        any(eco.classes %in% eco.bin.ids)
    }))
    eco.groups <- names(ancestral.evidence.code.bins)[bin.i]
    if (there.can.be.only.one) 
        eco.groups[[1]] else eco.groups
}

#' Enables appending to an existing RData file.
#'
#' @param ... the R objects to append
#' @param list same as in \code{base::save}
#' @param file A valid file path to an RData file. If it does not exist, it
#' will be created only containing arguments \code{...}.
#'
#' @return The result of invoking \code{base::save}.
#' @export
appendToRData <- function(..., list = character(), file) {
    previous <- load(file)
    var.names <- c(list, as.character(substitute(list(...)))[-1L])
    for (var in var.names) assign(var, get(var, envir = parent.frame()))
    save(list = unique(c(previous, var.names)), file = file)
}

#' Add alpha channel (transparency) to a color.
#'
#' @param col The color to add transparency to
#' @param alpha A value between zero and one, the higher the more transparency.
#' Default is \code{.25}.
#'
#' @return The modified color \code{col}.
#' @export
addAlpha <- function(col, alpha = 0.25) {
    apply(sapply(col, col2rgb)/255, 2, function(x) rgb(x[1], x[2], x[3], alpha = alpha))
}

#' Splits the concatenated Gene Ontology Term Annotations (GOA) assigned to a
#' MapManBin into a vector.
#'
#' @param map.man.bin.goas A characer vector of MapManBin-GOAs
#'
#' @return A character vector of sorted, unique GO Terms found in the comma
#' separated MapManBin-GOAs in argument \code{map.man.bin.goas}.
#' @export
splitMapManBinGOAs <- function(map.man.bin.goas) {
    sort(unique(unlist(strsplit(map.man.bin.goas, ","))))
}

#' Computes the normalized empirical Shannon Entropy for counts delivered in
#' the argument \code{counts.table}. Basis of the \code{log} function is
#' natural and normalization is done by division by the maximum entropy
#' \code{log(length(counts.table))} (see
#' \href{https://en.wikipedia.org/wiki/Entropy <- (information <- theory)#Efficiency)}{Normalized
#' Entropy}.
#'
#' @param counts.table result of invoking \code{base::table(count.vector)}
#' @param log.base The base of the logarithm. Default is the natural logarithm,
#' \code{getOption('MapMan2GO.entropy.log.base', base::exp(1))}.
#'
#' @export
#' @return A numeric value element [0,1], the normalized Shannon Entropy.
shannonEntropy <- function(counts.table, log.base = getOption("MapMan2GO.entropy.log.base", 
    base::exp(1))) {
    if (length(counts.table) <= 1) 
        return(0)
    c.t.s <- sum(counts.table)
    -sum(sapply(counts.table, function(x) x/c.t.s * log(x/c.t.s, base = log.base)))/log(length(counts.table), 
        base = log.base)
}

#' Testing function \code{MapMan2GO::shannonEntropy}
#'
#' @export
#' @return \code{TRUE} if and only if all tests are passed successfully.
testShannonEntropy <- function() {
    test.1 <- identical(shannonEntropy(table(c())), 0)
    test.2 <- identical(shannonEntropy(table(NULL)), 0)
    test.3 <- identical(shannonEntropy(table(NA)), 0)
    test.4 <- identical(shannonEntropy(table(1)), 0)
    test.5 <- identical(shannonEntropy(table(1:2)), 1)
    all(c(test.1, test.2, test.3, test.4, test.5))
}


