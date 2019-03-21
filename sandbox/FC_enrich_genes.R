require('gage')
require('plyr')
require('miceadds')
require('metap')
require('humanFC')
load("external_data/fc_space/sig.genes.RData")


sig.genes.lst <- apply(sig.genes, 2, function(fc) {rownames(sig.genes)[which(fc==1)]})
names(sig.genes.lst) <- sapply(c(1:ncol(sig.genes)), function(x) paste("FC", as.character(x), sep=""))
# these are all entrez genes


# map data to symbols
require('biomaRt')
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

#sig.genes2 <- c(metaObj2$filterResults$FDR0.05_es0.4_nStudies2_looaFALSE_hetero0$posGeneNames,
#metaObj2$filterResults$FDR0.05_es0.4_nStudies2_looaFALSE_hetero0$negGeneNames)
mapping.table <- getBM(c('hgnc_symbol', 'entrezgene'), "entrezgene",rownames(sig.genes), mart=ensembl)

sig.genes.lst2 <- lapply(sig.genes.lst, function(x) mapping.table[mapping.table$entrezgene %in% x,]$hgnc_symbol) 
sig.genes.lst3 <- lapply(sig.genes.lst2, function(x) x[!is.na(x)])
save(sig.genes.lst2, file="external_data/fc_sig_genes_symbol.RData")

# compute an enrichment score using GSEA

computeEnrich <- function(ds.obj, gsc){
  enrich <- gage(ds.obj$exp, gsets=gsc, ref=ds.obj$ref, samp=ds.obj$samp, compare="as.group", saaTest=gs.KSTest)
  return(enrich)
}
combine.p.vals <- function(my.dat, gsc){
  ## TODO: deduplicate
  greater <- data.frame(do.call(rbind, lapply(names(gsc), function(x) {
    list.p <- sapply(my.dat, function(ds) ds$greater[x,"p.val"]); 
    list.p2 <- list.p[!is.na(list.p)];
    if (length(list.p2) < 2){return("NA")}
    else{ 
      list.ps <- unlist(list.p2)[order(unlist(list.p2))];
      rth.rank <- list.ps[round(0.3*length(list.ps))] # 7th
      return(rth.rank )
    }
  })))
  colnames(greater) <- "combined.p"
  rownames(greater) <- names(gsc)
  greater$label <- rownames(greater)
  greater$combined.p <- sapply(greater$combined.p, function(x) as.numeric(as.character(x)))
  greater <- greater[order(greater$combined.p),]
  less <- data.frame(do.call(rbind, lapply(names(gsc), function(x) {
    list.p <- sapply(my.dat, function(ds) ds$less[x,"p.val"]); 
    list.p2 <- list.p[!is.na(list.p)];
    if (length(list.p2) < 2){return("NA")}
    else{ 
      list.ps <- unlist(list.p2)[order(unlist(list.p2))];
      rth.rank <- list.ps[round(0.3*length(list.ps))] # 7th
      return(rth.rank )
    }
  })))
  colnames(less) <- "combined.p"
  rownames(less) <- names(gsc)
  less$label <- rownames(less)
  
  less$combined.p <- sapply(less$combined.p, function(x) as.numeric(as.character(x)))
  less <- less[order(less$combined.p),]
  return(list("greater"=greater, "less"=less))
}


loadReformatDat <- function(idx){
  print(idx)
  load.Rdata(sprintf("data/processed/ds%sg.RData", idx), "my.ds")
  
  expDat <-data.frame(my.ds$expr)
  expDat$probe_id <- rownames(expDat)
  keys2 <- my.ds$keys[rownames(expDat)]
  probeGene <- stack(keys2)
  colnames(probeGene) <- c("gene", "probe_id")
  expDatv2 <- merge(expDat, probeGene, by="probe_id")
  res2 <- aggregate(expDatv2, list(expDatv2$gene), mean)
  res2$probe_id <- NULL
  rownames(res2) <- res2$Group.1
  res2$Group.1 <- NULL
  # this will lead to some many-to-many - mult rownames mapping to same gene
  # here I am taking the average - may want to replace
  my.list <- list(exp=res2, ref=which(my.ds$pheno=="male"), samp=which(my.ds$pheno=="female"))
  return(my.list)
}
allDat0 <- lapply(c(1:7), loadReformatDat) # slowwwww
allDat <- lapply(c(8:16), loadReformatDat)
enriched0 <- lapply(allDat0, function(ds) computeEnrich(ds, sig.genes.lst2))
enriched <- lapply(allDat, function(ds) computeEnrich(ds, sig.genes.lst2))


fc.enriched <- c(enriched0, enriched)
combined.p <- data.frame(do.call(rbind, lapply(names(sig.genes.lst2), function(x) {
  list.p <- sapply(fc.enriched, function(ds) ds$greater[x,"p.val"]); 
  list.p2 <- list.p[!is.na(list.p)];
  print(length(list.p2)); 
  if (length(list.p2) < 2){return("NA")}
  else{    list.ps <- unlist(list.p2)[order(unlist(list.p2))];
  print(list.ps)
  rth.rank <- list.ps[round(0.3*length(list.ps))] # 7th
  print(rth.rank)
  return(rth.rank )}
})))
rownames(combined.p) <- names(sig.genes.lst2)
head(combined.p)
combined.p.up <- data.frame(combined.p)
colnames(combined.p.up) <- c("p")
combined.p.up$FC <-rownames(combined.p.up)
combined.p.o.up <- combined.p.up[order(combined.p.up[,1]),]


### annotate - what is going on in these FCs?
## seems too good to be true, bleh


combined.p.down <- data.frame(do.call(rbind, lapply(names(sig.genes.lst2), function(x) {
  list.p <- sapply(c(enriched0, enriched), function(ds) ds$less[x,"p.val"]); 
  list.p2 <- list.p[!is.na(list.p)];
  print(length(list.p2)); 
  if (length(list.p2) < 2){return("NA")}
  else{
    list.ps <- unlist(list.p2)[order(unlist(list.p2))];
    print(list.ps)
    rth.rank <- list.ps[round(0.3*length(list.ps))] # 7th
    print(rth.rank)
    return(rth.rank )
  }
})))

rownames(combined.p.down) <- names(sig.genes.lst2)
combined.p.down <- data.frame(combined.p.down)
colnames(combined.p.down) <- c("p")
combined.p.down$FC <-rownames(combined.p.down)
combined.p.o.down <- combined.p.down[order(combined.p.down[,1]),]
head(combined.p.o.down)

#combined.p.o.down <- combined.p.down[order(sapply(combined.p.down$p, unlist)),1:3]

#combined.p.o.down$FC <- rownames(combined.p.o.down)

goCodes <- sapply(1:139, getGO)
names(goCodes) <- names(sig.genes.lst2)
fcGO <- stack(goCodes)
colnames(fcGO) <- c("GOID", "FC")

require('GO.db')
resGO <- select(GO.db, keys=unlist(goCodes), columns=c("GOID", "TERM"))
fcGOannot <- join(resGO, fcGO, type="full")
head(aggregate(fcGOannot[,1:2], by=list(sapply(fcGOannot$FC, as.character)), function(x) toString(unique(x))))
head(combined.p.o.down)

# make sure plyr is not loaded!
#detach(package:plyr)
require('dplyr')
goAnnotCollapsed <- data.frame(fcGOannot %>% dplyr::group_by(FC) %>% dplyr::summarise(
  terms = paste(unique(TERM), collapse=" // "), 
  goids = paste(unique(GOID), collapse=" // ")))
head(goAnnotCollapsed)
p.down.annot <- merge(combined.p.o.down, goAnnotCollapsed, by="FC" )
p.up.annot <- merge(combined.p.o.up, goAnnotCollapsed, by="FC" )

pval.cut <- 0.01/(139*2)
p.down.annot$p <- sapply(p.down.annot$p, unlist)
p.up.annot$p <- sapply(p.up.annot$p, unlist)

p.down.annot.o <- p.down.annot[order(p.down.annot$p),]
p.up.annot.o <- p.up.annot[order(p.up.annot$p),]
immune.idx.down <- which(lapply(p.down.annot.o$terms, function(x) grep("immune", x, value=FALSE))==1)
immune.idx.up <- which(lapply(p.up.annot.o$terms, function(x) grep("immune", x, value=FALSE))==1)
immune.fcs <- unique(c(p.down.annot.o[immune.idx.down,]$FC, p.up.annot.o[immune.idx.up,]$FC))

## way too many :/ 
p.keep.down <- (p.down.annot.o[p.down.annot.o$p < pval.cut,])
p.keep.up <- (p.up.annot.o[p.up.annot.o$p < pval.cut,])
write.csv(p.keep.down, file="data/fc_pvals_down.csv", row.names=FALSE)
write.csv(p.keep.up, file="data/fc_pvals_up.csv", row.names=FALSE)

### STICK THESE IN A PLOT
study.tab.up <- do.call(rbind, lapply(c(enriched0, enriched), function(x) x$greater[,"p.val"]))
studies <- read.xlsx("data/annot/LiverMAStudies.xlsx", 1)
rownames(study.tab.up) <- studies$Accession[1:16]

study.tab.down <- do.call(rbind, lapply(c(enriched0, enriched), function(x) x$less[,"p.val"]))
rownames(study.tab.down) <- studies$Accession[1:16]

log10p.up <- apply(study.tab.up, c(1,2), function(x) -log10(as.numeric(x)))
log10p.down <- apply(study.tab.down, c(1,2), function(x) -log10(as.numeric(x)))
log10.up.df <- log10p.up[,c(p.keep.up$FC)]
log10.down.df <- log10p.down[,c(p.keep.down$FC)]
log10.down.df <- apply(log10.down.df, c(1,2), function(x) x*(-1))
logdf.combined <- cbind(log10.up.df, log10.down.df)
logdf.combined2 <- logdf.combined
logdf.combined2 <- apply(logdf.combined, c(1,2), function(x) ifelse(x > 8, 8, ifelse(x < -8, -8, x)))
immune.fcs <- c("FC38", "FC2", "FC91", "FC3", "FC8", "FC81", "FC38", "FC85", "FC19", "FC12")

logdf.combined3 <- logdf.combined2[, apply(logdf.combined2, 1, function(x)  any(x>1) | any(x < -1))]
immuneFCs <- sapply(colnames(logdf.combined3), function(x) ifelse(x %in% immune.fcs, "purple", "black"))
heatmap.2(t(logdf.combined3), trace="none",col=colorRampPalette(c("red", "white", "orange"))(n=299), margins=c(15,10), colRow=immuneFCs, cexRow=2, cexCol=2)

### look at genes involved
fcs_enriched <- unique(colnames(logdf.combined2)) # note - some are only up or down

save(logdf.combined2, file="data/liverfcResults1115.RData")

function(fc, type=)
fc38 <- metaObj$metaAnalysis$pooledResults[sig.genes.lst2$FC47,] # 345 genes!
fc38 <- fc38[!is.na(fc38$effectSize),]
fc38 <- fc38[order(fc38$effectSize),]
f <- fc38[fc38$effectSizeFDR < 0.2 & fc38$effectSize < 0,]
f[setdiff(rownames(f), fValid),]
m <- fc38[fc38$effectSize > 0 & fc38$effectSizeFDR < 0.2,]
m[setdiff(rownames(m), mValid),]


fc15 <- metaObj$metaAnalysis$pooledResults[sig.genes.lst2$FC15,] # 380 genes!
fc15 <- fc15[!is.na(fc15$effectSize),]
fc15 <- fc15[order(fc15$effectSize),]
f15 <- fc15[fc15$effectSizeFDR < 0.2,]
f15[setdiff(rownames(f15), fValid),]

fc47 <- metaObj$metaAnalysis$pooledResults[sig.genes.lst2$FC47,] # 380 genes!
fc47 <- fc47[!is.na(fc47$effectSize),]
fc47 <- fc47[order(fc47$effectSize),]
f47 <- fc47[fc47$effectSizeFDR < 0.2,]
f47[setdiff(rownames(f47), fValid),]


hist(sapply(sig.genes.lst2, length), xlab="number of genes", main="FC sizes (n=139)")

# I need some sort of background