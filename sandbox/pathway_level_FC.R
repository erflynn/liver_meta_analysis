# pathway_level_FC.R
# E Flynn
# 4/16/2018
#
# Pathway-level analysis of sex-diff genes.
# Choices that I've made to think about:
#

require('gage')
require('plyr')
require('miceadds')
require('metap')
require('humanFC')
load("data/fc_space/sig.genes.RData")
data(kegg.gs)


# convert to a list
sig.genes.lst <- apply(sig.genes, 2, function(fc) {rownames(sig.genes)[which(fc==1)]})
names(sig.genes.lst) <- sapply(c(1:ncol(sig.genes)), function(x) paste("FC", as.character(x), sep=""))

loadReformatDat <- function(idx){
  print(idx)
  load.Rdata(sprintf("data/processed/ds%s.RData", idx), "my.ds")
  
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
allDat <- lapply(c(8:16), loadReformatDat)
allDat0 <- lapply(c(1:7), loadReformatDat)
save(allDat, allDat0, file="bcats_fc_input.RData")

load("bcats_fc_input.RData")
# run enrichment
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


kegg.enrichment <- combine.p.vals(c(enriched0.k, enriched.k), kegg.gs)
head(kegg.enrichment$greater)
pval.cut <- 0.05/(length(kegg.gs))

keep.greater.kegg <- rownames(kegg.enrichment$greater)[kegg.enrichment$greater$combined.p < pval.cut]
keep.less.kegg <- rownames(kegg.enrichment$less)[kegg.enrichment$less$combined.p < pval.cut]
keep.greater.kegg <- keep.greater.kegg[!is.na(keep.greater.kegg)]
keep.less.kegg <- keep.less.kegg[!is.na(keep.less.kegg)]


### MAKE A HEATMAP
study.tab.up <- do.call(rbind, lapply(c(enriched0.k, enriched.k), function(x) x$greater[,"p.val"]))
rownames(study.tab.up) <- studies$Accession[1:16]
  
study.tab.down <- do.call(rbind, lapply(c(enriched0.k, enriched.k), function(x) x$less[,"p.val"]))
rownames(study.tab.down) <- studies$Accession[1:16]



log10p.up <- apply(study.tab.up, c(1,2), function(x) -log10(as.numeric(x)))
log10p.down <- apply(study.tab.down, c(1,2), function(x) -log10(as.numeric(x)))
log10.up.df <- log10p.up[,c(keep.greater.kegg)]
log10.down.df <- log10p.down[,c(keep.less.kegg)]
log10.down.df <- apply(log10.down.df, c(1,2), function(x) x*(-1))
logdf.combined <- cbind(log10.up.df, log10.down.df)
logdf.combined2 <- logdf.combined
logdf.combined2 <- apply(logdf.combined, c(1,2), function(x) ifelse(x > 8, 8, ifelse(x < -8, -8, x)))
heatmap.2(t(logdf.combined2), trace="none", margins=c(9,20))
greater.p <- kegg.enrichment$greater$combined.p
names(greater.p) <- rownames(kegg.enrichment$greater)
p.log.greater <- sapply(greater.p[c(keep.greater.kegg)], function(x) -log10(as.numeric(as.character(x))))

less.p <- kegg.enrichment$less$combined.p
names(less.p) <- rownames(kegg.enrichment$less)

p.log.less <- sapply(less.p[c(keep.less.kegg)], function(x) log10(as.numeric(as.character(x))))
pooled.est <- data.frame(t(c(p.log.greater, p.log.less)))
colnames(pooled.est) <- c(keep.greater.kegg, keep.less.kegg)
logdf.combined3 <- rbind(logdf.combined2, pooled.est)
rownames(logdf.combined3)[17] <- "combined"

# enriched in females vs. males = greater
# enriched in males vs. females = less

heatmap.2(t(logdf.combined3[,c(keep.greater.kegg, keep.less.kegg)]), Rowv=FALSE, Colv=FALSE, dendrogram="none", trace="none", margins=c(9,30))

makeLogTab(t(greater.p), t(less.p),c(keep.greater.kegg, keep.less.kegg))

enriched0 <- lapply(allDat0, function(ds) computeEnrich(ds, sig.genes.lst))
enriched <- lapply(allDat, function(ds) computeEnrich(ds, sig.genes.lst))

head(enriched0.k[[1]]$greater)
enriched0.k <- lapply(allDat0, function(ds) computeEnrich(ds, kegg.gs))
enriched.k <- lapply(allDat, function(ds) computeEnrich(ds, kegg.gs))
combined.p.k <- data.frame(do.call(rbind, lapply(names(kegg.gs), function(x) {
  list.p <- sapply(c(enriched0.k, enriched.k), function(ds) ds$less[x,"p.val"]); 
  list.p2 <- list.p[!is.na(list.p)];
  print(length(list.p2)); 
  if (length(list.p2) < 2){return("NA")}
  else{return( sumlog(unlist(list.p2)) )}
})))
combined.p.k$pathway <- names(kegg.gs)
combined.p.k$p <- sapply(combined.p.k$p, unlist)
head(combined.p.k[order(combined.p.k$p),c(1:3, 5)])

combined.p.k.up <- data.frame(do.call(rbind, lapply(names(kegg.gs), function(x) {
  list.p <- sapply(c(enriched0.k, enriched.k), function(ds) ds$greater[x,"p.val"]); 
  list.p2 <- list.p[!is.na(list.p)];
  print(length(list.p2)); 
  if (length(list.p2) < 2){return("NA")}
  else{return( sumlog(unlist(list.p2)) )}
})))
combined.p.k.up$pathway <- names(kegg.gs)
combined.p.k.up$p <- sapply(combined.p.k.up$p, unlist)
head(combined.p.k.up[order(combined.p.k.up$p),c(1:3, 5)])
log10p <- sapply(combined.p.k.up$p, function(x) -log10(as.numeric(x)))
names(log10p) <- names(kegg.gs)
log10p.o <- log10p[order(-log10p)]
pcut.kegg <--log10(0.01/(2*length(kegg.gs)))
log10p.o2 <- log10p.o[!is.na(log10p.o)]
log10p.o2[unname(log10p.o2)>pcut.kegg]

fc.enriched <- c(enriched0, enriched)
combined.p <- data.frame(do.call(rbind, lapply(names(sig.genes.lst), function(x) {
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
rownames(combined.p) <- names(sig.genes.lst)
head(combined.p)
combined.p.up <- data.frame(combined.p)
colnames(combined.p.up) <- c("p")
combined.p.up$FC <-rownames(combined.p.up)
combined.p.o.up <- combined.p.up[order(combined.p.up[,1]),]


### annotate - what is going on in these FCs?
## seems too good to be true, bleh


combined.p.down <- data.frame(do.call(rbind, lapply(names(sig.genes.lst), function(x) {
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

rownames(combined.p.down) <- names(sig.genes.lst)
combined.p.down <- data.frame(combined.p.down)
colnames(combined.p.down) <- c("p")
combined.p.down$FC <-rownames(combined.p.down)
combined.p.o.down <- combined.p.down[order(combined.p.down[,1]),]
head(combined.p.o.down)

#combined.p.o.down <- combined.p.down[order(sapply(combined.p.down$p, unlist)),1:3]

#combined.p.o.down$FC <- rownames(combined.p.o.down)

goCodes <- sapply(1:139, getGO)
names(goCodes) <- names(sig.genes.lst)
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
View(p.keep.up)
p.keep.down <- (p.down.annot.o[p.down.annot.o$p < pval.cut,])
p.keep.up <- (p.up.annot.o[p.up.annot.o$p < pval.cut,])
write.csv(p.keep.down, file="fc_pvals_down.csv", row.names=FALSE)
write.csv(p.keep.up, file="fc_pvals_up.csv", row.names=FALSE)

### STICK THESE IN A PLOT
study.tab.up <- do.call(rbind, lapply(c(enriched0, enriched), function(x) x$greater[,"p.val"]))
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
immuneFCs <- sapply(colnames(logdf.combined2), function(x) ifelse(x %in% immune.fcs, "purple", "black"))
heatmap.2(t(logdf.combined2), trace="none", margins=c(9,20), colRow=immuneFCs )
greater.p <- kegg.enrichment$greater$combined.p
names(greater.p) <- rownames(kegg.enrichment$greater)
p.log.greater <- sapply(greater.p[c(keep.greater.kegg)], function(x) -log10(as.numeric(as.character(x))))

less.p <- kegg.enrichment$less$combined.p
names(less.p) <- rownames(kegg.enrichment$less)

p.log.less <- sapply(less.p[c(keep.less.kegg)], function(x) log10(as.numeric(as.character(x))))
pooled.est <- data.frame(t(c(p.log.greater, p.log.less)))
colnames(pooled.est) <- c(keep.greater.kegg, keep.less.kegg)
logdf.combined3 <- rbind(logdf.combined2, pooled.est)
rownames(logdf.combined3)[17] <- "combined"

# enriched in females vs. males = greater
# enriched in males vs. females = less

heatmap.2(t(logdf.combined3[,c(keep.greater.kegg, keep.less.kegg)]), Rowv=FALSE, Colv=FALSE, dendrogram="none", trace="none", margins=c(9,30))


rownames(enriched[[2]]$greater[enriched[[2]]$greater[,"q.val"] < 0.05/(139*2),])


### annotate - what is going on in these FCs?

#### leftover junk
combined.p <- data.frame(do.call(rbind, lapply(names(sig.genes.lst), function(x) sumlog(c(fc.p_ds4$greater[x,], fc.p_ds2$greater[x,]))))[,1:3])



fc.p_ds2 <- gage(myDat2$exp, gsets=sig.genes.lst, ref=myDat2$ref, samp=myDat2$samp, compare="as.group", saaTest=gs.KSTest)
head(fc.p_ds2$greater)
fc.p_ds2 <- gage(exp2, gsets=sig.genes.lst, ref=males, samp=females, compare="as.group", saaTest=gs.KSTest)

exp4 <- datasetInfo4$expr
rownames(exp4) <- sapply(rownames(exp4), function(x) datasetInfo4$keys[[x]])
fc.p_ds4 <- gage(exp4, gsets=sig.genes.lst, ref=which(datasetInfo4$pheno=="male"), samp=which(datasetInfo4$pheno=="female"), compare="as.group", saaTest=gs.KSTest)
head(fc.p_ds4$greater)
combined.p <- data.frame(do.call(rbind, lapply(names(sig.genes.lst), function(x) sumlog(c(fc.p_ds4$greater[x,], fc.p_ds2$greater[x,]))))[,1:3])
head(combined.p[order(sapply(combined.p$p, unlist)),])
