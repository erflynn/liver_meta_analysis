# prev_liver_results.R
# LiverMetaProject
# E Flynn
# Updated: 4/16/2018
#
# Code for looking at previous results from Mayne, Gershoni-Pietrokovski, etc..

require('xlsx')
setwd("~/Documents/EMILY/Stanford/Coursework/stats211/project/")

gp <- read.csv("alt_data/gershoni_pietr_sde.csv")
# map to chromosome
gp$Ensembl <- sapply(gp$Gene, function(x) 
  strsplit(strsplit(as.character(x), "_", fixed=TRUE)[[1]][[2]], 
           ".", fixed=TRUE)[[1]][[1]])

# get a count of number / size of genes per chromosome

# map to chromosome
chr.data <- getBM(attributes=c("chromosome_name", "ensembl_gene_id"),
                    filters=c("ensembl_gene_id"),
                    values=list(gp$Ensembl), mart=ensembl)
ensembl.to.chr <- split(chr.data$chromosome_name, chr.data$ensembl_gene_id)
gp$chr <- sapply(gp$Ensembl, function(x) {x <- as.character(x); y <- ifelse(x %in% names(ensembl.to.chr), ensembl.to.chr[[x]], NA); return(y)})

tissue.chr.breakdownF <- sapply(2:46, function(x) {my.chr <- gp[gp[,x]>0,]; my.chr$chr <- factor(sapply(my.chr$chr, as.character), levels=c(sapply(1:22, as.character), "X", "Y")); table(my.chr$chr)})
colnames(tissue.chr.breakdownF) <- colnames(gp[,2:46])

tissue.chr.breakdownM <- sapply(2:46, function(x) {my.chr <- gp[gp[,x]<0,]; my.chr$chr <- factor(sapply(my.chr$chr, as.character), levels=c(sapply(1:22, as.character), "X", "Y")); table(my.chr$chr)})
colnames(tissue.chr.breakdownM) <- colnames(gp[,2:46])

test.res <- do.call(rbind, lapply(c(1:24), function(my.row) {
  do.call(cbind, lapply(1:45, function(idx) {
  x <- tissue.chr.breakdownF[my.row,idx]; 
  y <- tissue.chr.breakdownM[my.row, idx];
  if(((x+y) < 10) ){return("NA")}
  single.res <- chisq.test(c("F"=x, "M"=y))$p.value; 
  return(single.res)}
  ))}))

colnames(test.res) <- colnames(gp[,2:46])
test.res
num.nas <- sum(sapply(1:45, function(x) sum(test.res[,x]=="NA")))
pval.cutoff <- 0.05/(45*24 - num.nas)
arr.ind <- which(apply(test.res, c(1,2), as.numeric) < pval.cutoff, arr.ind=T)

rat.tab <- tissue.chr.breakdownF/tissue.chr.breakdownM
my.list <- rat.tab[arr.ind[arr.ind[,1]=="19",]]
names(my.list) <- colnames(test.res)[arr.ind[arr.ind[,1]=="19","col"]]


## look into chromosome sizes
# chromosome 19 contains a larger number of genes!
# - Zhang et al. - what they found with chr 19 is a strong sex-bias
chr.counts <- table(chr.data$chromosome_name)[c(1:22, "X", "Y")] # number of genes per chromosome
apply(apply(tissue.chr.breakdownF, 2, function(x) unlist(x)/unlist(chr.counts))[1:22,], 2, which.max)

mayneACC <- read.xlsx2("alt_data/mayne_meta_results.xslx", sheetName="Anterior Cingulate Cortex")
table(sapply(mayneACC[mayneACC$Chromosome=="19" & sapply(mayneACC$adj.p.value, function(z) as.numeric(as.character(z)) <0.005),]$summary, function(x) as.numeric(as.character(x))) < 0) 
chisq.test(c("F"=7, "M"=17))

mayne <- read.xlsx2("alt_data/mayne_meta_results.xslx", sheetName="Liver")

gp.liver <- gp[,c("Gene", "Liver")]
gp.liver$Ensembl <- sapply(gp$Gene, function(x) 
  strsplit(strsplit(as.character(x), "_", fixed=TRUE)[[1]][[2]], 
           ".", fixed=TRUE)[[1]][[1]])
gp.liver$Symbol <- sapply(gp$Gene, function(x) strsplit(as.character(x), "_", fixed=TRUE)[[1]][[1]])
head(gp.liver)
head(mayne)
gp.liver.sig <- gp.liver[gp.liver$Liver!=0,]
nrow(gp.liver.sig) # 25
nrow(mayne) # 32
mayne$Symbol <- sapply(mayne$Symbol, as.character)
mayne$Ensembl <- sapply(mayne$Ensembl, as.character)

length(intersect(mayne$Symbol, gp.liver.sig$Symbol)) # 7
length(intersect(mayne$Ensembl, gp.liver.sig$Ensembl)) # 7
mayne.tab.filt <- mayne[mayne$Symbol %in% gp.liver.sig$Symbol,]
table(mayne.tab.filt$Chromosome %in% c("X", "Y")) # 5 are X or Y 

gp.liver.sig[gp.liver.sig$Symbol %in% mayne$Symbol,] # 2 up in males (negative, TMSB4Y, USP9Y), rest up in females
# same pattern in Mayne et al.
intersect(mayne$Symbol, gp.liver.sig$Symbol)


