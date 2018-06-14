

source("processing_utils.R")

gse17 <- getGEOData("GSE23766")
gse23766 <- gse17$originalData$GSE23766
# Transcriptional Profiling of Human Liver Identifies Sex-Biased Genes Associated with Polygenic Dyslipidemia and Coronary Artery Disease
# study on sex-biased genes in liver, pooled, but sex divided???
### see paper for more about how to analyze


head(gse23766$expr) # these are the fold change (already) 
dim(gse23766$expr)
gpl6480 <- getGEO(filename="data/GPL6480.annot.gz")
# MAP

colnames(Table(gpl6480))
feat.labels <- data.frame(apply(Table(gpl6480)[,c("ID", "Gene ID")], c(1,2), as.character))
feat.labelssep <- separate_rows(feat.labels[,c(1,2)], "Gene.ID", sep="///")
keys.vec <- sapply(feat.labelssep$`Gene.ID`, as.character)
names(keys.vec) <- sapply(feat.labelssep$ID, as.character)
keys.vec17 <- keys.vec

expData17 <- gse23766$expr
save(keys.vec17, expData17, file="data/processed/ds17.RData")

head(keys.vec17)
meanDat <- rowMeans(gse23766$expr)
sdDat <- apply(gse23766$expr, 1, sd)
compositeScore <- sapply(1:nrow(gse23766$expr), function(i) {my.row <- gse23766$expr[i,]; multi <- (my.row*meanDat[i]); return(length(which(multi>0)))})

## compute Z scores
## based on: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1907322/
# z-ratio = g1 - mean(g1...gn) / SD(g1...gn)
#
# is it ok to do:
# [mean(g1) - mean(mean(g1)...mean(gn)) ]/ SD(mean(g1)..mean(gn))
mean.of.means <- mean(meanDat)
plot(density(meanDat)) # maybe normal? right center but v tight
sd.means <- sd(meanDat)
zs <- sapply(meanDat, function(x) (x-mean.of.means)/sd.means)
names(zs) <- names(meanDat)

# now convert zs to hedges' g
# see conversion table, then did algebra
# n <- ncol(expData17)
# ds <- sapply(zs, function(z) { 
#   num <- 4*n*(z**2);
#   print(num)
#   denom <- (1-n*(z**2));
#   print(denom)
#   return(sqrt(num/denom))
#   })
# 
# head(ds)
n <- 224
ds2 <- sapply(zs, function(z)
  { # this is the formula provided in the Conversions table
  # slightly different
  num <-z*sqrt(n);
  denom <- (1- sqrt((z**2)*n**(-1)));
  
  # flip sign - hopefully ok?
  val <- sqrt(abs(num/denom));
  negBool <- ifelse(num/denom < 0, -1, 1)
  return(negBool*val);
})
computeSEg <- function(g, n1, n2){
  sqrt( (n1+n2)/(n1*n2) + 0.5*g^2 /(n1+n2-3.94) )
}
se.g <- sapply(ds2, function(g) computeSEg(g, n1, n2))
summ.df <- data.frame(g=ds2, se.g=se.g, keys=keys.vec17[names(ds2)])

# how do I calculate the se.d???

# worried about already normalized

### alternatively, what if I used the t-statistic???
# calculate with limma, convert to hedges g


# need to calculate mean, sd
# they calculated the mean differently...
# composite array score, ranging from 8–16, was also determined based on the number of arrays out of 16 that match direction

### no sex labels... 
load('data/processed/ds17.RData')
summary(datasetInfo17$summ)
datasetInfo17 <- list("expr"=expData17, "pheno"="NA", "keys"=keys.vec17, "ID"="GSE23766", summ=summ.df)
save(datasetInfo17, file="data/processed/ds17.RData")

