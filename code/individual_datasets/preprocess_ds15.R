

source("processing_utils.R")
gse25935obj <- getGEOData("GSE25935")
# Genetic identification, replication, and functional fine-mapping of expression quantitative trait loci in primary human liver tissue

gse25935 <- gse25935obj$originalData$GSE25935
pData_ds15 <- gse25935$pheno
# channel count = 1
head(pData_ds15) # 464... hmm this is more than I expected

summary(sapply(pData_ds15$"age:ch1", as.numeric ))
table(pData_ds15$"age:ch1" < 18)
pData_ds15.2 <- pData_ds15[pData_ds15$"age:ch1" > 18,]
pData_ds15.2rep1 <- pData_ds15.2[sapply(pData_ds15.2$title, function(x) strsplit(as.character(x), "_", fixed=TRUE)[[1]][[2]]=="rep1"),]
pData_ds15.2rep2 <- pData_ds15.2[sapply(pData_ds15.2$title, function(x) strsplit(as.character(x), "_", fixed=TRUE)[[1]][[2]]=="rep2"),]
levels(as.factor(sapply(pData_ds15.2$title, function(x) strsplit(as.character(x), "_", fixed=TRUE)[[1]][[2]])))
table(pData_ds15.2rep1$"gender:ch1") # 68f, 115m - 183 total, this makes sense (it says 206)

boxplot(gse25935$expr) ## too big, need to view bit by bit

# - collapse replicates - #
## will have to decide what to do abt reps...
# oh there are two reps per sample
#   either we can: 1) select only one per sample (randomly or just rep1)
#              or: 2) meta-analyze/average over reps... hmm
# samples have anywhere from 1 to 5 reps!!!
# collapse replicates - take the average

pData_ds15.names <- pData_ds15.2[,c("title", "geo_accession", "gender:ch1")]
pData_ds15.names2 <- separate(pData_ds15.names, title, c("sampleID", "replicate"), sep="_")
head(pData_ds15.names2)

## AVERAGE BASED ON THIS
sample.to.reps <- split(pData_ds15.names2$geo_accession, pData_ds15.names2$sampleID)
my.sample <- sample.to.reps[[1]]
sampleAvg <- (rowMeans(gse25935$expr[,my.sample], na.rm=TRUE))

calcSampleAvg <- function(y){
  if(length(y)==1) {
    return(gse25935$expr[,y])
  }
return(rowMeans(gse25935$expr[,y], na.rm=TRUE))
}
avgSample <- lapply(1:length(sample.to.reps), function(x) { y <- sample.to.reps[[x]]; print(y); calcSampleAvg(y)})
#print(y); print(length(y)); ifelse((length(y)==1), gse25935$expr[,y] ,rowMeans(gse25935$expr[,y], na.rm=TRUE))})
avgSampleDf <- do.call(cbind, avgSample)
colnames(avgSampleDf) <- names(sample.to.reps)
head(avgSampleDf[,1:5])
expData_ds15 <- avgSampleDf

#expData_ds15 <- gse25935$expr[,pData_ds15.2rep1$geo_accession] # has lots of NAs
boxplot(expData_ds15) # goes above 16 a little bit

# map to ENTREZ
keys.vec15 <- loadGPL4133EntrezAnnot()

# sex label
sexLabels_ds15 <- labelSex(expData_ds15, keys.vec15, ychr.genesEntrez$entrezgene, threshold=3) #, plot=FALSE)  
## this doesn't work...

sexLabels_ds15 <- simpleSexLabel(expData_ds15, keys.vec15)

# xist.probe <- names(keys.vec15)[as.vector(keys.vec15)=="7503"] # XIST, ENSG00000229807
# print(xist.probe)
# rps4y1.probe <- names(keys.vec15)[as.vector(keys.vec15)=="6192"] # RPS4Y1, ENSG00000129824
# print(rps4y1.probe)
# 
# plot(expData_ds15[xist.probe,], expData_ds15[rps4y1.probe,],  ylab="RPS4Y1", xlab="XIST")  # good separation
# my.df15 <- data.frame(t(expData_ds15[c(xist.probe, rps4y1.probe),]))
# colnames(my.df15) <- c("XIST", "RPS4Y1")
# my.df15.2 <- my.df15[complete.cases(my.df15),]  # none removed
# clus15 <- kmeans(my.df15.2, 2) # cluster 1 is male, cluster 2 is female - MAKE SURE!
# male.cluster <- which.max(clus15$centers[,"RPS4Y1"])
# sexLabels_ds15 <- ifelse(clus15$cluster==male.cluster, "male", "female")
# table(sexLabels_ds15)
# ds15imputed <- data.frame(sexLabels_ds15)
# ds15imputed.2 <- cbind(ds15imputed, names(sexLabels_ds15))
# colnames(ds15imputed.2) <- c("sex", "sampleID")
# df15labels <- pData_ds15.names2[,c("sampleID", "gender:ch1")]
# compareTab15 <- join(df15labels[!duplicated(df15labels),], ds15imputed.2) 
sanityCheckSexLabels(expData_ds15, keys.vec15, sexLabels_ds15)

mislabeled15 <- compareTab15[(compareTab15$sex!=compareTab15$`gender:ch1`),] # one doesn't match! ruh roh

plot(expData_ds15[xist.probe,], expData_ds15[rps4y1.probe,], 
     col=ifelse(colnames(expData_ds15)==mislabeled15$sampleID, "red", "black"), ylab="RPS4Y1", xlab="XIST") 

# this is convincing - labeled it male, it appears to be male

## look at these

# save
datasetInfo15 <- list("expr"=expData_ds15, "pheno"=sexLabels_ds15, "keys"=keys.vec15, "ID"="GSE25935")
save(datasetInfo15, file="data/processed/ds15.RData")


