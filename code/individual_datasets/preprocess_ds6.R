
source("code/processing_utils.R")


# ds6 - 
gse6 <- getGEOData("GSE55668")
gse55668 <- gse6$originalData$GSE55668
# Human liver, normal adult female vs male
# part of superseries GSE55669 - others are mouse, cells
# 7 male, 7 female - none need to be removed
# sex labels provided, no reference info
pData_ds6 <- gse55668$pheno
annotSexLabels_ds6 <- ifelse(pData_ds6$characteristics_ch1=="Sex: Female", "female", "male")
#head(pData_ds6) # - Agilent

#head((gse55668$expr)[,1:10]) # there are a lot of zeros... hmm, will try looking at the raw data as well

expData_ds6 <- gse55668$expr

# map
keys.vec6 <- loadGPL4133EntrezAnnot()

# sex label
sexLabels_ds6 <- labelSex(expData_ds6, keys.vec6, ychr.genesEntrez$entrezgene, threshold=3)

# this failed --> by hand
xist.probe <- names(keys_ds6.2)[as.vector(keys_ds6.2)=="7503"] # XIST, ENSG00000229807
print(xist.probe)
rps4y1.probe <- names(keys_ds6.2)[as.vector(keys_ds6.2)=="6192"] # RPS4Y1, ENSG00000129824
print(rps4y1.probe)
plot(expData_ds6[xist.probe,], expData_ds6[rps4y1.probe,], ylab="RPS4Y1", xlab="XIST")   ## good separation
my.df6 <- data.frame(t(expData_ds6[c(xist.probe, rps4y1.probe),]))
colnames(my.df6) <- c("XIST", "RPS4Y1")
clus6 <- kmeans(my.df6, 2) # cluster 1 is female, cluster 2 is male - MAKE SURE!
male.cluster <- which.max(clus6$centers[,"RPS4Y1"])
sexLabels_ds6 <- ifelse(clus6$cluster==male.cluster, "male", "female")
table(sexLabels_ds6) # 7 and 7
(annotSexLabels_ds6==sexLabels_ds6) ## ALL MATCH


# save
keys.vec6.2 <- keys.vec6[names(keys.vec6) %in% rownames(expData_ds6)]

datasetInfo6 <- list("expr"=expData_ds6, "pheno"=sexLabels_ds6, "keys"=keys.vec6.2, "ID"="GSE55668")
save(datasetInfo6, file="data/processed/ds6.RData")

