
# ds4 - GSE48452, GPL11532
source('processing_utils.R')

gse13 <- getGEOData("GSE48452")
gse48452 <- gse13$originalData$GSE48452
# Human liver biopsy of different phases from control to NASH
# 5m, 7f
# no reference, sex labels present

pData_ds13 <- gse48452$pheno
#"Sex:ch1"
#"age:ch1"
#"bmi:ch1"

# filter to get get normal-weight controls
pData_ds13.2 <- pData_ds13[pData_ds13$source_name_ch1=="Control",] # don't want patients after surgery
pData_ds13.2[,c("age:ch1",
                "bmi:ch1")] # all adult, normal bmi - though some "old" - note this as a caveat

#expData_ds13 <- gse48452$expr[,rownames(pData_ds13.2)]
#sexLabels_ds13 <- labelSex(expData_ds13, gse48452$keys, ychr.genes, threshold=4, plot=FALSE)  

cel.names_ds13 <- stripGSM(pData_ds13.2$supplementary_file)

#The affy package can process data from the Gene ST 1.x series of arrays,
#but you should consider using either the oligo or xps packages, which are specifically
#designed for these arrays
# --> use oligo

cel.names_ds13.2 <- sapply(cel.names_ds13, function(x) paste("data/affy/GSE48452", x, sep="/"))
rawData_ds13 <- read.celfiles(cel.names_ds13.2)
rawData_ds13 # there is no feature data

pmsLog13 <- log2(pm(rawData_ds13))
boxplot(pmsLog13) # looks fine
# Here, normalization is not necessary because the median of the samples is similar, and the data is already in log scale because the expression values are between 0 and 15. (If negative expression values would be observed e.g. the lowest expression value is -1, we recommend to shift all expression values of the dataset above 0 by adding +1 to each gene expression measurement in all samples.)
expData_ds13 <- rma(rawData_ds13) 

keys.vec13 <- loadGPL11532EntrezAnnot()
sexLabels_ds13 <- labelSex(exprs(expData_ds13), keys.vec13, ychr.genesEntrez$entrezgene, threshold=3, plot=FALSE) # 7f, 5m
(sexLabels_ds13$sex == pData_ds13.2$"Sex:ch1") # all match -- YAY!


datasetInfo13 <- list("expr"=exprs(expData_ds13), "pheno"=sexLabels_ds13$sex, "keys"=keys.vec13, "ID"="GSE48452")
save(datasetInfo13, file="data/processed/ds13.RData")
