source("processing_utils.R")

gse12 <- getGEOData("GSE33814")
gse33814 <- gse12$originalData$GSE33814
# Gene expression profiling unravels cancer-related hepatic molecular signatures in steatohepatitis but not in steatosis
# expected: 13 normal
# no sex labels, no reference
pData_ds12 <- gse33814$pheno
pData_ds12.2 <- pData_ds12[pData_ds12$`diagnosis:ch1`=="normal",]
dim(pData_ds12.2)
expData_ds12 <- gse33814$expr[,rownames(pData_ds12.2)]
boxplot(expData_ds12)

# map to entrez

keys.vec_12 <- getIlluminaKeys("GSE33814", "GPL6884_HumanWG-6_V3_0_R0_11282955_A.bgx.gz")

# sex label
sexLabels_ds12 <- labelSex(expData_ds12, keys.vec_12, ychr.genesEntrez$entrezgene, threshold=3, plot=FALSE)  # 7f, 6m
sanityCheckSexLabels(expData_ds12, keys.vec_12, sexLabels_ds12$sex) # looks good

# save
datasetInfo12 <- list("expr"=expData_ds12, "pheno"=sexLabels_ds12$sex, "keys"=keys.vec_12, "ID"="GSE33814")
save(datasetInfo12, file="data/processed/ds12.RData")
