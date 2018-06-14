source("processing_utils.R")

gse10 <- getGEOData("GSE32504")
gse32504 <- gse10$originalData$GSE32504
# Identification of expression quantitative trait loci (eQTL) in human liver
# Expectation: 149 samples! 
# patients were undergoing liver surgery

pData_ds10 <- gse32504$pheno

# other demographic info present - smoking status, alcohol, medication, ethnic background
# 109/149 on medications, 29 smokers, 95 drink
# if look at the supplementary table - some of them have cholestasis (24) and a few NAs -  may want to exclude these?
summary(apply(gse32504$pheno[,c("alcohol intake:ch1", "ethnic background:ch1", "gender:ch1", "presurgical medication:ch1", "smoking:ch1")], 2, as.factor))

sapply(gse32504$pheno[,"age:ch1"], function(x) as.numeric(as.character(x))) # one is young --> take out
gse32504$pheno[sapply(gse32504$pheno[,"age:ch1"], as.numeric) < 18, ] # two --> remove these
pData_ds10.2 <- pData_ds10[sapply(gse32504$pheno[,"age:ch1"], as.numeric) > 18,]
expData_ds10 <- gse32504$expr[,rownames(pData_ds10.2)]

# map to entrez
keys.vec10 <- getIlluminaKeys("GSE32504", "GPL13376_HumanWG-6_V2_0_R4_11223189_A.txt.gz")


# sex label
sexLabels_ds10 <- labelSex(expData_ds10, keys.vec10, ychr.genesEntrez$entrezgene, threshold=3)  # 77 female, 70 male  - good separation

# RUH ROH - exact opposite
sexLabels_ds10$sex
pData_ds10.2$"gender:ch1"
sanityCheckSexLabels(expData_ds10, keys.vec10, sexLabels_ds10$sex) ### looks great
# well prediction concurs with what I'd expect - one low XIST, one high xist and opp with rsp4y1
## they really messed up labeling!

# save
datasetInfo10 <- list("expr"=expData_ds10, "pheno"=sexLabels_ds10$sex, "keys"=keys.vec10, "ID"="GSE32504")
save(datasetInfo10, file="data/processed/ds10.RData")
