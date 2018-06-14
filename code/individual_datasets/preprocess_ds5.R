
source('processing_utils.R')



# ds5 - GSE14323, GPL571
gse <- getGEOData("GSE14323")
gse14323p1 <- gse$originalData$GSE14323_GPL571 # ignore other platform
#gse14323p2 <- gse$originalData$GSE14323_GPL96
# Transcription profiling of human liver samples from subjects with hepatitis C virus (HCV), HCV-HCC, or normal liver
# normal livers negative for HCV, from: https://www.pathology.umn.edu/research/liver-tissue-cell-distribution-system - have to be careful if overlap
# Notes: - 19 normal samples 
# no sex labels
pData_ds5p1 <- gse14323p1$pheno
#pData_ds5p2 <- gse14323p2$pheno # no non-disease --> ignore!
pData_ds5.2 <- pData_ds5p1[pData_ds5p1$"Tissue:ch1"=="Normal",]
cel.names_ds5 <- stripGSM(pData_ds5.2$supplementary_file)
res_ds5 <- loadAffyData(cel.names_ds5, "GSE14323")
expData_ds5 <- res_ds5$exp
boxplot(expData_ds5)
qc_ds5 <- res_ds5$qc
plot(qc_ds5) # numbers, not percents are red, one line has red ratios

#expData_ds5 <- gse14323p1$expr[,rownames(pData_ds5.2)]
#sexLabels_ds5 <- labelSex(expData_ds5, gse14323p1$keys, ychr.genes, threshold=3) # 7 female, 12 male, for either threshold


ds5.keys <- loadGPL571EntrezAnnot()
sexLabels_ds5 <- labelSex(expData_ds5, ds5.keys, ychr.genesEntrez$entrezgene, 3, plot=TRUE) # looks good!

table(sexLabels_ds5$sex) 
# --> sanity check: 
sanityCheckSexLabels(expData_ds5, ds5.keys, sexLabels_ds5$sex) # looks like good separation to me --> gonna stick with it



# save
datasetInfo5 <- list("expr"=expData_ds5, "pheno"=sexLabels_ds5$sex, "keys"=ds5.keys, "ID"="GSE14323")
save(datasetInfo5, file="data/processed/ds5.RData")
