
source("processing_utils.R")

# - ds11

gse11 <- getGEOData("GSE33116")
gse33116 <- gse11$originalData$GSE33116 
# Breast cancer and liver tissue biopsies used to develop and validate an index for liver contamination in metastatic breast cancer
# no paper listed
# no sex labels, no reference
pData_ds11 <- gse33116$pheno
# filter for liver tissue (rest is breast or liver/breast)
pData_ds11.2 <- pData_ds11[pData_ds11$"tissue:ch1"=="Liver biopsy",] # 31

# load and process the raw data
cel.names_ds11 <- stripGSM(pData_ds11.2$supplementary_file)
res_ds11 <- loadAffyData(cel.names_ds11, "GSE33116")
expData_ds11 <- res_ds11$exp
qc_ds11 <- res_ds11$qc
plot(qc_ds11) # red for ratios, percentages,  almost all lines are blue

# not sure what to do
boxplot(expData_ds11) # looks ok?

ds11.keys <- loadGPL96EntrezAnnot()

#expData_ds11 <- gse33116$expr[,rownames(pData_ds11.2)]
sexLabels_ds11 <- labelSex(expData_ds11, ds11.keys, ychr.genesEntrez$entrezgene, threshold=3)  # 18 female, 13 male - looks great separation wise

datasetInfo11 <- list("expr"=expData_ds11, "pheno"=sexLabels_ds11$sex, "keys"=ds11.keys, "ID"="GSE33116")
save(datasetInfo11, file="data/processed/ds11.RData")

