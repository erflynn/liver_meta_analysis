# preprocess_ds3.R
#
# Notes:
#  - GSE38941, GPL570
#  channel_count = 1 --> that's good
#  TODO: check relation - GSEM366070

source('code/processing_utils.R')

# load phenotype data, filter for normal
gse38941 <- getGEOData("GSE38941")
pData_ds3 <- gse38941$originalData$GSE38941$pheno
pData_ds3.2 <- pData_ds3[pData_ds3$`disease state:ch1` == "normal",]

## loading raw data
cel.names_ds3 <- stripGSM(pData_ds3.2$supplementary_file)
res_ds3 <- loadAffyData(cel.names_ds3, "GSE38941")
expData_ds3 <- res_ds3$exp
qc_ds3 <- res_ds3$qc
boxplot(expData_ds3) # looks fine

# map to entrez
ds3.keys <- loadEntrezAnnot("GPL570") # TODO: change so not hard-coded

# label sex
sexLabels_table3 <- labelSex(expData_ds3, ds3.keys, ychr.genesEntrez$entrezgene, 3, plot=FALSE)
sexLabels_ds3 <- sexLabels_table3$sex
names(sexLabels_ds3) <- sexLabels_table3$ID

# save
datasetInfo3 <- list("expr"=expData_ds3, "pheno"=sexLabels_ds3, "keys"=ds3.keys, "ID"="GSE38941")
save(datasetInfo3, file="data/processed/ds3.RData")
