source("processing_utils.R")


gse28893obj <- getGEOData("GSE28893") 
# more of the same as previous dataset - these were used for the eQTL analysis
gse28893 <- gse28893obj$originalData$GSE28893
pData_ds16 <- gse28893$pheno
pData_ds16.2 <- pData_ds16[(pData_ds16$"age:ch1" > 17),] # removes 4

table(pData_ds16.2$"gender:ch1") # F = 26, M = 30
expData_ds16 <- gse28893$expr[,pData_ds16.2$geo_accession]
boxplot(expData_ds16)
annotLabels_ds16 <- sapply(pData_ds16.2$"gender:ch1", function(x) ifelse(x=="F", "female", "male"))
# for many it says "this sample consists of two replicates" - the GSMs appear to be averaged? --> might want to go back to raw data

# map keys
keys_vec16 <- getIlluminaKeys("GSE28893","GPL6104_Illumina_HumanRef-8_V2_0_R1_11223162_A.bgx.gz")

# sex label
sexLabels_ds16 <- labelSex(expData_ds16, keys_vec16, ychr.genesEntrez$entrezgene, threshold=3)
(annotLabels_ds16==sexLabels_ds16$sex) # all match!

# save
datasetInfo16 <- list("expr"=expData_ds16, "pheno"=sexLabels_ds16$sex, "keys"=keys_vec16, "ID"="GSE28893")
save(datasetInfo16, file="data/processed/ds16.RData")
