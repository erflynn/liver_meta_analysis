# preprocess_ds1.R
#
# Notes: cannot sex label, use annotated

source('code/processing_utils.R')

gse <- getGEOData("GSE61276")
gse61276 <- gse$originalData$GSE

# filter based on phenotype
pTable <- gse61276$pheno
pTable2 <- pTable[pTable$`agegroup:ch1`=="Adult",]
unique(pTable2$'relation') # BE CAREFUL - it is used in another study - check for overlap


# look at annotated sex labels
sexLabels <- pTable2$`Sex:ch1`
table(sexLabels)
sexLabels <- sapply(sexLabels, function(x) tolower(as.character(x)))
names(sexLabels) <- pTable2$geo_accession

# plot exp data to check QC
boxplot(gse61276$expr[,rownames(pTable2)]) 
expData <- gse61276$expr[,rownames(pTable2)]


# map to entrez
keys.vec1 <- getIlluminaKeys("GSE61276", "GPL10558_HumanHT-12_V4_0_R2_15002873_B.txt")

# sex label - can't in this case - y chromosome probes are too poor!

# save - use annotated sex labels
datasetInfo1 <- list("expr"=expData, "pheno"=sexLabels, "keys"=keys.vec1, "ID"="GSE61267")
save(datasetInfo1, file="data/processed/ds1.RData")

