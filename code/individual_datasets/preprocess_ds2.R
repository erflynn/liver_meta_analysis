# preprocess_ds2.R
#
# Notes:
#  - no sex labels
#  - we want sample A (before cold perfusion) - not sample B (after portal reprofusion)
#  - this was previously done incorrectly and both samples were used
#  - types: living (LD); after cardiac death (DCD); after brain death, with subsequent post-implantation EAD in the recipient (DBD-EAD); 
#        and after brain death without EAD (DBD).

source("processing_utils.R")

gse2 <- getGEOData("GSE23649")
gse23649 <- gse2$originalData$GSE23649  

pTable_ds2 <- gse23649$pheno # - ILLUMINA

# controls, or time point A from not control
pTable_ds2.2 <- pTable_ds2[pTable_ds2$`characteristics_ch1.2`=="time: A" |
                             pTable_ds2$`characteristics_ch1.1`=="trasplant donor origin: Control",] # --> 36 samples
pTable_ds2.2$geo_accession

expData_ds2 <- gse23649$expr[,rownames(pTable_ds2.2)] 

is.matrix(expData_ds2)
boxplot(expData_ds2) # looks ok!

# map to entrez
keys.vec2 <- getIlluminaKeys("GSE23649", "GPL6947_HumanHT-12_V3_0_R1_11283641_A.bgx")

# sex label
sexLabels_table2 <- labelSex(expData_ds2, keys.vec2, ychr.genesEntrez$entrezgene, 3) # 14 female, 22 male - this is different than 26/42 - but I'd argue 13x2=26, 19x2=38+3=41 (not sure about last one?)
sexLabels_ds2 <- sexLabels_table2$sex
names(sexLabels_ds2) <- sexLabels_table2$ID ###!! check

# save
datasetInfo2 <- list("expr"=expData_ds2, "pheno"=sexLabels_ds2, "keys"=keys.vec2, "ID"="GSE23649")
save(datasetInfo2, file="data/processed/ds2.RData")
