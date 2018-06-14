# preprocess_ds4.R
# 
# ds4 - GSE12720, GPL570
# Transcription profiling of human adult to adult living donor liver grafts compared to deceased donor grafts	
# Expectation: 13m, 8f - from paper, not labeled
# Notes:
#   - 63 samples total
#   - use samples PRE, issue: some deceased - note, make sure these are the DONOR grafts, not the recipient. 
#   - Donor are all adult, disease-free, 13 deceased, 8 living
# can ignore "characteristics_ch1.2" == "HCV - " this info is about the patient receiving, not the donor!!
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2734955/table/T1/ for more demographic info
# no sex labels, no reference info

source('code/processing_utils.R')

gse12720 <- getGEOData("GSE12720")
pData_ds4 <- gse12720$originalData$GSE12720$pheno
pData_ds4.2 <- pData_ds4[pData_ds4$"characteristics_ch1.1" == "Biopsy Type = No Manipulation (Donor)" ,]
cel.names_ds4 <- stripGSM(pData_ds4.2$supplementary_file)
res_ds4 <- loadAffyData(cel.names_ds4, "GSE12720")
expData_ds4 <- res_ds4$exp
qc_ds4 <- res_ds4$qc
plot(qc_ds4) # red for ratios, percentages, and a couple lines :/
# not sure what to do
boxplot(expData_ds4) # looks ok?

# map to ENTREZ, label sex
ds4.keys <- loadEntrezAnnot("GPL570")
sexLabels_table4 <- labelSex(expData_ds4, ds4.keys, ychr.genesEntrez$entrezgene, 3, plot=TRUE) # looks good!

table(sexLabels_table4$sex) # 11 female, 10 male - 2 mislabeled?
# --> sanity check: 
sanityCheckSexLabels(expData_ds4, ds4.keys, sexLabels_table4$sex) # looks like good separation to me --> gonna stick with it
sexLabels_ds4 <- sexLabels_table4$sex
names(sexLabels_ds4) <- sexLabels_table4$ID

# save
datasetInfo4 <- list("expr"=expData_ds4, "pheno"=sexLabels_ds4, "keys"=ds4.keys, "ID"="GSE12720")
save(datasetInfo4, file="data/processed/ds4.RData")
