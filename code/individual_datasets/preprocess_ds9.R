source("processing_utils.R")


gse9 <- getGEOData("GSE89632")
gse89632 <- gse9$originalData$GSE89632
# Genome-wide analysis of hepatic gene expression in patients with non-alcoholic fatty liver disease  and in healthy donors in relation to hepatic fatty acid composition and other nutritional factors
# Expectation: 24 healthy, ignore rest
# PMID: https://www.ncbi.nlm.nih.gov/pubmed/25581263
# Description:
#  HCs were approached during their assessments for a live donor liver transplant. These participants underwent transient elastography, computed tomography, and/or magnetic resonance imaging.  Inclusion criteria were as follows: male or female; ≥18 years; for HCs, presence of a normal liver (no steatosis or cirrhosis) on imaging and/or histology;  Exclusion criteria were as follows: alcohol consumption >20 g/day; any other liver disease; use of medications that may cause steatohepatitis, ursodeoxycholic acid or any experimental drug, antioxidants, or PUFA supplements in the 6 months prior to entry; pregnancy or breast-feeding; Liver tissue was collected  as a wedge biopsy during hepatectomy (HC).

# NOTE: RAW DATA AVAILABLE

pData_ds9 <- gse89632$pheno
#"characteristics_ch1.8"== "gender: female"
#title=="liver_HC_HLD-47"
healthy.donors <- sapply(pData_ds9$title, function(x) strsplit(as.character(x), "_", fixed=TRUE)[[1]][[2]]=="HC")
pData_ds9.2 <- pData_ds9[healthy.donors,]

## issue: this pheno table looks misaligned :/ 
# has info on age, gender, BMI that we want to extract
# --> going to dump --> extract
ds9_sample_str <- apply(pData_ds9.2, 1, function(x) paste(x, collapse="\t"))
ages <- sapply(ds9_sample_str, function(x) strsplit(strsplit(x, "age (y): ", fixed=TRUE)[[1]][[2]], "\t", fixed=TRUE)[[1]][[1]])
genders <- sapply(ds9_sample_str, function(x) strsplit(strsplit(x, "gender: ", fixed=TRUE)[[1]][[2]], "\t", fixed=TRUE)[[1]][[1]])
# 13 female, 11 male
bmi <- sapply(ds9_sample_str, function(x) strsplit(strsplit(x, "body mass index (kg/m2): ", fixed=TRUE)[[1]][[2]], "\t", fixed=TRUE)[[1]][[1]])
pData_ds9.3 <- data.frame(cbind("ID"=rownames(ds9_sample_str), "Age"=ages, "Sex"=genders, "BMI"=bmi)) # most normal, a couple >30 - do we want to exclude?
expData_ds9 <- gse89632$expr[,rownames(pData_ds9.3)]
boxplot(expData_ds9) # all looks identical... over-normalized???? :/ 

# map - no bgx available --> use data file downloaded from GEO
gpl14951 <- read.delim("data/illumina/GPL14951-11332.txt", comment.char="#", header=TRUE, colClasses='character')
probe.gene_ds9 <- gpl14951[,c("ID", "Entrez_Gene_ID")]
key.vec_ds9 <-split(probe.gene_ds9$Entrez_Gene_ID, probe.gene_ds9$ID)
head(key.vec_ds9)
sexLabels_ds9 <- labelSex(expData_ds9, key.vec_ds9, ychr.genesEntrez$entrezgene, threshold=3, plot=FALSE)  # 13 female, 11 male
(sexLabels_ds9$sex==genders) # ALL match - yay!

# save
datasetInfo9 <- list("expr"=expData_ds9, "pheno"=sexLabels_ds9$sex, "keys"=key.vec_ds9, "ID"="GSE89632")
save(datasetInfo9, file="data/processed/ds9.RData")

