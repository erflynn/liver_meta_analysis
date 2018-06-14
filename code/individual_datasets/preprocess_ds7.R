
source("code/processing_utils.R")

# ds7 - GSE61260, GPL11532 ** issues with expression data - tried loading it and it was empty --> load from raw data

# Human liver gene expression data from subjects of varying ages
# part of superseries GSE61256 - other subseries are methylation or other tissue
# 40 normal BMI, 134 total
# double check on age, BMI info available
# disease status: : NAFLD (nonalcoholic fatty liver disease),
# NASH (nonalcoholic steatohepatitis), healthy obese, control (for healthy control subjects), primary biliary cirrhosis), and primary
# sclerosing cholangitis.
# Note from paper: Normal control samples were recruited from samples obtained for exclusion of liver malignancy
# during major oncological surgery. None of the normal control individuals underwent preoperative chemotherapy, and liver
# histology demonstrated absence of both cirrhosis and malignancy.
# Patients with evidence of viral hepatitis, hemochromatosis, or alcohol consumption greater than 20 g/d for women and 30 g/d for men were excluded


pData_ds7 <- read.delim("data/affy/E-GEOD-61260.sdrf.txt") # downloaded from array express
pData_ds7.2 <- pData_ds7[pData_ds7$Characteristics..diseasestatus.=="Normal Control",]
table(pData_ds7.2$Characteristics..Sex.) # 22 female, 16 male
summary(pData_ds7.2$Characteristics..age.) # good -- all adult

# make sure matches exactly!!

#expData_ds7 <- loadAffyData(cel.names_ds7, "GSE61260") # this fails
#The affy package can process data from the Gene ST 1.x series of arrays,
#but you should consider using either the oligo or xps packages, which are specifically
#designed for these arrays
# --> use oligo

cel.names_ds7 <- sapply(pData_ds7.2$Array.Data.File, function(x) paste("data/affy/GSE61260", paste(x, "gz", sep="."), sep="/"))
rawData_ds7 <- read.celfiles(cel.names_ds7)
rawData_ds7 # there is no feature data
tail(exprs(rawData_ds7)[,1:5]) # need to run RMA, log-transform (which order?)

#pmSeq <- pmSequence(rawData_ds7) # I guess these are present/missing
pmsLog <- log2(pm(rawData_ds7))
boxplot(pmsLog) # looks fine
# Here, normalization is not necessary because the median of the samples is similar, and the data is already in log scale because the expression values are between 0 and 15. (If negative expression values would be observed e.g. the lowest expression value is -1, we recommend to shift all expression values of the dataset above 0 by adding +1 to each gene expression measurement in all samples.)
expData_ds7 <- oligo::rma(rawData_ds7) 
#save(expData_ds7, file="data/robj/expData_ds7.RData")

# look at effect of sex


# load relevant table data
keys.vec7 <- loadGPL11532EntrezAnnot()


# sex-label
sexLabels_ds7 <- labelSex(exprs(expData_ds7), keys.vec7, ychr.genesEntrez$entrezgene, threshold=3, plot=FALSE) # 21f , 17 m - hmmm... one off
sexLabels_ds7
# look at the messed up one
annotLabels7 <- pData_ds7.2[,c("Array.Data.File", "Characteristics..Sex.")]
# reformat this
subDashDot <- function(x) gsub( "-", ".", as.character(x))
removeGz <- function(x) gsub(".gz", "", as.character(x))
annotLabels7$"Array.Data.File" <- sapply(annotLabels7$Array.Data.File, subDashDot)
sexLabels_ds7.2 <- sexLabels_ds7
sexLabels_ds7.2$oldID <- sexLabels_ds7.2$ID
sexLabels_ds7.2$ID <- sapply(sexLabels_ds7.2$ID, removeGz)
colnames(annotLabels7) <- c("ID", "AnnotSex")
annotLabels7$AnnotSex <- sapply(annotLabels7$AnnotSex, as.character)

compare.tab7 <- join(sexLabels_ds7.2[,c("ID", "oldID", "sex")], annotLabels7, "ID", type="full")
conflicting7 <- compare.tab7[compare.tab7$sex != compare.tab7$AnnotSex,] # one conflicts --> now let's look at this in detail


# sanityCheckSexLabels(expData_ds7, keys.vec7, sexLabels_ds7$sex) # no XIST --> FAILED. hmmm... try something else

# alternate sanity check
rps4y1.probe <- names(keys.vec7)[as.vector(keys.vec7)=="6192"]
kdm5d.probe <- names(keys.vec7)[as.vector(keys.vec7)=="8284"]

res2 <- exprs(expData_ds7)[c(as.character(rps4y1.probe), as.character(kdm5d.probe)),]
res3 <- data.frame(t(res2))
res3$ID <-sapply(rownames(res3), function(x) gsub( "-", ".", as.character(x)))
res4 <- join(res3, sexLabels_ds7)
plot(res4[,c(1,2)], col=ifelse(res4$sex=="female", "red", "blue"), pch=ifelse(res4$ID==conflicting7$oldID, "x", "o"), ylab="KDM5D", xlab="RPS4Y1") ### LOOKS GOOD
# there is no way this is female --> keeping this

colnames(exprs(expData_ds7))
sexLab <- sexLabels_ds7$sex
names(sexLab) <- sexLabels_ds7$ID
sexLab2 <- sexLab[sapply(colnames(exprs(expData_ds7)), subDashDot)]
datasetInfo7 <- list("expr"=exprs(expData_ds7), "pheno"=sexLab2, "keys"=keys.vec7, "ID"="GSE61260")
save(datasetInfo7, file="data/processed/ds7.RData")
