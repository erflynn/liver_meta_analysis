# Code for loading and extracting sex labels for each of the datasets.
# E Flynn
# Last Updated: November 2018

source('code/processing_utils.R')

### - load a dataset - ###

createMetaObj <- function(my.ds, selectedCols, sexLabels){
  # my.ds - gse object from meta integrator
  # selectedCols - the columns of the expression matrix to keep
  expMat <- my.ds$expr[rownames(my.ds$expr) %in% names(my.ds$keys),]
  
  expMat2 <- expMat[,selectedCols]
  sexLabels2 <- sexLabels[selectedCols] # make sure same order!
  
  # convert sexLabels to classes: 0 is female
  classes <- sapply(sexLabels2, function(x) ifelse(as.character(x)=="female", 0, 1))
  names(classes) <- colnames(expMat2) 
  
  # create a study object
  # remove rows that do not have any mapping
  my.ds$expr <- as.matrix(expMat2)
  my.ds$keys <- sapply(my.ds$keys[rownames(expMat2)], as.character)
  my.ds$pheno <- data.frame(sexLabels2)
  my.ds$class <- classes

  stopifnot(checkDataObject(my.ds, "Dataset"))
  return(my.ds)
}


##### -------- DS1 -------- #####
gse <- getGEOData("GSE61276")
gse61276 <- gse$originalData$GSE

# check keys
table(is.na(gse61276$keys)) # not too bad - 6k missing, 14k present

# extract the pheno data
pTable <- gse61276$pheno
pTable2 <- pTable[pTable$`agegroup:ch1`=="Adult",]

# filter expression
expData <- gse61276$expr[,rownames(pTable2)]

# sex label
sexLabels <- pTable2$`Sex:ch1`
table(sexLabels)
sexLabels <- sapply(sexLabels, function(x) tolower(as.character(x)))
names(sexLabels) <- pTable2$geo_accession

gse61276$add_comment <- "sex labels from metadata"

# update + save
datasetInfo1 <- createMetaObj(gse61276, pTable2$geo_accession, sexLabels)
save(datasetInfo1, file="data/processed/ds1g.RData")


##### -------- DS2 -------- #####
gse2 <- getGEOData("GSE23649")
gse23649 <- gse2$originalData$GSE23649  

# check keys
table(is.na(gse23649$keys))  # 24k NAs, 24k present...

# extract the pheno data
pTable_ds2 <- gse23649$pheno # - ILLUMINA

#   controls, or time point A from not control
pTable_ds2.2 <- pTable_ds2[pTable_ds2$`characteristics_ch1.2`=="time: A" |
                             pTable_ds2$`characteristics_ch1.1`=="trasplant donor origin: Control",] # --> 36 samples

# filter expression
expData_ds2 <- gse23649$expr[,rownames(pTable_ds2.2)] 

# plot expression
boxplot(expData_ds2) # looks ok!

# sex label
sexLabels_ds2 <- labelSex(expData_ds2, gse23649$keys, ychr.genes$hgnc_symbol, 3, plot=FALSE) # 14 female, 22 male - this is different than 26/42 - but I'd argue 13x2=26, 19x2=38+3=41 (not sure about last one?)
sexLabels_ds2t <- imputeSex(gse23649)
table(sexLabels_ds2==sexLabels_ds2t[names(sexLabels_ds2)]) # all match

# update + save
datasetInfo2 <- createMetaObj(gse23649, rownames(pTable_ds2.2), sexLabels_ds2)
save(datasetInfo2, file="data/processed/ds2g.RData")

##### -------- DS3 -------- #####

gse3 <- getGEOData("GSE38941")
gse38941 <- gse3$originalData$GSE38941
table(is.na(gse38941$keys)) #12k NA, 42k present

pData_ds3 <- gse38941$pheno
pData_ds3.2 <- pData_ds3[pData_ds3$`disease state:ch1` == "normal",]

expData_ds3 <- gse38941$expr[,rownames(pData_ds3.2)]
boxplot(expData_ds3) # looks fine

sexLabels_ds3 <- labelSex(expData_ds3, gse38941$keys, ychr.genes$hgnc_symbol, 3, plot=FALSE)
sexLabels_ds3t <- imputeSex(gse38941)
table(sexLabels_ds3==sexLabels_ds3t[names(sexLabels_ds3)])

# save
datasetInfo3 <- createMetaObj(gse38941, rownames(pData_ds3.2), sexLabels_ds3)
save(datasetInfo3, file="data/processed/ds3g.RData")

##### -------- DS4 -------- #####

gse4 <- getGEOData("GSE12720")
gse12720 <- gse4$originalData$GSE12720
table(is.na(gse12720$keys)) #42k present, 12k absent


pData_ds4 <- gse12720$pheno
pData_ds4.2 <- pData_ds4[pData_ds4$"characteristics_ch1.1" == "Biopsy Type = No Manipulation (Donor)" ,]
expData_ds4 <- gse12720$expr[,rownames(pData_ds4.2)]
boxplot(expData_ds4) # looks ok?


sexLabels_ds4 <- labelSex(expData_ds4, gse12720$keys, ychr.genes$hgnc_symbol, 3, plot=FALSE) # looks good!
table(sexLabels_ds4) # 11 female, 10 male - 2 mislabeled?
# --> sanity check: 
#sanityCheckSexLabels(expData_ds4, ds4.keys, sexLabels_table4$sex) # looks like good separation to me --> gonna stick with it
sexLabels_ds4t <- imputeSex(gse12720)
table(sexLabels_ds4==sexLabels_ds4t[names(sexLabels_ds4)])


# save
datasetInfo4 <- createMetaObj(gse12720, rownames(pData_ds4.2), sexLabels_ds4)
save(datasetInfo4, file="data/processed/ds4g.RData")


##### -------- DS5 -------- #####

gse5 <- getGEOData("GSE14323")
gse14323p1 <- gse5$originalData$GSE14323_GPL571 # ignore other platform
table(is.na(gse14323p1$keys)) # 1.7k NA, 20k map

pData_ds5p1 <- gse14323p1$pheno
pData_ds5.2 <- pData_ds5p1[pData_ds5p1$"Tissue:ch1"=="Normal",]

expData_ds5.2 <- gse14323p1$expr[,rownames(pData_ds5.2)]
boxplot(expData_ds5.2)

sexLabels_ds5 <- labelSex(expData_ds5.2, gse14323p1$keys, ychr.genes$hgnc_symbol, 3, plot=TRUE) # looks good!
sexLabels_ds5t <- imputeSex(gse14323p1)
table(sexLabels_ds5==sexLabels_ds5t[names(sexLabels_ds5)])

gse14323p1$expr <- expData_ds5.2
gse14323p1$pheno <- sexLabels_ds5
datasetInfo5 <- createMetaObj(gse14323p1, rownames(pData_ds5.2), sexLabels_ds5)

save(datasetInfo5, file="data/processed/ds5g.RData")


##### -------- DS6 -------- #####
gse6 <- getGEOData("GSE55668")
gse55668 <- gse6$originalData$GSE55668

table(is.na(gse55668$keys))

pData_ds6 <- gse55668$pheno
annotSexLabels_ds6 <- ifelse(pData_ds6$characteristics_ch1=="Sex: Female", "female", "male")

expData_ds6 <- gse55668$expr
expData_ds6.2 <- expData_ds6+3
boxplot(expData_ds6)
boxplot(expData_ds6.2)

# sex label
sexLabels_ds6 <- labelSex(expData_ds6.2, gse55668$keys, ychr.genes$hgnc_symbol, threshold=3, plot=FALSE)
(annotSexLabels_ds6==sexLabels_ds6) ## ALL MATCH


# save
gse55668$expr <- expData_ds6.2 # include adjusted data
gse55668$add_comment <- "Adjusted expr +3"
datasetInfo6 <- createMetaObj(gse55668, colnames(expData_ds6.2), sexLabels_ds6)

save(datasetInfo6, file="data/processed/ds6g.RData")


##### -------- DS7 -------- #####

gse7 <- getGEOData("GSE61260")
gse61260 <- gse7$OriginalData$GSE61260 # EMPTY

pData_ds7 <- read.delim("data/E-GEOD-61260.sdrf.txt") # downloaded from array express
pData_ds7.2 <- pData_ds7[pData_ds7$Characteristics..diseasestatus.=="Normal Control",]
table(pData_ds7.2$Characteristics..Sex.) # 22 female, 16 male
summary(pData_ds7.2$Characteristics..age.) # good -- all adult

# make sure matches exactly!!
cel.names_ds7 <- sapply(pData_ds7.2$Array.Data.File, function(x) paste("data/GSE61260", paste(x, "gz", sep="."), sep="/"))
rawData_ds7 <- read.celfiles(cel.names_ds7)
rawData_ds7 # there is no feature data
tail(exprs(rawData_ds7)[,1:5]) # need to run RMA, log-transform (which order?)

#pmSeq <- pmSequence(rawData_ds7) # I guess these are present/missing
pmsLog <- log2(pm(rawData_ds7))
boxplot(pmsLog) # medians are a little off
expData_ds7 <- oligo::rma(rawData_ds7) 
boxplot(expData_ds7)
# look at effect of sex

# load relevant table data
keys_vec7 <- loadGPL11532SymbolAnnot()

# sex-label
sexLabels_ds7 <- labelSex(exprs(expData_ds7), keys_vec7, ychr.genes$hgnc_symbol, threshold=3, plot=FALSE) # 21f , 17 m - hmmm... one off
sexLabels_ds7t <- imputeSex(gse61260) # doesn't work
table(sexLabels_ds7==sexLabels_ds7t[names(sexLabels_ds7)])

# look at the messed up one
annotLabels7 <- pData_ds7.2[,c("Array.Data.File", "Characteristics..Sex.")]
# reformat this
subDashDot <- function(x) gsub( "-", ".", as.character(x))
removeGz <- function(x) gsub(".gz", "", as.character(x))
annotLabels7$"Array.Data.File" <- sapply(annotLabels7$Array.Data.File, subDashDot)
names(sexLabels_ds7) <- sapply(names(sexLabels_ds7), removeGz)
colnames(annotLabels7) <- c("ID", "AnnotSex")
annotLabels7$AnnotSex <- sapply(annotLabels7$AnnotSex, as.character)
annotLabels7.2 <- annotLabels7$AnnotSex
names(annotLabels7.2) <- annotLabels7$ID
table(annotLabels7.2[names(sexLabels_ds7)]== sexLabels_ds7)
conflicting7 <- names(annotLabels7.2)[annotLabels7.2[names(sexLabels_ds7)] != sexLabels_ds7]
 # one conflicts --> now let's look at this in detail

# sanityCheckSexLabels(expData_ds7, keys_vec7, sexLabels_ds7$sex) # no XIST --> FAILED. hmmm... try something else

# alternate sanity check
rps4y1.probe <- names(keys_vec7)[as.vector(keys_vec7)=="RPS4Y1"]
kdm5d.probe <- names(keys_vec7)[as.vector(keys_vec7)=="KDM5D"]

res2 <- exprs(expData_ds7)[c(as.character(rps4y1.probe), as.character(kdm5d.probe)),]
res3 <- data.frame(t(res2))
res3$ID <-sapply(rownames(res3), function(x) gsub( "-", ".", as.character(x)))
res4 <- inner_join(res3, sexLabels_ds7)
plot(res4[,c(1,2)], col=ifelse(res4$sex=="female", "red", "blue"), pch=ifelse(res4$ID==conflicting7$oldID, "x", "o"), ylab="KDM5D", xlab="RPS4Y1") ### LOOKS GOOD
# there is no way this is female --> keeping this

expData_ds7.2 <- exprs(expData_ds7)
new_colnames <- sapply(colnames(expData_ds7.2), function(x) strsplit(x, "_", fixed=TRUE)[[1]][1])
names(new_colnames) <- NULL
colnames(expData_ds7.2) <- new_colnames
names(sexLabels_ds7) <- sapply(names(sexLabels_ds7), function(x)strsplit(x, "_", fixed=TRUE)[[1]][1])
#sexLab2 <- sexLabels_ds7[sapply(colnames(exprs(expData_ds7)), subDashDot)]
#keys_vec7.2 <- keys_vec7[rownames(exprs(expData_ds7))]
gse61260$"expr"<- expData_ds7.2
gse61260$keys <- keys_vec7
gse61260$pheno <- data.frame(colnames(expData_ds7.2))
rownames(gse61260$pheno) <-colnames(expData_ds7.2)
gse61260$formattedName <- "GSE61260"
gse61260$load_comment <- "Loaded by hand"
gse61260$exp_comment <- "Data was not in log scale originally"
gse61260$key_comment <- "Typical annotation"
gse61260$add_comment <- "Loaded everything by hand, from raw data"

# swap --> all of them are correct

datasetInfo7 <- createMetaObj(gse61260, colnames(gse61260$expr), sexLabels_ds7)
save(datasetInfo7, file="data/processed/ds7g.RData")


##### -------- DS8 -------- #####

# using data from Schadt
expData8 <- read.delim("data/human_liver_cohort/expression/expression.txt", header=TRUE)
phenoData8 <- read.delim("data/human_liver_cohort/phenotype/phenotype.txt", header=TRUE)
head(phenoData8)
phenoData8$ID <- sapply(phenoData8$Patient_ID, function(x) paste(c("X", x), collapse=""))

# filter by age
phenoData8.1 <- phenoData8[phenoData8$ID %in% colnames(expData8),]  # 466
phenoData8.2 <- phenoData8.1[phenoData8.1$AGE_.YRS.>=18,] # 438

# map keys
# map to ENTREZ IDs
gpl4372 <- getGEO(filename="data/annot/GPL4372.annot.gz")
colnames(Table(gpl4372))
feat.labels <- data.frame(apply(Table(gpl4372)[,c("ID", "Gene symbol")], c(1,2), as.character))
feat.labelssep <- separate_rows(feat.labels[,c(1,2)], "Gene.symbol", sep="///")
keys.vec <- sapply(feat.labelssep$`Gene.symbol`, as.character)
names(keys.vec) <- sapply(feat.labelssep$ID, as.character)
keys.vec8 <- keys.vec
head(keys.vec8)
head(expData8$reporterid)

# hmmm have to double check many-to-many...?? what is going on here - there ARE multiple keys with the same value
expData8.1 <- expData8[,colnames(expData8) %in% phenoData8.2$ID]
expData8.1$reporterid <- expData8$reporterid
expData8.2 <- expData8.1[expData8.1$reporterid %in% names(keys.vec8),]

rownames(expData8.2) <- expData8.2$reporterid
expData8.3 <- expData8.2
expData8.3$reporterid <- NULL

boxplot(expData8.3[,1:5])
min(expData8.3, na.rm=TRUE) # -2 is the minimum --> adjust upward
expData8.4 <- expData8.3+2

keys.vec8.2 <- (keys.vec8[rownames(expData8.3)])

# SEX LABEL
#sexLabels8 <- labelSex(expData8.4, keys.vec8.2, ychr.genes$gene, threshold=3) # didn't work - NAs

simpleSexLabel <- function(expData, keys.vec){
  xist.probe <- names(keys.vec)[as.vector(keys.vec)=="XIST"] # XIST, ENSG00000229807
  print(xist.probe)
  rps4y1.probe <- names(keys.vec)[as.vector(keys.vec)=="RPS4Y1"] # RPS4Y1, ENSG00000129824
  print(rps4y1.probe)
  
  xist.vals <- expData[xist.probe,]
  rps4y1.vals <- expData[rps4y1.probe,]
  ds <- data.frame(t(rbind(xist.vals, rps4y1.vals)))
  
  # remove incomplete cases, cluster
  clus <- kmeans(ds[complete.cases(ds),], 2) 
  print(clus$centers)
  male.cluster <- which.max(as.vector(clus$centers[,ncol(clus$centers)]))
  print(male.cluster)
  sexLabels <- ifelse(clus$cluster==male.cluster, "male", "female") 
  plot(apply(expData[xist.probe,], 2, mean), expData[rps4y1.probe,],  ylab="RPS4Y1", xlab="XIST", col=ifelse(sexLabels=="female", "red", "blue"))  # good separation - one big outlier
  
  return(sexLabels)
}
ss <- simpleSexLabel(expData8.4, keys.vec8.2)
length(ss)
table(ss) # 203 female, 235 male - none really appear ambiguous anymore!
names(ss) <- colnames(expData8.3)
ss.tab <- data.frame(ss, names(ss))
colnames(ss.tab) <- c("sex", "ID")
combinedTab8 <- full_join(phenoData8.2[,c("GENDER", "ID")], ss.tab, by="ID")
combinedTab8$GENDER <- sapply(combinedTab8$GENDER, tolower)
table(combinedTab8$GENDER==combinedTab8$sex) # only 18, 420 agree
combinedTab8[combinedTab8$GENDER!=combinedTab8$sex,] # CHECK THESE 18
mismatchedIDs <- combinedTab8[combinedTab8$GENDER!=combinedTab8$sex,]$ID

gse9588 <-list("expr"= expData8.4, "keys"=keys.vec8.2, "formattedName"="GSE9588")
gse9588$add_comment <- "Adjusted expr +2"
sexLabels_ds8 <-  ss.tab$sex
names(sexLabels_ds8) <- ss.tab$ID
datasetInfo8 <- createMetaObj(gse9588, colnames(expData8.4), sexLabels_ds8)
save(datasetInfo8, file="data/processed/ds8g.RData")


##### -------- DS9 -------- #####
gse9 <- getGEOData("GSE89632")
gse89632 <- gse9$originalData$GSE89632

table(is.na(gse89632$keys)) #25k map, 4k do not


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


sexLabels_ds9 <- labelSex(expData_ds9, gse89632$keys, ychr.genes$hgnc_symbol, threshold=3, plot=FALSE)  # 13 female, 11 male
(sexLabels_ds9==genders) # ALL match - yay!


# save
datasetInfo9 <- createMetaObj(gse89632, rownames(pData_ds9.3), sexLabels_ds9)
save(datasetInfo9, file="data/processed/ds9g.RData")

##### -------- DS10 -------- #####

gse10 <- getGEOData("GSE32504")
gse32504 <- gse10$originalData$GSE32504

table(is.na(gse32504$keys)) #30k missing, 20k present

pData_ds10 <- gse32504$pheno

# other demographic info present - smoking status, alcohol, medication, ethnic background
# 109/149 on medications, 29 smokers, 95 drink
# if look at the supplementary table - some of them have cholestasis (24) and a few NAs -  may want to exclude these?
summary(apply(gse32504$pheno[,c("alcohol intake:ch1", "ethnic background:ch1", "gender:ch1", "presurgical medication:ch1", "smoking:ch1")], 2, as.factor))

sapply(pData_ds10[,"age:ch1"], function(x) as.numeric(as.character(x))) # one is young --> take out
pData_ds10[sapply(pData_ds10[,"age:ch1"], as.numeric) < 18, ] # two --> remove these
pData_ds10.2 <- pData_ds10[sapply(pData_ds10[,"age:ch1"], as.numeric) > 18,]
expData_ds10 <- gse32504$expr[,rownames(pData_ds10.2)]
boxplot(expData_ds10[,1:15])


# sex label
sexLabels_ds10<- labelSex(expData_ds10, gse32504$keys, ychr.genes$hgnc_symbol, threshold=3)  # 77 female, 70 male  - good separation
sexLabels_ds10t <- imputeSex(gse32504)
# RUH ROH - exact opposite
table(sexLabels_ds10==sexLabels_ds10t[names(sexLabels_ds10)])
pData_ds10.2$"gender:ch1"
sanityCheckSexLabels(expData_ds10, gse32504$keys, sexLabels_ds10) ### looks great
# well prediction concurs with what I'd expect - one low XIST, one high xist and opp with rsp4y1
## they really messed up labeling!

# save
datasetInfo10 <- createMetaObj(gse32504, rownames(pData_ds10.2), sexLabels_ds10)
save(datasetInfo10, file="data/processed/ds10g.RData")


##### -------- DS11 -------- #####

gse11 <- getGEOData("GSE33116")
gse33116 <- gse11$originalData$GSE33116 

table(is.na(gse33116$keys)) # 20k present, 1.6k missing

# Breast cancer and liver tissue biopsies used to develop and validate an index for liver contamination in metastatic breast cancer
# no paper listed
# no sex labels, no reference
pData_ds11 <- gse33116$pheno
# filter for liver tissue (rest is breast or liver/breast)
pData_ds11.2 <- pData_ds11[pData_ds11$"tissue:ch1"=="Liver biopsy",] # 31

expData_ds11 <- gse33116$expr[,rownames(pData_ds11.2)]
boxplot(expData_ds11) # looks ok?

# some are negative! 
expData_ds11.2 <- expData_ds11+5

sexLabels_ds11 <- labelSex(expData_ds11.2, gse33116$keys, ychr.genes$hgnc_symbol, threshold=3)  # 18 female, 13 male - looks great separation wise

gse33116$add_comment <- "Adjusted expr +5"
gse33116$expr <- expData_ds11.2 
datasetInfo11 <- createMetaObj(gse33116, rownames(pData_ds11.2), sexLabels_ds11)
save(datasetInfo11, file="data/processed/ds11g.RData")


##### -------- DS12 -------- #####

gse12 <- getGEOData("GSE33814")
gse33814 <- gse12$originalData$GSE33814
# Gene expression profiling unravels cancer-related hepatic molecular signatures in steatohepatitis but not in steatosis
# expected: 13 normal
# no sex labels, no reference
pData_ds12 <- gse33814$pheno
pData_ds12.2 <- pData_ds12[pData_ds12$`diagnosis:ch1`=="normal",]
dim(pData_ds12.2)
expData_ds12 <- gse33814$expr[,rownames(pData_ds12.2)]
boxplot(expData_ds12) # looks fine, all positive

# sex label
sexLabels_ds12 <- labelSex(expData_ds12, gse33814$keys, ychr.genes$hgnc_symbol, threshold=3, plot=FALSE)  # 7f, 6m
sexLabels_ds12t <- imputeSex(gse33814)
table(sexLabels_ds12==sexLabels_ds12t[names(sexLabels_ds12)]) # all match
sanityCheckSexLabels(expData_ds12, gse33814$keys, sexLabels_ds12) # not working :/ 

gse33814$expr <- expData_ds12
gse33814$pheno <- sexLabels_ds12

# save
datasetInfo12 <- createMetaObj(gse33814, rownames(pData_ds12.2), sexLabels_ds12)
save(datasetInfo12, file="data/processed/ds12g.RData")


##### -------- DS13 -------- #####

source('processing_utils.R')

gse13 <- getGEOData("GSE48452")
gse48452 <- gse13$originalData$GSE48452
# Human liver biopsy of different phases from control to NASH
# 5m, 7f
# no reference, sex labels present

pData_ds13 <- gse48452$pheno

# filter to get get normal-weight controls
pData_ds13.2 <- pData_ds13[pData_ds13$source_name_ch1=="Control",] # don't want patients after surgery
pData_ds13.2[,c("age:ch1",
                "bmi:ch1")] # all adult, normal bmi - though some "old" - note this as a caveat
boxplot(gse48452$expr[, rownames(pData_ds13.2)]) # looks fine

sexLabels_ds13 <- labelSex(gse48452$expr[, rownames(pData_ds13.2)], gse48452$keys, ychr.genes$hgnc_symbol, threshold=3, plot=FALSE) # 7f, 5m
(sexLabels_ds13 == pData_ds13.2$"Sex:ch1") # all match -- YAY!
sexLabels_ds13t <- imputeSex(gse48452)
table(sexLabels_ds13==sexLabels_ds13t[names(sexLabels_ds13)])

datasetInfo13 <- createMetaObj(gse48452, rownames(pData_ds13.2), sexLabels_ds13)
save(datasetInfo13, file="data/processed/ds13g.RData")

##### -------- DS14 -------- #####
# (originally ds18)
load("data/E-MEXP-3291.eSet.r")
mexp3291 <- study 
pData_ds18 <- pData(mexp3291)

head(pData_ds18)
pData_ds18.2 <- pData_ds18[pData_ds18$Characteristics.DiseaseState.=="normal",]

table(pData_ds18.2$Characteristics.Sex. ) # 9f, 10m - one more woman than listed

### load raw data
cel.names_ds18 <- sapply(pData_ds18.2$Array.Data.File, function(x) paste("data/E-MEXP-3291/", x,  sep="/"))

# this is also a hugene --> use oligo
rawData_ds18 <- read.celfiles(cel.names_ds18)
rawData_ds18 # there is no feature data
tail(exprs(rawData_ds18)[,1:5]) # need to run RMA, log-transform (which order?)

pmsLog <- log2(pm(rawData_ds18))
boxplot(pmsLog) # medians are a bit move-y --> normalize
expData_ds18 <- rma(rawData_ds18) 
boxplot(expData_ds18) # better

# downloaded additional file .zip this from ArrayExpress: https://www.ebi.ac.uk/arrayexpress/files/A-AFFY-183/?ref=E-MEXP-3291
affy183 <- read.csv("data/annot/HuGene-1_0-st-v1.na29.hg18.transcript.csv", comment.char = "#")

# extract symbols
affy183$symbol <- sapply(affy183$gene_assignment, function(x) {
  gene.info <- strsplit(as.character(x), " /// ", fixed=TRUE)[[1]]; 
  gene.ids <- sapply(gene.info, function(y) trimws(strsplit(y, "//", fixed=TRUE)[[1]][2])); 
  return(paste(unique(gene.ids), collapse=" // "))})

affy183tab <- affy183[,c("probeset_id", "symbol")]
affy183tab.2 <- separate_rows(affy183tab, symbol, sep=" // ")
keys.vec18 <- split(affy183tab.2$symbol, affy183tab.2$probeset_id)

# label sex
sexLabels_ds18 <- labelSex(exprs(expData_ds18), keys.vec18, ychr.genes$hgnc_symbol, threshold=4, plot=TRUE)  # 10 male, 9 female
# note - this switched from before when I wasn't using the raw data... hmm :/ 

# check for concordance
annotSex18 <- pData_ds18.2[,c("Characteristics.Sex.", "Array.Data.File")]
colnames(annotSex18) <- c("annotSex", "ID")
annotSex18$ID <- sapply(annotSex18$ID, as.character)
sexLabels_ds18.2 <- sexLabels_ds18
names(sexLabels_ds18.2) <- sapply(names(sexLabels_ds18.2), function(x) gsub("X", "", as.character(x)))
(annotSex18$annotSex==sexLabels_ds18.2[annotSex18$ID]) # YAY - all match

# expand data because probes map to multiple genes
#head(datasetInfo18$keys[(sapply(datasetInfo18$keys, length)> 1)])

# separate rows effectively
expData_ds18 <- exprs(expData_ds18)
keys.vec18.collapsed <- lapply(keys.vec18, function(x) paste(x, collapse=" // "))
expData_ds18.2 <- data.frame(expData_ds18[names(keys.vec18.collapsed),])
expData_ds18.2$keys <- unlist(sapply(keys.vec18.collapsed,unname))
expData_ds18.2$row.info <- rownames(expData_ds18.2)
expData_ds18.3 <- separate_rows(expData_ds18.2, "keys", sep=" // ")

# rename the rows
keys.vec18.updated <- expData_ds18.3$keys
names(keys.vec18.updated) <- rownames(expData_ds18.3)
expData_ds18.updated <- expData_ds18.3[,1:(ncol(expData_ds18.3)-2)]

# save
emexp3291 <- list("expr"=expData_ds18.updated, "keys"=keys.vec18.updated, "formattedName"="E-MEXP-3291")
datasetInfo14 <- createMetaObj(emexp3291, names(sexLabels_ds18), sexLabels_ds18)
save(datasetInfo14, file="data/processed/ds14g.RData")


##### -------- DS15 -------- #####

gse25935obj <- getGEOData("GSE25935")
# Genetic identification, replication, and functional fine-mapping of expression quantitative trait loci in primary human liver tissue

gse25935 <- gse25935obj$originalData$GSE25935
pData_ds15 <- gse25935$pheno

summary(sapply(pData_ds15$"age:ch1", as.numeric ))
table(pData_ds15$"age:ch1" < 18)
pData_ds15.2 <- pData_ds15[pData_ds15$"age:ch1" >= 18,]
levels(as.factor(sapply(pData_ds15.2$title, function(x) strsplit(as.character(x), "_", fixed=TRUE)[[1]][[2]])))
table(pData_ds15.2rep1$"gender:ch1") # 68f, 115m - 183 total, this makes sense (it says 206)

# - collapse replicates - #
## will have to decide what to do abt reps...
# oh there are two reps per sample
#   either we can: 1) select only one per sample (randomly or just rep1)
#              or: 2) meta-analyze/average over reps... hmm
# samples have anywhere from 1 to 5 reps!!!
# collapse replicates - take the average

pData_ds15.names <- pData_ds15.2[,c("title", "geo_accession", "gender:ch1")]
pData_ds15.names2 <- separate(pData_ds15.names, title, c("sampleID", "replicate"), sep="_")

## AVERAGE BASED ON THIS
sample.to.reps <- split(pData_ds15.names2$geo_accession, pData_ds15.names2$sampleID)

calcSampleAvg <- function(y){
  if(length(y)==1) {
    return(gse25935$expr[,y])
  }
  return(rowMeans(gse25935$expr[,y], na.rm=TRUE))
}
avgSample <- lapply(1:length(sample.to.reps), function(x) { y <- sample.to.reps[[x]]; print(y); calcSampleAvg(y)})
avgSampleDf <- do.call(cbind, avgSample)
colnames(avgSampleDf) <- sapply(names(sample.to.reps), function(x) paste("sample", x, sep=""))
head(avgSampleDf[,1:5])
expData_ds15 <- avgSampleDf

boxplot(expData_ds15) # goes above 16 a little bit
# construct the proper pheno table
pData_ds15.3 <- unique(pData_ds15.names2[,c("sampleID", "gender:ch1")])
(table(pData_ds15.3$sampleID)[table(pData_ds15.3$sampleID) > 1])
pData_ds15.3[pData_ds15.3$sampleID=="325",] # ruh roh, sample is swapping sex within a dataset... eeps!
pData_ds15.3[pData_ds15.3$sampleID=="325",]$`gender:ch1` <- NA # set to NA
pData_ds15.3 <- unique(pData_ds15.3)
rownames(pData_ds15.3) <- sapply(pData_ds15.3$sampleID, function(x) paste("sample", x, sep=""))

# sex label
sexLabels_ds15 <- labelSex(expData_ds15, gse25935$keys, ychr.genes$hgnc_symbol, threshold=3) # doesn't work

gse25935$expr <- expData_ds15
gse25935$pheno <- pData_ds15.3
annotSex_ds15 <- pData_ds15.3$`gender:ch1`
names(annotSex_ds15) <- rownames(pData_ds15.3)
sexLabels_ds15t <- imputeSex(gse25935)
table(sexLabels_ds15t==annotSex_ds15[names(sexLabels_ds15t)]) # all match except swapped one 
# labeled male, appears to be male

## look at these
gse25935$add_comment <- "Averaged over replicates"
gse25935$class <- NULL # remove so this works

# save
datasetInfo15 <- createMetaObj(gse25935, colnames(gse25935$expr), sexLabels_ds15t)
save(datasetInfo15, file="data/processed/ds15g.RData")


##### -------- DS16 -------- #####

gse28893obj <- getGEOData("GSE28893") 
# more of the same as previous dataset - these were used for the eQTL analysis
gse28893 <- gse28893obj$originalData$GSE28893
pData_ds16 <- gse28893$pheno
pData_ds16.2 <- pData_ds16[(pData_ds16$"age:ch1" >= 18),] # removes 4

table(pData_ds16.2$"gender:ch1") # F = 26, M = 30
expData_ds16 <- gse28893$expr[,pData_ds16.2$geo_accession]
boxplot(expData_ds16)
annotLabels_ds16 <- sapply(pData_ds16.2$"gender:ch1", function(x) ifelse(x=="F", "female", "male"))
# for many it says "this sample consists of two replicates" - the GSMs appear to be averaged? --> might want to go back to raw data

# sex label
sexLabels_ds16 <- labelSex(expData_ds16, gse28893$keys, ychr.genes$hgnc_symbol, threshold=3)
(annotLabels_ds16==sexLabels_ds16) # all match!

# save
datasetInfo16 <- createMetaObj(gse28893, pData_ds16.2$geo_accession, sexLabels_ds16)
save(datasetInfo16, file="data/processed/ds16g.RData")


##### -------- DS17 -------- #####
# cannot set this up the same way... 
