
# ds18 - E-MEXP-3291 --> NOW DS14
source("processing_utils.R")


# downloaded .r directly from array express - using this for pheno data, will load raw data
load("data/affy/E-MEXP-3291.eSet.r")
mexp3291 <- study 
pData_ds18 <- pData(mexp3291)

head(pData_ds18)
pData_ds18.2 <- pData_ds18[pData_ds18$Characteristics.DiseaseState.=="normal",]

table(pData_ds18.2$Characteristics.Sex. ) # 9f, 10m - one more woman than listed

### load raw data
cel.names_ds18 <- sapply(pData_ds18.2$Array.Data.File, function(x) paste("data/affy/E-MEXP-3291", x,  sep="/"))

# this is also a hugene --> use oligo
rawData_ds18 <- read.celfiles(cel.names_ds18)
rawData_ds18 # there is no feature data
tail(exprs(rawData_ds18)[,1:5]) # need to run RMA, log-transform (which order?)

pmsLog <- log2(pm(rawData_ds18))
boxplot(pmsLog) # looks fine
# Here, normalization is not necessary because the median of the samples is similar, and the data is already in log scale because the expression values are between 0 and 15. (If negative expression values would be observed e.g. the lowest expression value is -1, we recommend to shift all expression values of the dataset above 0 by adding +1 to each gene expression measurement in all samples.)
expData_ds18 <- rma(rawData_ds18) 

#checked: Characteristics.DevelopmentalStage., Characteristics.Age. # --> all adult
#length(unique(pData_ds19.2[,"Characteristics.Individual."])) # all different people! but check for overlap b/c dist system
#expData_ds19 <- exprs(mexp3291)
#colnames(expData_ds19) <- sapply(colnames(expData_ds19), function(x) strsplit(x, ".", fixed=TRUE)[[1]][1])
#expData_ds19.2 <- expData_ds19[,sapply(pData_ds19.2$Source.Name, as.character)]
# NEED TO GET FVAR DATA - but it's empty
#fData(mexp3291)

# downloaded additional file .zip this from ArrayExpress: https://www.ebi.ac.uk/arrayexpress/files/A-AFFY-183/?ref=E-MEXP-3291
affy183 <- read.csv("HuGene-1_0-st-v1.na29.hg18.transcript.csv", comment.char = "#")

# extract entrez IDs
affy183$entrez_ids <- sapply(affy183$gene_assignment, function(x) {
  gene.info <- strsplit(as.character(x), " /// ", fixed=TRUE)[[1]]; 
  gene.ids <- sapply(gene.info, function(y) trimws(strsplit(y, "//", fixed=TRUE)[[1]][5])); 
  return(paste(unique(gene.ids), collapse=" // "))})

affy183tab <- affy183[,c("probeset_id", "entrez_ids")]
affy183tab.2 <- separate_rows(affy183tab, entrez_ids, sep=" // ")
keys.vec18 <- split(affy183tab.2$entrez_ids, affy183tab.2$probeset_id)

# label sex
sexLabels_ds18 <- labelSex(exprs(expData_ds18), keys.vec18, ychr.genesEntrez$entrezgene, threshold=4, plot=TRUE)  # 10 male, 9 female
# note - this switched from before when I wasn't using the raw data... hmm :/ 

# check for concordance
annotSex18 <- pData_ds18.2[,c("Characteristics.Sex.", "Array.Data.File")]
colnames(annotSex18) <- c("annotSex", "ID")
annotSex18$ID <- sapply(annotSex18$ID, as.character)
sexLabels_ds18$ID <- sapply(sexLabels_ds18$ID, function(x) gsub("X", "", as.character(x)))
compareTable18 <- join(annotSex18, sexLabels_ds18[,c("ID","sex")], "ID")
(compareTable18$annotSex==compareTable18$sex) # YAY - all match


# expand data because probes map to multiple genes
head(datasetInfo18$keys[(sapply(datasetInfo18$keys, length)> 1)])

# separate rows effectively
expData_ds18 <- exprs(expData_ds18)
keys.vec18.collapsed <- lapply(keys.vec18, function(x) paste(x, collapse=" // "))
expData_ds18.2 <- data.frame(expData_ds18 )
expData_ds18.2$keys <- unlist(sapply(keys.vec18.collapsed,unname))
expData_ds18.2$row.info <- rownames(expData_ds18.2)
expData_ds18.3 <- separate_rows(expData_ds18.2, "keys", sep=" // ")

# rename the rows
keys.vec18.updated <- expData_ds18.3$keys
names(keys.vec18.updated) <- rownames(expData_ds18.3)
expData_ds18.updated <- expData_ds18.3[,1:(ncol(expData_ds18.3)-2)]
datasetInfo14 <- list("expr"=expData_ds18.updated, "pheno"=sexLabels_ds18$sex, "keys"=keys.vec18.updated, "ID"="E-MEXP-3291")
save(datasetInfo14, file="data/processed/ds14.RData")


