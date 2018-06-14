

source("code/processing_utils.R")


# ds8 - GSE9588
# Mapping the Genetic Architecture of Gene Expression in Human Liver
# Note: I'm confused about what might be pooled
#   it looks like samples vs. universal pooled --> this might be tricky to analyze??
# http://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.0060107
#   HLC was assembled from a total of 780 liver samples (1–2 g) that were acquired from Caucasian individuals from three independent liver collections at tissue resource centers at Vanderbilt University, the University of Pittsburgh, and Merck Research Laboratories (Table S1). The Vanderbilt samples (n = 504) included both postmortem tissue and surgical resections from organ donors and were obtained from the Nashville Regional Organ Procurement Agency (Nashville, Tennessee), the National Disease Research Interchange (Philadelphia, Pennsylvania), and the Cooperative Human Tissue Network (University of Pennsylvania, Ohio State University, and University of Alabama at Birmingham). The Pittsburgh samples were normal postmortem human liver and were obtained through the Liver Tissue Procurement and Distribution System (Dr. Stephen Strom, University of Pittsburgh, Pittsburgh, Pennsylvania). The University of Pittsburgh samples (n = 211) were all postmortem, as were the Merck samples (n = 65), which collected by the Drug Metabolism Department and reported previously 
# Demographics table https://doi.org/10.1371/journal.pbio.0060107.st001
#   It's unclear if I can re-identify from the demographics table - I can't (tried this)
# Raw data is NOT available

gse8 <- getGEOData("GSE9588")
gse9588 <- gse8$originalData$GSE9588
head(gse9588$pheno)
gse9588$pheno$title # do not seem to line up with anything :/ 

expData_ds8 <- gse9588$expr
boxplot(expData_ds8)

# labeling is weird too... can we link to demographic information?
#ds8_background <- read.xlsx("data/journal.pbio.0060107.st001.XLS", 1)
#colnames(ds8_background) <- sapply(ds8_background[2,], as.character)
#ds8_background.2 <- ds8_background[3:nrow(ds8_background),] # NOPE, will not link

## - do we have replicates ? - ##
# HH-646-1, HH-646-2
# HH-881-1, HH-881-2 , HH-927-2, HH-952-2  

# sexLabels_ds8 <- labelSex(gse9588$expr, gse9588$keys, ychr.genes$hgnc_symbol, 3 ) # does not work b/c there are multiple NAs

# map to ENTREZ IDs
gpl4372 <- getGEO(filename="data/GPL4372.annot.gz")
colnames(Table(gpl4372))
feat.labels <- data.frame(apply(Table(gpl4372)[,c("ID", "Gene ID")], c(1,2), as.character))
feat.labelssep <- separate_rows(feat.labels[,c(1,2)], "Gene.ID", sep="///")
keys.vec <- sapply(feat.labelssep$`Gene.ID`, as.character)
names(keys.vec) <- sapply(feat.labelssep$ID, as.character)
keys.vec8 <- keys.vec

simpleSexLabel <- function(expData, keys.vec){
  xist.probe <- names(keys.vec)[as.vector(keys.vec)=="7503"] # XIST, ENSG00000229807
  print(xist.probe)
  rps4y1.probe <- names(keys.vec)[as.vector(keys.vec)=="6192"] # RPS4Y1, ENSG00000129824
  print(rps4y1.probe)
  
  plot(apply(expData[xist.probe,], 2, mean), expData[rps4y1.probe,],  ylab="RPS4Y1", xlab="XIST")  # good separation - one big outlier

  xist.vals <- expData[xist.probe,]
  rps4y1.vals <- expData[rps4y1.probe,]
  ds <- data.frame(t(rbind(xist.vals, rps4y1.vals)))
  
  # remove incomplete cases, cluster
  clus <- kmeans(ds[complete.cases(ds),], 2) 
  male.cluster <- which.max(as.vector(clus$centers[,"rps4y1.vals"]))

  sexLabels <- ifelse(clus$cluster==male.cluster, "male", "female") 
  return(sexLabels)
}

sexLabels_ds8 <- simpleSexLabel(expData_ds8, keys.vec8)
table(sexLabels_ds8)

missing.info <- colnames(expData_ds8)[!(colnames(expData_ds8) %in% names(sexLabels_ds8))]


ds.keys2 <- keys.vec8[!is.na(keys.vec8)]
ychr.probes8 <- names(ds.keys2)[as.vector(ds.keys2) %in% ychr.genesEntrez$entrezgene]

missing.ychr <- expData_ds8[ychr.probes8,missing.info]
missing.ychr2 <- missing.ychr[complete.cases(missing.ychr),]
ychr.probes <- row.names(missing.ychr2)
ychr.df <- data.frame(row.names=unique(ychr.probes))
ychr.all <- expData_ds8[ychr.probes, ]
ychr.filt <- data.frame(ychr.all[,!colSums((is.na(ychr.all)))])
min(ychr.filt)
ychr.filt <- data.frame(apply(ychr.filt, c(1,2), function(x) x+2)) # have to adjust up! doesn't work with negative data
ds.y.out <- massi_y(ychr.filt, ychr.df)
massi_y_plot(ds.y.out) # -- selecting threshold

massi.select.out <-
  + massi_select(ychr.filt, ychr.df, threshold=4) # default threshold 3 - may want to adjust

head(massi.select.out)[,1:5]
results <- massi_cluster(massi.select.out)
sample.results <- data.frame(results[[2]])
head(sample.results)
print(table(sample.results$sex)) # looks better
# 164f, 199m

massi_cluster_plot(massi.select.out, results) 

### now look at overlap w other labels
head(sample.results)
sample.results$ID <- sapply(sample.results$ID, as.character)
simpleLabels8 <- data.frame(sexLabels_ds8, names(sexLabels_ds8))
colnames(simpleLabels8) <- c("sex.s", "ID")
compareTab <- join(simpleLabels8, sample.results, "ID", type="full")
labeledDat <- compareTab[!is.na(compareTab$sex ) & (!is.na(compareTab$sex.s)),]
conflicting8 <- labeledDat[labeledDat$sex.s!=labeledDat$sex,] ### one is conflicted

xist.probe <- names(keys.vec8)[as.vector(keys.vec8)=="7503"] # XIST, ENSG00000229807
print(xist.probe)
rps4y1.probe <- names(keys.vec8)[as.vector(keys.vec8)=="6192"] # RPS4Y1, ENSG00000129824
print(rps4y1.probe)
expData_ds8[c(xist.probe, rps4y1.probe), conflicting8$ID] ### AMBIGUOUS - discard!!! probably multiple X chromosomes

compareTab2 <- compareTab[compareTab$ID!=conflicting8$ID,]
selectSexID <- function(x1, x2){ 
  label.list <- c(x1, x2)
  my.label <- unique(label.list[!is.na(label.list)])
  return(my.label)
}
compareTab2$sexLabel <- sapply(1:nrow(compareTab2), function(i) selectSexID(compareTab2$sex.s[[i]], compareTab2$sex[[i]]))
table(compareTab2$sexLabel)
expData_ds8.2 <- data.frame(apply(expData_ds8, c(1,2), function(x) x+2)) # this is artificial - all versus a reference...
# I am not sure what makes sense here
head(expData_ds8.2[,1:5])
boxplot(expData_ds8.2)
# remove the conflicting ID from BOTH sex labels and the expr matrix

expData_ds8.3 <- expData_ds8.2[,colnames(expData_ds8.2)!=as.character(conflicting8$ID)] 
compareTab2.1 <- compareTab2[compareTab2$ID != as.character(conflicting8$ID),]
sexLabels_ds8 <- compareTab2$sexLabel
table(sexLabels_ds8) # 186 female, 225 male

# are there other ambiguous ones we should remove???
datasetInfo8 <- list("expr"=expData_ds8.2, "pheno"=sexLabels_ds8, "keys"=keys.vec8, "ID"="GSE9588")
save(datasetInfo8, file="data/processed/ds8.RData")




# apply(t(expData_ds8[ychr.probes8,missing.info]), 1, function(x) which.max(x))
# 
# 
# ds8imputed <- data.frame(sexLabels_ds8)
# table(ds8imputed$sexLabels_ds8)
# inconclusive <- rownames(ds8[ds8[,1] > 0 & ds8$rps4y1.vals8 > 0,])[1] ## discard this one
# #rownames(ds8[ds8[,2] > 0 & ds8$rps4y1.vals8 > 0,])[1] ## same result with either probe
# 
# table(ds8imputed[rownames(ds8imputed)!=inconclusive,]) #185f    226m = 411 total
# 
# 
# # -- what about the 15 others? -- #
# missing.labels <- rownames(ds8[!complete.cases(ds8),])
# 
# # pull in some extra y chromosome probes
# 
# sex.chr8.df <- data.frame(gse9588$expr[c(ychr8.probes, xist.probes),])
# sex.chr8.df$"GSM242269" <- NULL
# sex.chr8.df[,missing.labels]
# 
# # what are the possible issues:
# # - y chromosome probes don't do a great job discriminating in this case
# # - messy data
# # - XIST weighted too heavily?
# # - weighting of y chromosome
# 
# 
# y.chr8.df <- data.frame(gse9588$expr[c(ychr8.probes),])
# 
# # select probes that DO discriminate
# female.samples <- rownames(ds8imputed)[ds8imputed$sexLabels_ds8=="female"]
# male.samples <- rownames(ds8imputed)[ds8imputed$sexLabels_ds8=="male"]
# max.diff <- which(rowMeans(y.chr8.df[,female.samples], na.rm=TRUE) 
#                   - rowMeans(y.chr8.df[,male.samples], na.rm=TRUE) < -0.1)
# 
# head(y.chr8.df[max.diff,male.samples])
# dim(y.chr8.df[max.diff,female.samples])
# 
# head(y.chr8.df[,1:5])
# 
# simpleLabel <- function(missing.sample){
#   ychr.mean <- mean(ychr8.df[max.diff,missing.sample], na.rm=TRUE)
#   xist.mean <- mean(ychr8.df[xist.probes,missing.sample], na.rm=TRUE)
#   print(ychr.mean) 
#   print(xist.mean)
#   ### todo - think about these cutoffs - should really be clustering BUT missing data
#   sexlabel <- ifelse(xist.mean < -0.1, ifelse(ychr.mean > 0, "male", "NA"), ifelse(ychr.mean < -0.1, "female", "NA"))
#   print(sexlabel)
#   return(list("x"=xist.mean, "y"=ychr.mean, "label"=sexlabel))
# }
# 
# labels.for.missing <- lapply(missing.labels, simpleLabel)
# plotLabelsSimple <- function(new.label.dat){
#   missing.df <- data.frame(do.call(rbind, new.label.dat))
#   
#   missing.df <- data.frame(apply(missing.df, c(1,2), function(x) unlist(x)))
#   missing.df$x <- sapply(missing.df$x, function(x) as.numeric(as.character(x)))
#   missing.df$y <- sapply(missing.df$y, function(x) as.numeric(as.character(x)))
#   print(table(missing.df$label)  )
#   plot(missing.df$y ~ missing.df$x, col=ifelse(missing.df$label=="female", "red", "blue"))  
#   return(missing.df)
# }
# 
#  # missing one - appears to be female, but low XIST
# missing.df <- plotLabelsSimple(labels.for.missing)
# rownames(missing.df) <- missing.labels
# 
# # plot these the same way
# simpleLabels <- lapply(ds8imputed.df$X1, simpleLabel)
# simpleLabelDf <- plotLabelsSimple(c(simpleLabels))
# rownames(simpleLabelDf) <- ds8imputed.df$X1
# table(sapply(simpleLabelDf$label, as.character)==sapply(ds8imputed.df$X2, as.character)) # 181 do not match
# ### THERE IS A PROBLEM HERE!!
# 
# fullDf <- plotLabelsSimple(c(simpleLabels, labels.for.missing)) # solid separation 
# 
# 


### look at the human liver cohort data

expData8 <- read.delim("data/human_liver_cohort/expression/expression.txt", header=TRUE)
phenoData8 <- read.delim("data/human_liver_cohort/phenotype/phenotype.txt", header=TRUE)
head(phenoData8)
phenoData8$ID <- sapply(phenoData8$Patient_ID, function(x) paste(c("X", x), collapse=""))
head(phenoData8$AGE_.YRS.)
males <- phenoData8$ID[phenoData8$GENDER=="Male"]
females <- phenoData8$ID[phenoData8$GENDER=="Female"]
fwithExp <- females[females %in% colnames(expData8)]
fExpr <- expData8[,fwithExp]
meanFemales <- rowMeans(fExpr, na.rm=TRUE)
mwithExp <- males[males %in% colnames(expData8)]
mExpr <- expData8[,mwithExp]
meanMales <- rowMeans(mExpr, na.rm=TRUE)

gender.labels <- sapply(colnames(expData8), function(x) phenoData8$GENDER[phenoData8$ID==x])
colnames(expData8.2) <- expData8.2["reporterid",]
expData8.3 <- data.frame(expData8.2[-c(1:4),])
expData8.3$gender <-unname(gender.labels[5:length(gender.labels)])
head(expData8[,1:5])

xist <- which(expData8$genesymbol=="XIST")
my.df <- data.frame(apply(expData8.3[,1:(ncol(expData8.3)-1)], c(1,2), as.numeric))
head(my.df[,1:5])
# remove columns with NAs
na_count <-sapply(my.df, function(y) sum(length(which(is.na(y)))))
my.df.trim <- my.df[,na_count == 0]
pcs <- prcomp(my.df.trim, retx=TRUE)
gender.labels2 <- unlist(unname(gender.labels)[5:length(gender.labels)])
#colors <- c("red", "blue")[as.factor(gender.labels2)]
#biplot(pcs, col=colors)

#require('pvca')
#pvcaObj <- pvcaBatchAssess(as.matrix(my.df.trim), gender.labels2)

library(ggbiplot)
g <- ggbiplot(pcs, obs.scale = 1, var.scale = 1, 
              groups = gender.labels2, ellipse = TRUE, 
              circle = TRUE)
g <- g + scale_color_discrete(name = '')
g <- g + theme(legend.direction = 'horizontal', 
               legend.position = 'top')
print(g)

require('reshape2')
my.df.trim$ID <- rownames(my.df.trim)
my.df.trim$gender <- gender.labels2


### DO THIS ON PC'd data :/ -- too slow otherwise

plot(pcs, type="l")
cum.sum.var <- (cumsum(pcs$sdev^2 / sum(pcs$sdev^2)))
numpcs <- which(cum.sum.var > 0.95)[[1]] - 1
numpcs
# melt the PC matrix

expData8.filt <- expData8[expData8$genesymbol!="Pro25G",]

pc.mat <- data.frame(pcs$x[,1:numpcs])
pc.mat$gender <- gender.labels2
pc.mat$ID <- rownames(pc.mat)
pc.mat$Age <- phenoData8[phenoData8$ID %in% colnames(expData8),]$AGE_.YRS.
pc.mat$ethnicity <- phenoData8[phenoData8$ID %in% colnames(expData8),]$self_reported_ethnicity
molten <- melt(pc.mat, id.vars=c("gender", "ID", "Age", "ethnicity"))
molten$gender <- as.factor(molten$gender)
fit <- lm(value ~ variable + gender, data=molten)
anova(fit) # SIGNIFICANT
molten2 <- molten
molten2$variable <- as.factor(molten2$variable)

require('ggplot2')
ggplot(molten, aes(variable, value))+geom_boxplot(aes(colour=gender))
boxplot(value ~ variabl, data=molten2, col=c("red", "blue"))
cor(molten[molten$variable=="PC1" & molten$gender=="Male",]$value, molten[molten$variable=="PC1" & molten$gender=="Female",]$value)

pDat <- phenoData8[,c("ID", "GENDER")]

load("data/processed/ds13.RData")

pcs4 <- prcomp(t(datasetInfo13$expr[complete.cases(datasetInfo13$expr),]))
cum.sum.var2 <- (cumsum(pcs4$sdev^2 / sum(pcs4$sdev^2)))
numpcs2 <- which(cum.sum.var2 > 0.8)[[1]] - 1
numpcs2
pc.mat2 <- data.frame(pcs4$x[,1:numpcs2])
pc.mat2$gender <- datasetInfo13$pheno
pc.mat2$ID <- rownames(pc.mat2)
head(pc.mat2)
molten_v2 <- melt(pc.mat2, varnames=colnames(pc.mat2)[1:numpcs2])
fit2 <- lm(value ~ variable + gender, data=molten_v2)
anova(fit2) # SIGNIFICANT
ggplot(molten_v2, aes(variable, value))+geom_boxplot(aes(colour=gender))


load("data/processed/ds8.RData")
ds8.geo <- datasetInfo8$expr[,names(datasetInfo8$pheno)]
pcs8 <- prcomp(t(ds8.geo[complete.cases(ds8.geo),]))
cum.sum.var8 <- (cumsum(pcs8$sdev^2 / sum(pcs8$sdev^2)))
numpcs8 <- which(cum.sum.var8 > 0.8)[[1]] - 1
numpcs8
pc.mat8 <- data.frame(pcs8$x[,1:numpcs8])
pc.mat8$gender <- datasetInfo8$pheno
pc.mat8$ID <- rownames(pc.mat8)
head(molten_v8)
molten_v8 <- melt(pc.mat8, varnames=colnames(pc.mat8)[1:numpcs8])
fit8 <- lm(value ~ variable + gender, data=molten_v8)
anova(fit8) # SIGNIFICANT
ggplot(molten_v8, aes(variable, value))+geom_boxplot(aes(colour=gender))

### IGNORE ABOVE


# process data from Schadt et al. 
expData8 <- read.delim("data/human_liver_cohort/expression/expression.txt", header=TRUE)
phenoData8 <- read.delim("data/human_liver_cohort/phenotype/phenotype.txt", header=TRUE)
head(phenoData8)
phenoData8$ID <- sapply(phenoData8$Patient_ID, function(x) paste(c("X", x), collapse=""))

# filter by age
phenoData8.1 <- phenoData8[phenoData8$ID %in% colnames(expData8),]  # 466
phenoData8.2 <- phenoData8.1[phenoData8.1$AGE_.YRS.>=18,] # 438

# map keys
# map to ENTREZ IDs
gpl4372 <- getGEO(filename="data/GPL4372.annot.gz")
colnames(Table(gpl4372))
feat.labels <- data.frame(apply(Table(gpl4372)[,c("ID", "Gene ID")], c(1,2), as.character))
feat.labelssep <- separate_rows(feat.labels[,c(1,2)], "Gene.ID", sep="///")
keys.vec <- sapply(feat.labelssep$`Gene.ID`, as.character)
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

keys.vec8.2 <- (keys.vec8[rownames(expData8.3)])

# SEX LABEL
sexLabels8 <- labelSex(expData8.3, keys.vec8.2, ychr.genesEntrez$entrezgene, threshold=3) # didn't work - NAs

expData <- expData8.3[,mismatchedIDs]
keys.vec <- keys.vec8.2
simpleSexLabel <- function(expData, keys.vec){
  xist.probe <- names(keys.vec)[as.vector(keys.vec)=="7503"] # XIST, ENSG00000229807
  print(xist.probe)
  rps4y1.probe <- names(keys.vec)[as.vector(keys.vec)=="6192"] # RPS4Y1, ENSG00000129824
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
ss <- simpleSexLabel(expData8.3, keys.vec8.2)
length(ss)
table(ss) # 203 female, 235 male - none really appear ambiguous anymore!
names(ss) <- colnames(expData8.3)
ss.tab <- data.frame(ss, names(ss))
colnames(ss.tab) <- c("sex", "ID")
combinedTab8 <- join(phenoData8.2[,c("GENDER", "ID")], ss.tab, by="ID", type="full")
combinedTab8$GENDER <- sapply(combinedTab8$GENDER, tolower)
table(combinedTab8$GENDER==combinedTab8$sex) # only 18, 420 agree
combinedTab8[combinedTab8$GENDER!=combinedTab8$sex,] # CHECK THESE 18
mismatchedIDs <- combinedTab8[combinedTab8$GENDER!=combinedTab8$sex,]$ID

mismatchedIDs

# this convincing to me based on the data
ds2 <- data.frame(ds, "ID"=rownames(ds))
join(ds2, combinedTab8, by="ID")

ychr8 <- keys.vec8.2[keys.vec8.2 %in% ychr.genesEntrez$entrezgene]
ychrDf <- expData8.3[names(ychr8),]
ychrDf.2 <- ychrDf[complete.cases(ychrDf),]
sexLabels8.2 <- labelSex(apply(expData8.3, c(1,2), function(x) x+2), keys.vec8.2[names(keys.vec8.2) %in% rownames(ychrDf.2)], ychr.genesEntrez$entrezgene, 3) # visualization looks terrible
# only works if you add to the expression vlaues so they're >=0
colnames(sexLabels8.2)[colnames(sexLabels8.2)=="sex"] <- "sex2"
combinedTab8.2 <- join( combinedTab8, sexLabels8.2, by="ID", type="full")
head(combinedTab8.2)
combinedTab8.2[combinedTab8.2$ID %in% mismatchedIDs,] ## this agrees with us --> we keep!!

datasetInfo8 <- list("expr"= expData8.3, "pheno"=ss.tab$sex, "keys"=keys.vec8.2, "ID"="GSE9588")
save(datasetInfo8, file="data/processed/ds8.RData")
