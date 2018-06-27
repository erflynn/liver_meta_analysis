


# Qs before generalizing:
# - is MetaIntegrator doing good enough pre-processing?
# - should we do something different?
# - really should meta-analyze probes to genes...

manual_annot <- read.delim("~/Downloads/manual_annotations_ursa (1).csv") # not actually a csv!!
table(manual_annot$BrendaTerm)[order(table(manual_annot$BrendaTerm))]
# --> 215 categories 
# GSE3526 and GSE7307 to label with most similar tissue type

# what about the commas?

length(table(manual_annot$BrendaTerm)[table(manual_annot$BrendaTerm) >=10]) # 140 greater than 10

geoToBatchInput <- function(acc, platform){
  gse <- getGEOData(acc)
  expData <- gse$originalData[[1]]$expr
  pheData <- gse$originalData[[1]]$pheno
  keys <- loadEntrezAnnot(platform)
  
  # use the entrez as keys
  keys2 <- keys[names(keys) %in% rownames(expData)]
  expData2 <- data.frame(cbind(unlist(keys2), expData[names(keys2),]))
  expData2$V1 <-sapply(expData2$V1, as.character)
  expData3 <- separate_rows(expData2, V1, sep=",")
  
  # make it numeric
  expData3[,1:ncol(expData)] <- apply(expData3[,1:ncol(expData)], c(1,2), as.numeric)

  # group by gene
  by_gene <- group_by(expData3, V1)
  mult.cols <- summarise_all(by_gene, .funs=c(mean)) # summarize all
  mult.cols2 <- mult.cols[!is.na(mult.cols$V1),]
  
  # now write out for a batch
  mult.cols.rot <- data.frame(t(mult.cols2))
  colnames(mult.cols.rot) <- sapply(mult.cols.rot[1,], as.character)
  mult.cols.rot <- mult.cols.rot[-1,]
  
  write.table(pheData, file=sprintf("tissue_data/%s_phe.txt", acc), quote=FALSE, sep="\t")
  write.table(mult.cols.rot, file=sprintf("tissue_data/%s_batch.input",acc), 
              row.names=TRUE, col.names=TRUE, sep="\t", quote=FALSE) 
}

geoToBatchInput("GSE3526", "GPL570")

########## - WHAT I DID FOR ONE DATASET - #########

gse38941 <- getGEOData("GSE38941")
pData_ds3 <- gse38941$originalData$GSE38941$pheno
pData_ds3.2 <- pData_ds3[pData_ds3$`disease state:ch1` == "normal",]
normal3 <- rownames(pData_ds3.2)[1]
expData_ds3 <- gse38941$originalData$GSE38941$expr
ds3.keys <- loadEntrezAnnot("GPL570") # use these keys instead - needs to be ENTREZ

#keys3 <- gse38941$originalData$GSE38941$keys

expData_ds3.2 <- data.frame(cbind(unlist(ds3.keys), expData_ds3[names(ds3.keys),]))
expData_ds3.2$V1 <-sapply(expData_ds3.2$V1, as.character)
expData_ds3.3 <- separate_rows(expData_ds3.2, V1, sep=",")

# make it numeric
expData_ds3.3[,1:ncol(expData_ds3)] <- apply(expData_ds3.3[,1:ncol(expData_ds3)], c(1,2), as.numeric)

# turn it into lists

# turn it into batch

# group by gene
by_gene <- group_by(expData_ds3.3, V1)
mult.cols <- summarise_all(by_gene, .funs=c(mean)) # summarize all
mult.cols2 <- mult.cols[!is.na(mult.cols$V1),]

# now write out for a batch
mult.cols.rot <- data.frame(t(mult.cols2))
colnames(mult.cols.rot) <- sapply(mult.cols.rot[1,], as.character)
mult.cols.rot <- mult.cols.rot[-1,]


pData_ds3$`disease state:ch1`
write.table(mult.cols.rot[c(1:5, 20:25),], file="../tissue_data/liver_d3_batch.input", row.names=TRUE, col.names=TRUE, sep="\t", quote=FALSE) # five normal, five dz