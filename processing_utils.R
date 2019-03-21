# processing_utils.R
# E Flynn 
# Last Updated: 4/14/2018
#
# Code for preprocessing expression data for meta-analysis
#
# TODO:
#   - many to many mapping
#   - create common utils - so much repetitive code!
#   - write annotation mapper

# packages
require('MetaIntegrator')
require('massiR')
require('biomaRt')
require('beadarray')
require('GEOquery')
require('rlist')
require('simpleaffy')
require('oligo') # for gene/oligo arrays - some of these are affy
require('pd.hugene.1.1.st.v1')
require('tidyverse')
require('illuminaio')


# identify the y chromosome genes using Entrez
ensembl <- useMart("ensembl", dataset="hsapiens_gene_ensembl") 
# # 
# ychr.genesEntrez <- getBM(attributes= "entrezgene",
#                     filters=c("chromosome_name"),
#                      values=list("Y"), mart=ensembl) 
ychr.genes <- getBM(attributes= "hgnc_symbol",
                         filters=c("chromosome_name"), values=list("Y"), mart=ensembl) 

# Load the y chromosome genes in Entrez --> save these permanently


loadAffyData <- function(cel.names, ID){
  covdesc <- data.frame(row.names = cel.names)
  covdesc$cov <- "n"
  write.table(covdesc, file=sprintf("../data/affy/%s/covdesc", ID), quote=FALSE)
  raw.data <- read.affy(path = sprintf("../data/affy/%s", ID))
  x.rma <- call.exprs(raw.data, "rma") # always returns logged data --> I'm ok!
  myExpData <- exprs(x.rma)
  x.mas5 <- call.exprs(raw.data, "mas5")
  qcs <- qc(raw.data, x.mas5)   
  return(list("exp"=myExpData, "qc"=qcs))
}

stripGSM <- function(list.cel){
  sapply(list.cel, function(x) { y<- strsplit(as.character(x), "/", fixed=TRUE)[[1]]; y[[length(y)]]})
}

# load annotations - TODO: not done!
loadEntrezAnnot <- function(gpl){
  if (gpl %in% c("GPL570", "GPL571", "GPL96")){
    loadBioconductorEntrez(gpl)
  }
  else {
    # FILL IN
    print(sprintf("Error, %s not implemented yet", gpl))
  }
}

loadBioconductorEntrez <- function(gpl){
  if (gpl == "GPL570"){
    require('hgu133plus2.db')
    x <- hgu133plus2ENTREZID
  }
  if (gpl == "GPL571"){
    require('hgu133a2.db')
    x <- hgu133a2ENTREZID
  }
  if (gpl == "GPL96"){
    require('hgu133a.db')
    x <- hgu133aENTREZID
  }

  # Get the probe identifiers that are mapped to an ENTREZ Gene ID
  mapped_probes <- mappedkeys(x)
  # Convert to a list
  xx <- as.list(x[mapped_probes])
  return(xx)
}


loadGPL11532EntrezAnnot <- function(){
  gpl11532 <- getGEO(filename="../data/affy/GPL11532.annot.gz")
  colnames(Table(gpl11532))
  feat.labels <- data.frame(apply(Table(gpl11532)[,c("ID", "Gene ID")], c(1,2), as.character))
  feat.labelssep <- separate_rows(feat.labels[,c(1,2)], "Gene.ID", sep="///")
  keys.vec <- sapply(feat.labelssep$`Gene.ID`, as.character)
  names(keys.vec) <- sapply(feat.labelssep$ID, as.character)
  return(keys.vec)
}

loadGPL11532SymbolAnnot <- function(){
  gpl11532 <- getGEO(filename="data/annot/GPL11532.annot.gz")
  colnames(Table(gpl11532))
  feat.labels <- data.frame(apply(Table(gpl11532)[,c("ID", "Gene symbol")], c(1,2), as.character))
  feat.labelssep <- separate_rows(feat.labels[,c(1,2)], "Gene.symbol", sep="///")
  keys.vec <- sapply(feat.labelssep$`Gene.symbol`, as.character)
  names(keys.vec) <- sapply(feat.labelssep$ID, as.character)
  return(keys.vec)
}


loadGPL4133EntrezAnnot <- function(){
  res <- read.delim("../data/agilent/GPL4133-12599.txt", comment.char="#", header=TRUE, colClasses='character')
  res.df <- res[,c("ID", "GENE")]
  id.to.gene <- split(res.df$GENE, res.df$ID)
  return(id.to.gene)
}

getIlluminaKeys <- function(geo.id, file.name){
  bgx <- readBGX(sprintf("data/illumina/%s/%s", geo.id, file.name))
  map.to.entrez <- (bgx$probes[,c("Probe_Id", "Entrez_Gene_ID")])
  map.to.entrez <- data.frame(apply(map.to.entrez, c(1,2), as.character))
  my.vec <- split(probe.to.entrez$Entrez_Gene_ID, probe.to.entrez$Probe_Id)
  return(my.vec)
}




# SEX LABEL
labelSex <- function(expData, ds.keys, ychr.genes, threshold=3, plot=TRUE){
  
  # check that there are y chromosome probes, if there are none, return an error
  ds.keys2 <- ds.keys[!is.na(ds.keys)]
  ychr.probes <- names(ds.keys2)[as.vector(ds.keys2) %in% ychr.genes]
  if (length(ychr.probes)==0){
    print("error, no y chromosome probes")
    return(NULL)
  }
  
  # label sex
  ychr.df <- data.frame(row.names=unique(ychr.probes))
  ds.y.out <- massi_y(data.frame(expData), ychr.df)
  massi_y_plot(ds.y.out) # -- selecting threshold
  
  massi.select.out <-
    + massi_select(data.frame(expData), ychr.df, threshold=threshold) # default threshold 3 - may want to adjust
  
  head(massi.select.out)[,1:5]
  results <- massi_cluster(massi.select.out)
  sample.results <- data.frame(results[[2]])
  head(sample.results)
  print(table(sample.results$sex))
  if (plot==TRUE){
    massi_cluster_plot(massi.select.out, results) # looks good!
  }
  
  sexLabels <- sample.results$sex
  names(sexLabels) <- sample.results$ID
  
  return(sexLabels)
}

sanityCheckSexLabels <- function(expData, keys.vec, sex.labels){
  keys.vec <- keys.vec[!is.na(keys.vec)]
  xist.probe <- names(keys.vec)[as.vector(keys.vec)=="XIST"] # XIST, ENSG00000229807
  print(xist.probe)
  xist.vals <- expData[xist.probe,]
  if (length(xist.probe)> 1){ xist.vals <- rowMeans(xist.vals, na.rm=TRUE)}
  rps4y1.probe <- names(keys.vec)[as.vector(keys.vec)=="RPS4Y1"] # RPS4Y1, ENSG00000129824
  print(rps4y1.probe)
  rpys4y1.vals <- expData[rpsy41.probe,]
  if (length(rps4y1.probe)> 1){ rps4y1.vals <- rowMeans(rps4y1.vals, na.rm=TRUE)}
  
  plot(xist.vals, rpsy41.vals,
       col=ifelse(sex.labels=="female", "red", "blue"), ylab="RPS4Y1", xlab="XIST")  
}



extractStudyMI <- function(idx){
  print(idx)
  load.Rdata(sprintf("data/processed/ds%s.RData", idx), "my.ds")
  
  # this is based on MetaIntegrator
  # dataset object contains:
  #   expr = expression matrix, rownames are probe IDs, colnames are samples
  #   pheno = list of sex labels, names are column names
  #   keys = named list of key values 
  #   summ = effect summaries (if calculated)
  #   ID = geo/arrayexpress identifier
  
  # TODO: check that this matches labels on data  
  sexLabels <- my.ds$pheno$sexLabels 
  names(sexLabels) <- colnames(my.ds$expr)
  
  # convert sexLabels to classes: 0 is female
  classes <- sapply(sexLabels, function(x) ifelse(as.character(x)=="female", 0, 1))
  names(classes) <- colnames(my.ds$expr) 
  
  # create a study object
  # remove rows that do not have any mapping
  expMat <- my.ds$expr[rownames(my.ds$expr) %in% names(my.ds$keys),]
  my.ds$expr <- as.matrix(expMat)
  my.ds$keys <- sapply(my.ds$keys[rownames(expMat)], as.character)
  my.ds$pheno <- data.frame(sexLabels)
  my.ds$class <- classes
  #study <- list("expr"=as.matrix(expMat), "class"=classes, "keys"=sapply(my.ds$keys[rownames(expMat)], as.character), 
  #              "formattedName"=my.ds$ID, "pheno"=data.frame(sexLabels))
  stopifnot(checkDataObject(my.ds, "Dataset"))
  return(my.ds)
}
