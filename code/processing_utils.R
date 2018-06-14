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
# ensembl.m <- useMart("ensembl") #ENSEMBL_MART_ENSEMBL", host="asia.ensembl.org", verbose=T) #, dataset = "hsapiens_gene_ensembl") #, host="archive.ensembl.org")
# listDatasets(ensembl.m)
#ensembl <- useDataset("hsapiens_gene_ensembl", mart=ensembl.m) # strange - did not load at first :/ --> sometimes have to load this twice
# 
# ychr.genesEntrez <- getBM(attributes= "entrezgene",
#                      filters=c("chromosome_name"),
#                      values=list("Y"), mart=ensembl) 

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
  }
}

loadBioconductorEntrez <- function(gpl){
  if (gpl == "GPL570"){
    bioPackage <- 'hgu133plus2.db'
    require(bioPackage)
    x <- hgu133plus2ENTREZID
  }
  if (gpl == "GPL571"){
    bioPackage <- 'hgu133a2.db'
    require(bioPackage)
    x <- hgu133a2ENTREZID
  }
  if (gpl == "GPL96"){
    bioPackage <- 'hgu133a.db'
    require(bioPackage)
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


# maybe pre-process?

# FIND A WAY TO MAP TO ENSEMBL



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
  
  return(sample.results)
}

sanityCheckSexLabels <- function(expData, keys.vec, sex.labels){
  xist.probe <- names(keys.vec)[as.vector(keys.vec)=="7503"] # XIST, ENSG00000229807
  print(xist.probe)
  xist.vals <- expData[xist.probe,]
  if (length(xist.probe)> 1){
    xist.vals <- rowMeans(xist.vals, na.rm=TRUE)
  }
  rps4y1.probe <- names(keys.vec)[as.vector(keys.vec)=="6192"] # RPS4Y1, ENSG00000129824
  print(rps4y1.probe)
  plot(xist.vals, expData[rps4y1.probe,], # TAKE MEANS!
       col=ifelse(sex.labels=="female", "red", "blue"), ylab="RPS4Y1", xlab="XIST")  
}

sexLabelMissing <- function(expData, keys.vec){
  xist.probe <- names(keys.vec)[as.vector(keys.vec)=="7503"] # XIST, ENSG00000229807
  print(xist.probe)
  rps4y1.probe <- names(keys.vec)[as.vector(keys.vec)=="6192"] # RPS4Y1, ENSG00000129824
  print(rps4y1.probe)
  
  # y chromosome probes
  
}


# format into a dataset object

# META-ANALYSIS
# - using Ensembl IDs (??)
# - fixed effects within gene



summarizeGene <- function(key, es){
  es.g <- es[!is.na(es$keys) & es$keys==key,]
  if (nrow(es.g) > 1){
    g.summ <- rma(yi =g, sei=se.g, dat=es.g, measure="SMD", method="FE")
    res <- c("g"=g.summ$b, "se.g"=g.summ$se)
  }
  else {
    res <- c("g"=es.g$g, "se.g"=es.g$se.g)
  }
  return(res)
}


metaAnalyzeGene <- function(key, ds.list, method="DL"){
  
  # check that gene is in all!
  
  gene.sum <- data.frame(do.call(rbind, lapply(ds.list, function(x) summarizeGene(key, x))))
  sum.res <- rma( yi=g, sei=se.g, dat=gene.sum, measure="SMD", method=method)
  return(sum.res)
  #return(list("key"=key, "g"=sum.res$b, "se.g"=sum.res$se, "p"=sum.res$pval,
  #            "ci.lb"=sum.res$ci.lb, "ci.ub"=sum.res$ci.ub, "tau2"=sum.res$tau2, "se.tau2"=sum.res$se.tau2, 
  #            "I2"=sum.res$I2, "H2"=sum.res$H2, "R2"=sum.res$R2)) ### TODO EXTRACT MORE - INCLUDING pval, heterogeneity!!!
  
}