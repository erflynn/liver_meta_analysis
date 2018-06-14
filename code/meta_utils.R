# meta_utils.R
# E Flynn
# Last Updated: 4/14/2018
#
# Utilities for running meta-analysis.
# Note: this is code adapted from the Khatri lab bitbucket to do what I need.
#       This is not entirely my own.


require('metafor')
require('miceadds')
require('rjson')

json_data <- fromJSON(paste(readLines("external_data/symbols-human.json"), collapse=""))
json.df <- stack(json_data)
head(json.df)
colnames(json.df) <- c("key", "gene")

### ------ code from Khatri lab bitbucket ------ ####
# https://bitbucket.org/khatrilab/meta-analysis-comparison-code/
cleanNA <- function(x) return( x[!is.na(x) & is.finite(x) ] )
getES <- function(v, g){
  stopifnot( identical( length(v), length(g) ) )
  
  x <- cleanNA( v[ which(g==1) ] )
  y <- cleanNA( v[ which(g==0) ] )
  
  n1 <- length(x); n2 <- length(y)
  if( n1 < 2 | n2 < 2 )
    return( c(n1=NA, m1=NA, sd1=NA,
              n2=NA, m2=NA, sd2=NA,
              diff=NA, pooled.sd=NA,
              g=NA, se.g=NA) )
  
  m1   <- mean(x); 
  m2 <- mean(y)
  diff <- m1 - m2
  
  sd1  <- sd(x);  
  sd2 <- sd(y);
  
  sp   <- sqrt( ( (n1-1)*sd1^2 + (n2-1)*sd2^2 )/( n1 + n2 - 2 ) )
  
  cf   <- 1 - 3/( 4*(n1 + n2) - 9 )
  g    <- cf * diff/sp
  se.g <- sqrt( (n1+n2)/(n1*n2) + 0.5*g^2 /(n1+n2-3.94) )
  
  return( c(n1=n1, m1=m1, sd1=sd1,
            n2=n2, m2=m2, sd2=sd2,
            diff=diff, pooled.sd=sp,
            g=g, se.g=se.g) )
}

### ---- adapted code ---- ###

summarizeGene <- function(key, es){
  # use fixed-effects meta-analysis to summarize probes to gene
  #  key = Entrez ID for gene
  #  es = table containing effect size info for a gene
  #         note: this table contains the follow columns: 
  #                 key, g, se.g
  es.g <- es[!is.na(es$keys) & es$keys==key,]
  #print(es.g)
  if (nrow(es.g[!is.na(es.g$g),]) == 0){
    return(NA)
  }
  
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
  #print(key)
  # check that gene is in all!
  gene.sum <- data.frame(do.call(rbind, lapply(1:length(ds.list), 
                                               function(x) { 
                                                 #print(x); 
                                                 summarizeGene(key, ds.list[[x]])})))
  sum.res <- rma( yi=g, sei=se.g, dat=gene.sum, measure="SMD", method=method)
  return(sum.res)
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
  sexLabels <- my.ds$pheno 
  names(sexLabels) <- colnames(my.ds$expr)
  
  # convert sexLabels to classes: 0 is female
  classes <- sapply(sexLabels, function(x) ifelse(as.character(x)=="female", 0, 1))
  names(classes) <- colnames(my.ds$expr) 
    
  # create a study object
  # remove rows that do not have any mapping
  expMat <- my.ds$expr[rownames(my.ds$expr) %in% names(my.ds$keys),]
  study <- list("expr"=as.matrix(expMat), "class"=classes, "keys"=sapply(my.ds$keys[rownames(expMat)], as.character), 
                "formattedName"=my.ds$ID, "pheno"=data.frame(sexLabels))
  stopifnot(checkDataObject(study, "Dataset"))
  return(study)
}

extractStudy <- function(idx){
  # Extracts the study with the given assigned number and reformats and
  # computes the effect sizes for that study

  print(idx)
  load.Rdata(sprintf("data/processed/ds%s.RData", idx), "my.ds")
  # this is based on MetaIntegrator
  # dataset object contains:
  #   expr = expression matrix, rownames are probe IDs, colnames are samples
  #   pheno = list of sex labels, names are column names
  #   keys = named list of key values 
  #   summ = effect summaries (if calculated)
  #   ID = geo/arrayexpress identifier


  # use effect sizes if it already is there
  if("summ" %in% names(my.ds)){
    return(my.ds$summ)
  }
  
  # TODO: check that this matches labels on data  
  sexLabels <- my.ds$pheno 

  # convert sexLabels to classes: 0 is female
  classes <- sapply(sexLabels, function(x) ifelse(as.character(x)=="female", 0, 1))
  
  # create a study object
  study <- list("expr"=my.ds$expr, "class"=classes, "keys"=my.ds$keys)

  # compute effect sizes
  summ <- t( apply( study$expr, 1, getES, g=study$class ) )

  # TODO: make sure the labels on the keys match the data
  res <- sapply(study$keys, as.character)
  names(res) <- sapply(names(res), as.character)
  res2 <- res[rownames(study$expr)]

  # create an object with the summarized effect sizes
  summ2 <- data.frame(summ, keys=unname(res2))
  return(summ2)
}