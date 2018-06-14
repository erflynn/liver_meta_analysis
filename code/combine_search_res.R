# combine_search_res.R
# LiverMetaProject
# E Flynn
# Updated: 4/16/2018
#
# Combine search results and de-duplicate by ID.

setwd("~/Documents/EMILY/Stanford/Coursework/stats211/project/")

# load datasets
array_exp <- read.delim("ArrayExpress-Experiments-180203-011011.txt")
geo <- read.delim("liver_meta.tsv", colClasses = 'character') 

# extract GEO IDs
gses <- sapply(geo$gse, function(x) length(strsplit(as.character(x), "E", fixed=TRUE)[[1]])==2)
geo.ids <- sapply(geo$gse[gses], function(x) strsplit(as.character(x), "E", fixed=TRUE)[[1]][[2]])                  

# convert to GEO ID
array.geo.rows <- sapply(array_exp$Accession, function(x) length(strsplit(as.character(x), "-GEOD-", fixed=TRUE)[[1]])==2)
array.geo.ids <- sapply(array_exp$Accession[array.geo.rows], function(x) strsplit(as.character(x), "-GEOD-", fixed=TRUE)[[1]][[2]])
geo.not.arrayexp <- (setdiff(geo.ids, array.geo.ids)) # 206
length(intersect(geo.ids, array.geo.ids)) # 185
geo.not.arrayexp.gses <- sapply(geo.not.arrayexp, function(x) paste("GSE", x, sep=""))

# put it all together! --> geo file is files only in GEO, ArrayExp contains remainder
geo.filt <- geo[geo$gse %in% geo.not.arrayexp.gses,]
gse.to.gpl <- aggregate(geo.filt[,3], by=list(geo.filt$gse), paste, collapse=",")
colnames(gse.to.gpl) <- c("gse", "gpl")

geo.filt2 <- unique(geo.filt[,1:2])
geo.filt3 <- merge(geo.filt2, gse.to.gpl)

write.table(geo.filt3, file="liver_meta_geo_only.tsv", sep="\t", quote=FALSE, row.names=FALSE)
