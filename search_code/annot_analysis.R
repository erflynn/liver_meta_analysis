# annotat_analysis.R
# LiverMetaProject
# E Flynn
# Updated: 6/14/218
#
# Code for examining annotation results
# Notes:
#   - breakdown of samples with multiple criteria?
#   - pull down abstracts + PMIDs too! there is a tool that does this
#   - first pass filtering based on array? title keywords? something else?
#   - count number filtered by title, number filtered by abstract, paper, etc... - when was this done

require('xlsx')

# load annotated datasets
array_exp <- read.xlsx("data/annot/ArrayExpress-Experiments_annot.xlsx", 1)
geo <- read.delim("data/annot/liver_meta_geo_only.tsv", header=TRUE)

# put these together
convertArrayExpToGeoID <- function(acc){
  acc_keys <- strsplit(as.character(acc), "-GEOD-", fixed=TRUE)[[1]]
  if (length(acc_keys) !=2){
    return(as.character(acc))
  }
  id <- acc_keys[[2]]
  return(paste("GSE", id, sep=""))
}
array_exp$Accession2 <- sapply(array_exp$Accession, function(x) convertArrayExpToGeoID(x))

final.col.names <- c("Accession", "Title", "Keep", "Sample.count", "Platform", "Exclusion.criteria", "Notes")
array_exp_sm <- array_exp[,c("Accession2", "Title", "Keep", "Sample.count", "Platform", "Exclusion.criteria", "Notes")]
colnames(array_exp_sm) <- final.col.names
geo_sm <- geo[,c("gse", "title", "keep", "number.of.samples", "gpl", "exclusion.criteria", "notes")]
colnames(geo_sm) <- final.col.names
combined_ds <- rbind(array_exp_sm, geo_sm)



### now upload the updated annotations
keep.reannot <- read.xlsx("data/annot/studies_to_keep.xlsx", 1)
keep.reannot$Accession <- sapply(keep.reannot$Accession, function(x) convertArrayExpToGeoID(x))

# update based on this info
keep.reannot.sm <- keep.reannot[,colnames(combined_ds)]
keep.reannot.sm <- keep.reannot.sm[!is.na(keep.reannot.sm$Accession),]
nrow(keep.reannot.sm)
reannot.ids <- keep.reannot.sm$Accession
table(combined_ds$Accession %in% reannot.ids)
reannotated <- combined_ds[combined_ds$Accession %in% reannot.ids,]
not.reannotated <- combined_ds[!combined_ds$Accession %in% reannot.ids,]
fixed.tab <- rbind( keep.reannot.sm, not.reannotated)
fixed.tab <- fixed.tab[!duplicated(fixed.tab$Accession),]

## make a nicely formatted output

# convert Keep to lowercase
fixed.tab$Keep <- sapply(fixed.tab$Keep, function(x) trimws(tolower(as.character(x))))

## write this out --> this is the properly annotated data
write.table(fixed.tab, file="data/annot/final_annot.tsv", sep="\t", row.names=FALSE, quote=FALSE)


# look at the breakdown of excluded entries
table(fixed.tab$Keep)
excluded.entries <- fixed.tab[fixed.tab$Keep=="no", ]

# format the exclusion criteria
excluded.entries[excluded.entries$Exclusion.criteria=="diseae" & !is.na(excluded.entries$Exclusion.criteria),]$Exclusion.criteria <- "disease"
excl.crit <- sapply(excluded.entries$Exclusion.criteria, function(x) trimws(as.character(x)))

# reorder pairs so they are aggregated
excl.crit2 <- sapply(excl.crit, function(x) {
  my.list <- strsplit(x, ", ", fixed=TRUE)[[1]]; 
  ifelse(length(my.list) >1, paste(my.list[order(my.list)], collapse=","), my.list)})
table(as.factor(excl.crit2))

# fix a couple labels
excl.crit2[excl.crit2=="treated"]  <- "disease"
excl.crit2[excl.crit2=="mixed sex"] <- "pooled"
exclusion <- excl.crit2
excl.tab <-table(exclusion)
excl.tab <- excl.tab[order(-excl.tab)]
par(mar=c(14, 4, 4, 2))
barplot(excl.tab[excl.tab>1], las=2, main="Excluded samples")

# remove multi-entries
multi.excl <- sapply(names(excl.tab), function(x) length(strsplit(x, ",", fixed=TRUE)[[1]]) >1)
multi <- excl.tab[multi.excl] # 147 - meeting multiple criteria
single <- excl.tab[!multi.excl] # 877
par(mar=c(10, 4, 4, 2))

# make a bar plot
barplot(single[single>1], las=2, main="Excluded sample - single criteria")
barplot(excl.tab[excl.tab>2], las=2, main="Excluded samples")

# make a pie chart
other.tot <- sum(excl.tab[excl.tab<10])
simplified.tab <- c(excl.tab[excl.tab>10], "other (<10)"=other.tot)
slices <- simplified.tab
lbls <- names(simplified.tab)
lbls <- sapply(1:length(lbls), function(i) 
  paste(c(lbls[i], "(",slices[i], ")"), collapse="")) 
pie(slices, labels=lbls, col=rainbow(length(lbls)), cex=0.75)


# look at drug stuff
a.keep <- array_exp[grepl("^Yes", array_exp$Keep),]
g.keep <- geo[grepl("^yes", geo$keep),]

a.maybe <- array_exp[grepl("^Maybe", array_exp$Keep),]
g.maybe <- geo[grepl("^maybe", geo$keep),]

a.no <- array_exp[grepl("^No", array_exp$Keep),]
g.no <- geo[grepl("^no", geo$keep),]


levels(a.no$Exclusion.criteria)
levels(g.no$exclusion.criteria)
exclusion <- c(sapply(g.no$exclusion.criteria, as.character), sapply(a.no$Exclusion.criteria, as.character))
trim <- function (x) gsub("^\\s+|\\s+$", "", x)

exclusion <- sapply(exclusion, trim)
table(exclusion)

# count which ones are not included and why
par(mar=c(8, 4, 4, 2))
barplot(table(exclusion)[order(-table(exclusion))][table(exclusion)>2], las=2, main="Excluded samples")

a.drug <- array_exp[grepl("^drug", array_exp$Notes),]
g.drug <- geo[grepl("^drug", geo$notes),]

colnames(g.drug[,1:7])
colnames(a.drug[,1:7])
g.drug2 <- g.drug[,c(1:4, 7, 5:6)]
colnames(g.drug2) <- colnames(a.drug[,1:7])
drug.df <- rbind(a.drug[,1:7], g.drug2)

## write out the drug data 
write.table(drug.df, file="data/annot/drug_exposures.tsv", sep="\t", row.names=FALSE, quote=FALSE)

# -- add a column to the annotation data -- #
#drug_data$Accession <- sapply(drug_data$Accession, convertArrayExpToGeoID)
#annot_data$drug_exposure <- (annot_data$Accession %in% drug_data$Accession)
#write.table(annot_data, file="annot_data_drug.tsv", sep="\t", quote=FALSE, row.names=FALSE)
