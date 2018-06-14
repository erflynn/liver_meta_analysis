# check_compiled_data.R
# E Flynn
# Updated 3/20/2018
#
# Code for checking compiled datasets to see if they have any new samples


# compiled datasets to check:
# E-MTAB-62
# E-MTAB-950
# E-TABM-185
# E-MTAB-3732

# E-MTAB-3732
liver.compiled <- read.csv("compiled_info/liver_mtab3732", header=FALSE) # E-MTAB-62 or 3732 ??? which was this
liver.compiled$study <- sapply(liver.compiled$V2, function(x) strsplit(as.character(x), "_", fixed=TRUE)[[1]][[1]])
table(liver.compiled$study)
levels(liver.compiled$V3)
converted <- sapply(liver.compiled$study, convertArrayExpToGeoID)
fixed.tab[fixed.tab$Accession %in% converted,] # only one is a yes, we are already using it
#  GSE12720

setdiff(converted, fixed.tab$Accession) # "GSE19495" - but cells, only 8 --> DISCARD

liver_62 <- read.delim("compiled_info/liver_mtab62.txt", header=FALSE)
liver_950 <- read.delim("compiled_info/liver_mtab950.txt", header=FALSE)
liver_185 <- read.delim("compiled_info/liver_tabm185.txt", header=FALSE)

liver_62f <- liver_62[liver_62$V2 != "Oliver,,Frank",] # remove misannot (none of his are liver - checked) --> 66
levels(liver_62f$V5)
liver_62f2 <- liver_62f[(liver_62f$V5=="normal solid tissue"),] # only 1
liver_62f2$V3 # E-MEXP-216 - insufficient sample, non-human, non-liver :/ 
## --> discard

head(liver_185)
liver_185f <- liver_185[!liver_185$V6 %in% c("fetal", "fetus" ),] # remove fetal --> 11
View(liver_185f)  ### unclear which samples these are!! 
table(liver_185f$V11) # however, 8male, 2 female, 1 ?? --> insufficient sample

head(liver_950) # 280
levels(liver_950$V12)
levels(liver_950$V13)
levels(liver_950$V19)
liver_950f <- liver_950[liver_950$V19 %in% c(" ", "normal"),]
studies_950 <- unique(liver_950f$V4)
converted_ids2 <- sapply(studies_950, function(x) convertArrayExpToGeoID(x))
converted_ids2

fixed.tab[fixed.tab$Accession %in% converted_ids2,] # only one is a yes, we are already using it
#  GSE12720

setdiff(converted_ids2, fixed.tab$Accession) # none
