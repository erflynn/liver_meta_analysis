

require('MetaIntegrator')
source("code/processing_utils.R")
source("code/meta_utils.R")


# load in all the data

loadStudy <- function(idx){
  load.Rdata(sprintf("data/processed/ds%sg.RData", idx), "my.ds")
  return(my.ds)
}
list_MI_studies <- lapply(1:16, loadStudy)

discoveryStudies <- c(1:7, 9:11) # 434
validationStudies <- c(8, 12:16) # put 8 in validation because it is one of the largest!

# create a meta-object

metaObj <- list()
metaObj$originalData <- list_MI_studies[discoveryStudies]
names(metaObj$originalData) <- sapply(metaObj$originalData, function(x) x$formattedName)
checkDataObject(metaObj, "Meta", "Pre-Analysis")

metaObjv <- list()
metaObjv$originalData <- list_MI_studies[validationStudies]
names(metaObjv$originalData) <- sapply(metaObjv$originalData, function(x) x$formattedName)

checkDataObject(metaObjv, "Meta", "Pre-Analysis")

# run meta-analysis
metaObj <- runMetaAnalysis(metaObj)
metaObjv <- runMetaAnalysis(metaObjv)
save(metaObj, metaObjv, file="data/metaResults_1108.RData")




#Check Distribution of Effect Sizes
# (Separated out + zoomed in for clarity!)

require('reshape2')
ggplot(melt(metaObj$metaAnalysis$datasetEffectSizes[,1:3], 
            varnames = c("Gene", "Study")), aes_string(x = "value", 
                                                       colour = "Study")) + geom_density(size = 1.1) + theme_bw() + 
  scale_color_discrete(name = "Dataset") + xlim(-3,3) # concerned about GSE38941... Wider than normal
ggplot(melt(metaObj$metaAnalysis$datasetEffectSizes[,4:6], 
            varnames = c("Gene", "Study")), aes_string(x = "value", 
                                                       colour = "Study")) + geom_density(size = 1.1) + theme_bw() + 
  scale_color_discrete(name = "Dataset")+ xlim(-3,3)
ggplot(melt(metaObj$metaAnalysis$datasetEffectSizes[,7:10], 
            varnames = c("Gene", "Study")), aes_string(x = "value", 
                                                       colour = "Study")) + geom_density(size = 1.1) + theme_bw() + 
  scale_color_discrete(name = "Dataset")+ xlim(-3,3)

ggplot(melt(metaObjv$metaAnalysis$datasetEffectSizes[,1:3], 
            varnames = c("Gene", "Study")), aes_string(x = "value", 
                                                       colour = "Study")) + geom_density(size = 1.1) + theme_bw() + 
  scale_color_discrete(name = "Dataset")+ xlim(-3,3) # GSE9588 looks narrow, other two look a little wide...
ggplot(melt(metaObjv$metaAnalysis$datasetEffectSizes[,4:6], 
            varnames = c("Gene", "Study")), aes_string(x = "value", 
                                                       colour = "Study")) + geom_density(size = 1.1) + theme_bw() + 
  scale_color_discrete(name = "Dataset")+ xlim(-3,3) # E-MEXP-3291 looks a little wide

### Most of the effect sizes look mostly normal and symmetric, but a few are a little wide, and one is a little narrow


#### --- What are the results? --- ####

metaObj <- filterGenes(metaObj, isLeaveOneOut = TRUE, effectSizeThresh=0.3, FDRThresh = 0.2, numberStudiesThresh = 2)
#metaObj <- filterGenes(metaObj, isLeaveOneOut = FALSE, effectSizeThresh=0.4, FDRThresh = 0.05, numberStudiesThresh = 2)
res <- summarizeFilterResults(metaObj, "FDR0.2_es0.3_nStudies2_looaTRUE_hetero0")
dim(res$neg)


p1 <- violinPlot(metaObj$filterResults$FDR0.2_es0.3_nStudies2_looaTRUE_hetero0, metaObjv$originalData$GSE9588, labelColumn = 'sexLabels2')
p2 <- violinPlot(metaObj$filterResults$FDR0.2_es0.3_nStudies2_looaTRUE_hetero0, metaObjv$originalData$GSE33814, labelColumn = 'sexLabels')
p3 <- violinPlot(metaObj$filterResults$FDR0.2_es0.3_nStudies2_looaTRUE_hetero0, metaObjv$originalData$GSE48452, labelColumn = 'sexLabels')
p4 <- violinPlot(metaObj$filterResults$FDR0.2_es0.3_nStudies2_looaTRUE_hetero0, metaObjv$originalData$E.MEXP.3291, labelColumn = 'sexLabels2')
p5 <- violinPlot(metaObj$filterResults$FDR0.2_es0.3_nStudies2_looaTRUE_hetero0, metaObjv$originalData$GSE25935, labelColumn = 'sexLabels')
p6 <- violinPlot(metaObj$filterResults$FDR0.2_es0.3_nStudies2_looaTRUE_hetero0, metaObjv$originalData$GSE28893, labelColumn = 'sexLabels')
require('gridExtra')
grid.arrange(p1, p2, p3, p4, p5, p6, nrow=2)


# look at chromosome location

require('biomaRt')
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
pos.genes <- metaObj$filterResults$FDR0.2_es0.3_nStudies2_looaTRUE_hetero0$posGeneNames
neg.genes <- metaObj$filterResults$FDR0.2_es0.3_nStudies2_looaTRUE_hetero0$negGeneNames
sig.genes <- c(pos.genes, neg.genes)
#sig.genes2 <- c(metaObj2$filterResults$FDR0.05_es0.4_nStudies2_looaFALSE_hetero0$posGeneNames,
#metaObj2$filterResults$FDR0.05_es0.4_nStudies2_looaFALSE_hetero0$negGeneNames)
mapping.table <- getBM(c('hgnc_symbol','chromosome_name', 'entrezgene'), "hgnc_symbol",sig.genes, mart=ensembl)
mapping.table2 <- mapping.table[mapping.table$chromosome_name %in% c(1:22, "X", "Y"),] # remove the two patch ones
table(mapping.table2$chromosome_name)


# -- plot effect size vs. chromosome! -- #
sigGenes <- rbind(res$pos, res$neg) # 32 up in m, 42 up in f
sigGenes$hgnc_symbol <- rownames(sigGenes)
geneChrES <- full_join(select(sigGenes, effectSize, hgnc_symbol), mapping.table2)
head(geneChrES)
geneChrES$chr_type <- sapply(geneChrES$chromosome_name, function(x) ifelse(x == "X", "X", ifelse(x=="Y", "Y", "A")))

# divide into autosomal vs not
geneChrES2 <- geneChrES[order(geneChrES$effectSize),]
x$name <- factor(x$name, levels = x$name[order(x$val)])

geneChrES2$hgnc_symbol <-factor(geneChrES2$hgnc_symbol, levels=geneChrES2$hgnc_symbol[order(geneChrES2$effectSize)])
ggplot(geneChrES2, aes(y=effectSize, x=hgnc_symbol, fill=chr_type))+ geom_bar(stat="identity") + theme(axis.text.x = element_text(angle = 90, hjust = 1))

y.genes <- mapping.table2[mapping.table2$chromosome_name =="Y",] # 8
x.genes <- mapping.table2[mapping.table2$chromosome_name =="X",] # 13
xy.genes <- mapping.table2[(mapping.table2$chromosome_name %in% c("X","Y")),]
auto.genes <- mapping.table2[!(mapping.table2$chromosome_name %in% c("X","Y")),] # 53
res_xy <- list()
res_xy$posGeneNames <- sapply(geneChrES2[geneChrES2$hgnc_symbol %in% xy.genes$hgnc_symbol & geneChrES2$effectSize > 0,]$hgnc_symbol, as.character) # positive effect size, in x or y chromosome
res_xy$negGeneNames <- sapply(geneChrES2[geneChrES2$hgnc_symbol %in% xy.genes$hgnc_symbol & geneChrES2$effectSize < 0,]$hgnc_symbol, as.character)

res_auto <- list()
res_auto$posGeneNames <- sapply(geneChrES2[geneChrES2$hgnc_symbol %in% auto.genes$hgnc_symbol & geneChrES2$effectSize > 0,]$hgnc_symbol, as.character) 
res_auto$negGeneNames <- sapply(geneChrES2[geneChrES2$hgnc_symbol %in% auto.genes$hgnc_symbol & geneChrES2$effectSize < 0,]$hgnc_symbol, as.character)


### --- look at overlap with validation --- ###
sigGenesValid <- metaObjv$metaAnalysis$pooledResults[sigGenes$hgnc_symbol,]
sigGenesValid2 <- sigGenesValid[sigGenesValid$effectSize > 0.3 | sigGenesValid$effectSize < -0.3,] # 55 of 74
sigGenesValid2[order(sigGenesValid2$effectSize),1:7]

# plot heatmap #
heatmapPlot(metaObj, metaObj$filterResults$FDR0.2_es0.3_nStudies2_looaTRUE_hetero0 )
heatmapPlot(metaObjv, metaObj$filterResults$FDR0.2_es0.3_nStudies2_looaTRUE_hetero0 )


### --- look at overlap with Erika's --- ###
load("../erikaResults/sexMetaObj.RData")
intersect(res_xy$posGeneNames,xySig$negGeneNames) # 8 - Y chromosome + CD99
# "CD99"   "UTY"    "ZFY"    "DDX3Y"  "USP9Y"  "EIF1AY" "KDM5D"  "RPS4Y1"
intersect(res_xy$negGeneNames,xySig$posGeneNames) # 4 - X chromosome
# "XIST"   "RPS4X"  "EIF1AX" "SMC1A" 
intersect(res_auto$posGeneNames,autoSig$negGeneNames) # 1 - LGALS1
#"LGALS1"
intersect(res_auto$negGeneNames,autoSig$posGeneNames) # none



### --- run immunostates --- ###

## try this with an example...
erikaRes <- immunoStatesDecov(sexMetaObj)
erikaISMeta <- immunoStatesMeta(erikaRes)
erikaISMeta <- runMetaAnalysis(erikaISMeta)

liverIS <- immunoStatesDecov(metaObj)
sapply(liverIS$immunoStates, function(x) summary(x$Correlation)) 
# low correlations as expected 
#  range 0.13 to 0.3ish

liverIS2 <- immunoStatesMeta(metaObj)
liverISmeta <- runMetaAnalysis(liverIS2)
# not really sure how to filter

save(liverIS, liverISmeta, file="data/liverISresults.RData")

