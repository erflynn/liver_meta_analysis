
require('ontologyIndex')
brendaOnt <- get_ontology("../../tissue_data/BrendaTissue.obo",  propagate_relationships=c("is_a", "part_of"))
brenda <- read.table("../../tissue_data/liver_d3_batch.brenda.output")
brenda.rot <- data.frame(t(brenda))
brenda.o <- brenda.rot[order(-brenda.rot[,1]),]

# reformat
rownames(brenda.o) <- sapply(rownames(brenda.o), function(x) paste(strsplit(x, "O")[[1]], collapse="O:"))

# select top hits
top.hits <- rownames(brenda.o)[brenda.o[,1]>0.9] # 42
top.hits.c <- minimal_set(brendaOnt, top.hits) # condenses to 16
  # this may not be what I want... - these are all the most specific nodes

top.hits.label <- sapply(top.hits.c, function(x) get_term_property(ontology=brendaOnt, property="name", term=x))
names(top.hits.label) <- top.hits.c
cbind(brenda.o[names(top.hits.label),1:2], data.frame(top.hits.label))

addLabels <- function(hits, df){
  hits.label <- sapply(hits, function(x) get_term_property(ontology=brendaOnt, property="name", term=x))
  names(hits.label) <-hits
  final.df <- cbind(df[names(hits.label),1:2], data.frame(hits.label))
  return(final.df)
}

# what if instead we look at the shallowest nodes?


# then check which are most child + most parent