# "README.md"
### author: "E Flynn"
### date: "4/27/2018"


#### Human liver meta-analysis project. 

Performs Meta-analysis of liver expression data.


The directory `search_code` contains the code that was used for searching GEO and processing the search results.
* The results of this are located in `data/search_res`

The data was pre-processed using the file: `load_extract_sex_label.R`. Briefly, each dataset was downloaded, filtered for normal liver samples, expression data normalized and mapped to genes, and sex labels extracted. 
* This results in the files ds1g - ds17g.RData located in `data/processed/`
* Each of these files contains an .RData object with a particular processed dataset.

The bulk of the analysis is then done using: `LiverMA.Rmd`.
