library(openxlsx)
library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(gridExtra)

## File downloaded from toxo db
## Searches -> Taxonomy -> ME49 -> download -> tab delimited -> select desired columns
gene.GO <- read.csv('../Input/toxo_genomics/gene_function/ME49_Genes_GO_terms_Summary.csv')

GO.cols <- colnames(gene.GO)[grep ('GO', colnames(gene.GO))]

GO.terms <- lapply(1:(length(GO.cols)/2), function(i){
  varname1 <- GO.cols[(2 * i - 1)]
  varname2 <- GO.cols[(2 * i)]
  
  tmp <- gene.GO %>% dplyr::select(Gene.ID, !!varname1, !!varname2)
  
  ids <- strsplit(c(tmp[,2]), split =";")
  names <- strsplit(c(tmp[,3]), ";")
  ll <- lapply(1:length(ids), function(j){
    length(ids[j][[1]]) == length(names[j][[1]])
  })
  tmp <- tmp[unlist(ll), ]
  tmp <- tmp %>% mutate(!!varname1 := strsplit(tmp[,2], split =";"), !!varname2 := strsplit(tmp[,3], ";")) %>% 
    unnest(-Gene.ID) %>% dplyr::filter(!!as.symbol(varname1) != 'N/A' & !!as.symbol(varname2) != 'N/A')
  colnames(tmp) <- c('GeneID', 'GO.ID', 'GO.name')
  tmp$category <- strsplit(varname1, split = '\\.')[[1]][3]
  tmp$type <- strsplit(varname1, split = '\\.')[[1]][1]
  
  return(tmp)
})


GO.terms <- bind_rows(GO.terms)

write.xlsx(GO.terms, '../Input/toxo_genomics/gene_function/ME49_toxo_db_go_terms.xlsx')
