library(openxlsx)
library(tidyverse)
library(splines)
library(parallel)
library(ggplot2)
library(tidytext)
library(ggrepel)
library(geomtextpath)




source('./util_funcs.R')

## Read in microarray data from Michael White. 
count.file <-  "../Input/toxo_cdc/micorarray_counts/Cyclic.Genes.table.xlsx"
expmnt.file <- "../Input/toxo_cdc/micorarray_counts/Samples.xlsx"

count <- read.xlsx(count.file)
expmnt <- read.xlsx(expmnt.file)
GeneID <- count$GeneName
x <- count[,colnames(count) %in% expmnt$Sample[5:28]]
rownames(x) <- GeneID

## Generate time-course data
Samples <-  expmnt[5:28, ]  %>% transmute(Sample = Sample, time_point = as.numeric(gsub('h', '', Time))) %>% distinct()

# ## Readjusting the start time
# ## Microarray starts at S phsase. 6h corresponds to 0 (!!!double check this with figures in the paper)
# Samples$time_point <- (Samples$time_point + 6) %% 12

tc_expr <- x %>% as.data.frame() %>%
  mutate(GeneID = rownames(x)) %>%
  pivot_longer(-GeneID, names_to = "Sample", values_to = "expr")

tc_expr <- right_join(tc_expr, Samples, by = 'Sample')

tc_expr.rep <- tc_expr %>% group_by(GeneID, time_point) %>% summarise(rep = 1:n())
tc_expr$rep <- tc_expr.rep$rep


#saveRDS(tc_expr, '../Input/toxo_cdc/rds/tg_microarray_tc_expr.rds')


num.cores <- detectCores(all.tests = FALSE, logical = TRUE)

micoarray.genes <- unique(tc_expr$GeneID)
microarray.rna.spline.fits <- mclapply(1:length(micoarray.genes), function(i){
  tmp <- tc_expr %>% dplyr::filter(GeneID == micoarray.genes[i]) %>%
    transmute(GeneID = GeneID, x = time_point, y = expr) %>% arrange(x)
  
  y <- tmp$y
  t <- tmp$x
    
  sc.rna.sp <- smooth.spline(t, y, lambda = 0.001)
  
  sc.rna.sp <- predict(sc.rna.sp, seq(0, 12, by = 1/3)) 
  
  ## Filtering the data to set the start at G
  ## According to the paper, G starts at 4.3 (!!! double check)
  sc.rna.sp.filt <- sc.rna.sp %>% data.frame() %>% dplyr::filter(x > 4.3)
  sc.rna.sp.filt$x <- ((sc.rna.sp.filt$x - min(sc.rna.sp.filt$x))/(max(sc.rna.sp.filt$x) - min(sc.rna.sp.filt$x))) * 6
  #plot(sc.rna.sp.filt, type = 'l')
  mu <- data.frame(x = sc.rna.sp.filt$x, y = sc.rna.sp.filt$y) 
  mu <- data.frame(GeneID = rep(tmp$GeneID[1], length(mu[,1])), x = mu[,1], y = mu[,2])
  return(mu)
}, mc.cores = num.cores)

microarray.rna.spline.fits <- bind_rows(microarray.rna.spline.fits)


## Map the IDs to ME49
GT1_ME49 <- read.xlsx('../Input/toxo_genomics/Orthologs/TGGT1_ME49 Orthologs.xlsx')
microarray.rna.spline.fits <- left_join(microarray.rna.spline.fits, GT1_ME49, by = c('GeneID' = 'TGGT1'))

saveRDS(object = microarray.rna.spline.fits ,file = "../Input/toxo_cdc/rds/microarray_splines.rds")


## Calculate alignment and correlation with single cell data
sc.rna.spline.fits <- readRDS('../Input/toxo_cdc/rds/sc_rna_spline_fits_all_genes.rds')

microarray.rna.spline.fits <- microarray.rna.spline.fits %>% 
  transmute(GeneID = gsub('_', '-', TGME49), x = x, y = y)


## Turn the data into wide format (time by gene) and center & scale each gene
sc.rna.wide <- sc.rna.spline.fits %>% 
  pivot_wider(names_from = 'GeneID', values_from = 'y') %>% 
  mutate_at(vars(matches("TGME")), scale) %>%
  as.data.frame()

microarray.rna.wide <- microarray.rna.spline.fits %>% 
  pivot_wider(names_from = 'GeneID', values_from = 'y') %>% 
  mutate_at(vars(matches("TGME")), scale) %>%
  as.data.frame()


## Convertback to long format
sc.rna.mu.scale <- sc.rna.wide %>% 
  pivot_longer(-x, names_to = 'GeneID', values_to = 'expr')

microarray.mu.scale <- microarray.rna.wide %>% 
  pivot_longer(-x, names_to = 'GeneID', values_to = 'expr')


## Remove NA's in sc data. Double check
sc.rm <- unique(sc.rna.mu.scale$GeneID[which(is.nan(sc.rna.mu.scale$expr))])
sc.rna.mu.scale <- sc.rna.mu.scale %>% dplyr::filter(!(GeneID %in% sc.rm))

## Filter to include markers only
common.genes <- intersect(unique(sc.rna.mu.scale$GeneID), unique(microarray.mu.scale$GeneID))

sc.rna.mu.scale <- sc.rna.mu.scale %>% dplyr::filter(GeneID %in% common.genes)
microarray.mu.scale <- microarray.mu.scale %>% dplyr::filter(GeneID %in% common.genes)


cc.dat <- mclapply(1:length(common.genes), function(i){
  microarray.g <- microarray.mu.scale %>% dplyr::filter(GeneID == common.genes[i]) %>% 
    transmute(t = x, y = expr) %>% arrange(t)
  sc.g <- sc.rna.mu.scale %>% dplyr::filter(GeneID == common.genes[i]) %>% 
    transmute(t = x, y = expr) %>% arrange(t)
  #plot(microarray.g$t, microarray.g$y, col = 'red', type = 'l')
  #points(sc.g$t, sc.g$y, col = 'blue', type = 'l')
  tmp <- ccf(c(microarray.g$y), c(sc.g$y), plot = F)
  L <- list(Lag = tmp$lag[which.max(tmp$acf)], cc = max(tmp$acf))
  return(L)
}, mc.cores = num.cores)

cc.dat <- data.frame(GeneID = common.genes, Lags = unlist(lapply(cc.dat, `[[`, 1)), ccs = unlist(lapply(cc.dat, `[[`, 2)))

saveRDS(cc.dat, '../Input/toxo_cdc/rds/sc_rna_microarray_cross_cor_lag.rds')
ggplot(cc.dat, aes(x = ccs)) + 
  geom_histogram(aes(y = ..density..),
                 colour = 1, fill = "white") +
  geom_density(lwd = 1.2,
               linetype = 2,
               colour = 2)



