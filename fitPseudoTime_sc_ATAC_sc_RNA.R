library(Seurat)
library(ggplot2)
library(ggfortify)
library(jcolors)
require(gridExtra)
library(grid)
library(slingshot)
library(gam)
library(princurve)
library(parallel)
library(dtwclust)
library(doParallel)
library(tidyverse)
library(tidytext)
library(fda)
library(sme)
library(MyEllipsefit)
library(openxlsx)


#library(sctransform)

source('./util_funcs.R')

num.cores <- detectCores(all.tests = FALSE, logical = TRUE)

## Fit a pseudo-time curve and align using sync data
S.O.integrated <- readRDS('../Input/toxo_cdc/rds/S.O.intra_atac_integrated.rds')
S.O.integrated@meta.data$Sample <- rownames(S.O.integrated@meta.data)
Idents(S.O.integrated) <- 'orig.ident'

atac_sub <- subset(S.O.integrated, ident = 'scATAC')
rna_sub <- subset(S.O.integrated, ident = 'scRNA')


## Pseudo-time analysis with SLingshot

fitTime <- function(S.O, method = 'rna', reverse.t = F){
  pc.tg <- getPCA(S.O)
  sds.data <- getPrinCurve(pc.tg)
  ind <- match(sds.data$Sample, pc.tg$Sample)
  sds.data$PC_1 <- pc.tg$PC_1[ind]
  sds.data$PC_2 <- pc.tg$PC_2[ind]
  sds.data$phase <- S.O@meta.data$phase[match(sds.data$Sample, rownames(S.O@meta.data))]
  sds.data$phase <- factor(sds.data$phase, levels = c('G1.a', 'G1.b', 'S', 'M', 'C'))
  
  Idents(S.O) <- 'phase'
  DimPlot(S.O, reduction = 'pca', label = T)


  pt <- sds.data$pt
  
  sds.data$pt <- 6 * ((as.numeric(pt) - min(as.numeric(pt)))/(max(as.numeric(pt)) - min(as.numeric(pt))))
  
  plot(sds.data$phase, sds.data$pt)
  
  if(reverse.t){
    sds.data$pt <- 6 - sds.data$pt
  }
  
  # Shift the time to start at G1.a
  tmp <- sds.data %>% dplyr::filter(phase == 'G1.a') %>% arrange(pt)
  lag.time <- tmp$pt[which.max((tmp$pt[2:length(tmp$pt)] - tmp$pt[1:(length(tmp$pt) - 1)])) + 1]
  #lag.time <- quantile(tmp$pt, p = 0.205) ## excluse the ones overlapping with G1.b
  
  sds.data$pt <- (sds.data$pt - lag.time + 6) %% 6
  

  plot(sds.data$phase, sds.data$pt)
  
 
  ind.G1.a <- which(sds.data$phase == 'G1.a')
  ind.G1.b <- which(sds.data$phase == 'G1.b')
  ind.S <- which(sds.data$phase == 'S')
  ind.M <- which(sds.data$phase == 'M')
  ind.C <- which(sds.data$phase == 'C')
  
  
  L <- wiskerPlot(S.O)

  par(mar = c(5, 5, 4, 4) + 0.1)
  plot(x = -35:10, y = -25:20, type = 'n', xlab = 'PC1', ylab = 'PC2',  lwd = 2, cex.lab = 1.5, cex.main = 2, cex.axis = 1.5)
  whiskers(as.matrix(L$pc[,c(1,2)]), L$fit$s, col = "gray")
  points(sds.data$PC_1, sds.data$PC_2, cex = 0.6, col = sds.data$phase, pch = 20)
  points(sds.data$sc1[sds.data$cell.ord],sds.data$sc2[sds.data$cell.ord], cex = 0.2, col = 'red')

  
  # Scale the pt based on know biology: Radke et. al 2000
  G <- c(0, 3) # 3h
  S <- c(3, 4.7) # 1.7h
  M <- c(4.7, 5) # ~20 min
  C <- c(5, 6) # 1h
  
  t1 <- (quantile(sds.data$pt[ind.G1.a], prob=0.75) + 
    quantile(sds.data$pt[ind.G1.b], prob=0.25))/2 ## Start of G1.b
  t2 <- (quantile(sds.data$pt[ind.G1.b], prob=0.75) + 
    quantile(sds.data$pt[ind.S], prob=0.25)) / 2 ## End of G1.b, Start of S
  t3 <- (quantile(sds.data$pt[ind.S], prob=0.75) + 
    quantile(sds.data$pt[ind.M], prob=0.25)) / 2## End of S, Start of M
  t4 <- (quantile(sds.data$pt[ind.M], prob=0.75) + 
    quantile(sds.data$pt[ind.C], prob=0.25)) / 2## End of M, Start of C

  t0 <- 0
  t5 <- 6
  
  # sds.data$pt.shift <- (sds.data$pt + t.shift) %% 6
  # sds.data$pt.shift <- 6 * ((sds.data$pt.shift - min(sds.data$pt.shift)) / 
  #                             (max(sds.data$pt.shift) - min(sds.data$pt.shift)))
  # plot(sds.data$phase, sds.data$pt.shift)
  
  slp.g <- (G[2] - G[1]) / (t2 - t0)
  inc.g  <- c(t0, G[1])
  
  slp.s <- (S[2] - S[1]) / (t3 - t2)
  inc.s  <- c(t2, S[1])
  
  slp.m <- (M[2] - M[1]) / (t4 - t3)
  inc.m  <- c(t3, M[1])
  
  slp.c <- (C[2] - C[1]) / (t5 - t4)
  inc.c  <- c(t4, C[1])
  
  sds.data <- sds.data %>%  
    mutate(pt.shifted.scaled = case_when(phase %in% c('G1.a', 'G1.b') ~ inc.g[2] + slp.g * (pt - inc.g[1]),
                                         phase == 'S' ~ inc.s[2] + slp.s * (pt - inc.s[1]),
                                         phase == 'M' ~ inc.m[2] + slp.m * (pt - inc.m[1]),
                                         phase == 'C' ~ inc.c[2] + slp.c * (pt - inc.c[1])))
  
  plot(sds.data$phase, sds.data$pt.shifted.scaled)
  
  ## Exclude outlier samples
  #q.ex <- quantile(sds.data$pt.shifted.scaled, p = 0.998)
  q.ex <- 7
  sds.data <- sds.data %>% dplyr::filter(pt.shifted.scaled <= q.ex)
  ## Rescale to [0, 6]
  sds.data$pt.shifted.scaled <- 6 * ((sds.data$pt.shifted.scaled - min(sds.data$pt.shifted.scaled))/
                                       (max(sds.data$pt.shifted.scaled) - min(sds.data$pt.shifted.scaled)))
  plot(sds.data$phase, sds.data$pt.shifted.scaled)
  # plot(sort(sds.data$pt.shifted.scaled))
  # 

  S.O.filt <- S.O
  S.O.filt@meta.data$Sample <- rownames(S.O.filt@meta.data)
  Idents(S.O.filt) <- 'Sample'
  S.O.filt <- subset(S.O.filt, idents = sds.data$Sample)


  genes.expr <- as.matrix(S.O.filt@assays$RNA@data)
  #genes.expr <- as.matrix(S.O@assays$RNA@data)
  

  genes.df <- data.frame(GeneID = rownames(genes.expr),
                                    genes.expr) %>%
    pivot_longer(-c(GeneID), names_to = 'Sample', values_to = 'log2.expr')

  if(method == 'atac'){
    genes.df$Sample <- gsub('\\.', '-', genes.df$Sample)
  }

  genes.df <- inner_join(genes.df, sds.data, by = 'Sample')
 
  sds.data <- as.data.frame(sds.data)
  rownames(sds.data) <- sds.data$Sample
  
  ## Add the new clusters as meta-data
  S.O.filt <- AddMetaData(S.O.filt, sds.data)
  #S.O <- AddMetaData(S.O, sds.data)

  L <- list(sds.data = sds.data, genes.df = genes.df,
            S.O = S.O.filt)
  
  return(L)
}


DefaultAssay(atac_sub) <- 'RNA'
L.atac <- fitTime(atac_sub, 'atac', reverse.t = T)

DefaultAssay(rna_sub) <- 'RNA'
L.rna <- fitTime(rna_sub, 'rna', reverse.t = T)


saveRDS(L.rna$genes.df, '../Input/toxo_cdc/rds/sc_rna_genes_expr_pt.rds')
saveRDS(L.atac$genes.df, '../Input/toxo_cdc/rds/sc_atac_genes_expr_pt.rds')

saveRDS(L.rna$sds.data, '../Input/toxo_cdc/rds/sc_rna_sds_data.rds')
saveRDS(L.atac$sds.data, '../Input/toxo_cdc/rds/sc_atac_sds_data.rds')

saveRDS(L.rna$S.O, '../Input/toxo_cdc/rds/S.O_intra_lables_pt.rds')
saveRDS(L.atac$S.O, '../Input/toxo_cdc/rds/S.O_intra_atac_lables_pt.rds')


