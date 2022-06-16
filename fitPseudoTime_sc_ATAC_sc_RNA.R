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
Idents(S.O.integrated) <- 'orig.ident'

atac_sub <- subset(S.O.integrated, ident = 'scATAC')
rna_sub <- subset(S.O.integrated, ident = 'scRNA')


## IDs
prod.desc  <- read.xlsx('../Input/toxo_genomics/genes/ProductDescription_GT1.xlsx')
TGGT1_ME49 <- read.xlsx('../Input/toxo_genomics/Orthologs/TGGT1_ME49 Orthologs.xlsx')
prod.desc <- left_join(prod.desc, TGGT1_ME49, by = c('GeneID' = 'TGGT1'))


## Pseudo-time analysis with SLingshot

fitTime <- function(S.O, method = 'rna'){
  pc.tg <- getPCA(S.O)
  sds.data <- getPrinCurve(pc.tg)
  pc.sds.tg <- left_join(pc.tg, sds.data, by = "Sample")
  pc.sds.tg$phase <- S.O@meta.data$phase[match(pc.sds.tg$Sample, rownames(S.O@meta.data))]
  pc.sds.tg$phase <- factor(pc.sds.tg$phase, levels = c('G1.a', 'G1.b', 'S', 'M', 'C'))
  
  Y <- log2(S.O@assays$RNA@data + 1)
  var.genes <- names(sort(apply(Y, 1, var),decreasing = TRUE))#[1:1000] 
  Y <- Y[var.genes, ]
  
  pt <- sds.data$pt
  
  ## Map the pseudo-time to 0-6:20 hours 
  #t <- (12 + 1/3) * ((as.numeric(pt) - min(as.numeric(pt)))/(max(as.numeric(pt)) - min(as.numeric(pt))))
  t <- (6 + 1/6) * ((as.numeric(pt) - min(as.numeric(pt)))/(max(as.numeric(pt)) - min(as.numeric(pt))))
  
  sds.data$t <- t
  
  ## Update pseudo-time to 0:6h
  sds.data$pt <- t <- 6 * ((as.numeric(pt) - min(as.numeric(pt)))/(max(as.numeric(pt)) - min(as.numeric(pt))))
  
  ## time-index cells in 20 min intervals and identify cells in each partition
  ## They will be considered as replicates
  #time.breaks <- seq(1/3, 12 + 1/3, by = 1/3) 
  time.breaks <- seq(1/6, 6 + 1/6, by = 1/6) 
  time.idx <- rep(0, nrow(sds.data))
  
  ind <- which(sds.data$t <= time.breaks[1])
  time.idx[ind] <- 0
  
  for(i in 2:length(time.breaks)){
    ind <- which(sds.data$t > time.breaks[(i-1)] & sds.data$t <= time.breaks[i])
    time.idx[ind] <- i - 1
  }
  
  sds.data$time.idx <- time.idx
  
  ## Update the time to 20 min increments
  #sds.data$t <- (time.idx) * (1/3)
  sds.data$t <- (time.idx) * (1/6)
  
  # sds.data <- sds.data %>%  
  #   group_by(time.idx) %>% mutate(rep = seq(1:n()))
  # 
  # rownames(sds.data) <- sds.data$Sample
  
  
  
  ## Run a GAM regression of expression on the pseudo-time
  ## Use parallel computation to speed things up. 16 cores
  gam.pval <- mclapply(1:nrow(Y), function(z){
    d <- data.frame(z=as.numeric(Y[z,]), t=as.numeric(pt))
    tmp <- gam(z ~ lo(t), data=d)
    #p <- summary(tmp)[4][[1]][1,5]
    p <- summary(tmp)$anova$`Pr(F)`[2]
    p
  }, mc.cores = num.cores)
  
  gam.pval <- unlist(gam.pval)
  names(gam.pval) <- rownames(Y)
  ## Remove the NA's and get the best fits
  if(any(is.na(gam.pval))){
    gam.pval <- gam.pval[-which(is.na(gam.pval))]
  }
  
  gam.pval.adj <- p.adjust(gam.pval, method = 'fdr', n = length(gam.pval))
  gam.pval.sig <- gam.pval[gam.pval.adj < 0.01] 
  print(length(gam.pval.sig)) ## number of correlating genes
  
  ## Sort the cells on the pt
  cell.ord <- sds.data$cell.ord
  
  topgenes <- names(sort(gam.pval.sig, decreasing = FALSE))  
  cell.cycle.genes.expr <- as.matrix(S.O@assays$RNA@data[topgenes, cell.ord])
  #cell.cycle.genes.expr <- as.matrix(S.O.bd.filt@assays$smooth@data[topgenes, cell.ord]) ## smoothed version
  
  
  cell.cycle.genes.df <- data.frame(GeneID = rownames(cell.cycle.genes.expr),
                                    cell.cycle.genes.expr) %>% 
    pivot_longer(-c(GeneID), names_to = 'Sample', values_to = 'log2.expr')
  
  if(method == 'atac'){
    cell.cycle.genes.df$Sample <- gsub('\\.', '-', cell.cycle.genes.df$Sample)
  }
  
  cell.cycle.genes.df$GeneID <- gsub('-', '_', cell.cycle.genes.df$GeneID)
  #cell.cycle.genes.df <- left_join(cell.cycle.genes.df, sds.data, by = 'Sample')
  cell.cycle.genes.df$cluster <- S.O@meta.data$seurat_clusters[match(cell.cycle.genes.df$Sample, 
                                                                        rownames(S.O@meta.data))]
  
  cell.cycle.genes.df$phase <- S.O@meta.data$phase[match(cell.cycle.genes.df$Sample, 
                                                            rownames(S.O@meta.data))]
  cell.cycle.genes.df$phase <- factor(cell.cycle.genes.df$phase, 
                                      levels = c('G1.a', 'G1.b', 'S', 'M', 'C'))
  
  
  
  tmp <- left_join(cell.cycle.genes.df, sds.data, by = 'Sample')
  start.times <- tmp %>% group_by(GeneID) %>% arrange(phase, pt) %>% summarise(st = t[1])
  lag.time <- which(sort(unique(tmp$t)) == unique(start.times$st)) + 2 ## all are the same
  
  print(paste('lag time:', lag.time))
  
  
  adjusted.time <- (sds.data$time.idx * 1/6) -  sort(unique(sds.data$time.idx) * 1/6)[lag.time]
  neg.ind <- ifelse(adjusted.time < 0, T, F)
  adjusted.time[neg.ind] <- adjusted.time[neg.ind] + (6 + 1/6)
  sds.data$adj.time <-  adjusted.time
  
  
  start.times <- tmp %>% group_by(GeneID) %>% arrange(phase, pt) %>% summarise(spt = pt[1])
  lag.time <- which(sort(unique(tmp$pt)) == unique(start.times$spt))
  
  adjusted.pstime <- sds.data$pt -  sort(unique(sds.data$pt))[lag.time]
  neg.ind <- ifelse(adjusted.pstime < 0, T, F)
  adjusted.pstime[neg.ind] <- adjusted.pstime[neg.ind] + 6
  sds.data$adj.pt <-  adjusted.pstime
  
  
  ## create fine-resolution 20 min cluster of cells
  clusters <- paste('C', 1:length(unique(sds.data$adj.time)), sep = '')
  sds.data$cluster <- clusters[as.integer((sds.data$adj.time) * 6 + 1)]
  
  
  ## Updat time index
  time.breaks <- seq(1/6, 6 + 1/6, by = 1/6) 
  time.idx <- rep(0, nrow(sds.data))
  
  ind <- which(sds.data$adj.time <= time.breaks[1])
  time.idx[ind] <- 0
  
  for(i in 2:length(time.breaks)){
    ind <- which(sds.data$adj.time > time.breaks[(i-1)] & sds.data$adj.time <= time.breaks[i])
    time.idx[ind] <- i - 1
  }
  
  sds.data$adj.time.idx <- time.idx
  

  
  sds.data <- sds.data %>%  ungroup() %>%
    group_by(adj.time.idx) %>% mutate(rep = seq(1:n()))
  
  sds.data <- as.data.frame(sds.data)
  rownames(sds.data) <- sds.data$Sample
  
  ## Add the new clusters as meta-data
  S.O <- AddMetaData(S.O, sds.data)
  
  pc <- S.O[['pca']]@cell.embeddings
  pc <- data.frame(pc) %>% dplyr::mutate(Sample = rownames(pc))
  pc$cluster <- S.O$cluster
  
  
  
  pc.sds.adj <- left_join(pc, sds.data, by = "Sample")
  
  lvs <- paste('C', unique(sort(as.numeric(gsub('C', '', pc.sds.adj$cluster.y)))), sep = '')
  pc.sds.adj$cluster.y <- factor(pc.sds.adj$cluster.y, levels = lvs)
  
  
  cell.cycle.genes.df.adj <- left_join(cell.cycle.genes.df, sds.data[,c('Sample','adj.pt', 'adj.time', 
                                                                        'adj.time.idx', 
                                                                        'rep', 'cluster')], by = 'Sample')
  sc.tc.df.adj <- cell.cycle.genes.df.adj %>% 
    transmute(y = log2.expr, tme = adj.time, ind = rep, variable = GeneID)
  
  L <- list(cell.cycle.genes.df = cell.cycle.genes.df, 
            cell.cycle.genes.df.adj = cell.cycle.genes.df.adj,
            sc.tc.df.adj = sc.tc.df.adj,
            sds.data = sds.data,
            S.O = S.O,
            gam.genes = names(gam.pval.sig))
  return(L)
}


DefaultAssay(atac_sub) <- 'RNA'
L.atac <- fitTime(atac_sub, 'atac')

DefaultAssay(rna_sub) <- 'RNA'
L.rna <- fitTime(rna_sub)


## Scale the pt based on know biology: Radke et. al 2000
scalePT <- function(pt.phase){
  G <- c(0, 3) # 3h
  S <- c(3, 4.7) # 1.7h
  M <- c(4.7, 5) # ~20 min
  C <- c(5, 6) # 1h
  
  G.adj.pt <- c(min(pt.phase$pt[pt.phase$phase %in% c('G1.a', 'G1.b')]),
                max(pt.phase$pt[pt.phase$phase %in% c('G1.a', 'G1.b')]))
  
}



saveRDS(L.atac$cell.cycle.genes.df.adj, '../Input/toxo_cdc/rds/atac_seq_gene_access_pseudo_time.rds')
saveRDS(L.rna$cell.cycle.genes.df.adj, '../Input/toxo_cdc/rds/rna_seq_gene_expr_pseudo_time.rds')


# Scale the pt based on know biology: Radke et. al 2000
scalePT <- function(pt.phase){
  G <- c(0, 3) # 3h
  S <- c(3, 4.7) # 1.7h
  M <- c(4.7, 5) # ~20 min
  C <- c(5, 6) # 1h
  
  G.adj.pt <- c(min(pt.phase$pt[pt.phase$phase %in% c('G1.a', 'G1.b')]),
                max(pt.phase$pt[pt.phase$phase %in% c('G1.a', 'G1.b')]))
  
}

pt.phase <- sc.rna.cell.cycle.genes.df.adj %>% transmute(Sample = Sample, pt = adj.pt, phase = phase) %>% distinct()

par(mfrow = c(4,1))
plot(pt.phase$pt[pt.phase$phase %in% c('G1.a', 'G1.b')])
plot(pt.phase$pt[pt.phase$phase == 'S'])
plot(pt.phase$pt[pt.phase$phase == 'M'])
plot(pt.phase$pt[pt.phase$phase == 'C'])

