library(tidyverse)
library(openxlsx)
library(doParallel)
library(npreg)

## For parallel calculations
num.cores <- detectCores(all.tests = FALSE, logical = TRUE)


## Read in the data.
sc.rna.cell.cycle.genes.df.adj  <- readRDS('../Input/toxo_cdc/rds/rna_seq_gene_expr_pseudo_time.rds')
sc.atac.cell.cycle.genes.df.adj <- readRDS('../Input/toxo_cdc/rds/atac_seq_gene_access_pseudo_time.rds')


sc.rna.genes.expr.pt.df <- readRDS('../Input/toxo_cdc/rds/sc_rna_genes_expr_pt.rds')
sc.atac.genes.expr.pt.df <- readRDS('../Input/toxo_cdc/rds/sc_atac_genes_expr_pt.rds')

filter.low.var <- function(genes.expr.pt.df){
  genes.var <- genes.expr.pt.df %>% group_by(GeneID) %>% 
    summarise(mean.expr = mean(log2.expr), var.expr = var(log2.expr))
  q.ex <- quantile(genes.var$var.expr, prob = 0.5)
  gene.ex <- genes.var$GeneID[which(genes.var$var.expr <= q.ex)]
  genes.expr.pt.df.filt <- genes.expr.pt.df %>% filter(!(GeneID %in% gene.ex))
  
  return(genes.expr.pt.df.filt)
}

sc.rna.genes.expr.pt.df.filt <- filter.low.var(sc.rna.genes.expr.pt.df)
sc.atac.genes.expr.pt.df.filt <- filter.low.var(sc.atac.genes.expr.pt.df)

## Common genes between the data sets
comm.genes <- intersect(unique(sc.rna.genes.expr.pt.df.filt$GeneID), unique(sc.atac.genes.expr.pt.df.filt$GeneID))

## Filtering 0's by sliding thresholding
smoothFilter <- function(y){
  window.length <- floor(length(y) / 30)
  filt.y <- rep(F, length(y))
  for(i in 1:30){
    #vd <- density(y[((i - 1) * window.length + 1):(i * window.length)])
    #vplot(d)
    # 
    seg <- ((i - 1) * window.length + 1):(i * window.length)
    prop.non.zero <- sum(y[seg] > 0) / 
      length(seg)
    
    #print(paste(t[seg][1], t[seg][length(seg)]))
    print(prop.non.zero)
    
    if(prop.non.zero > 0.05){
      filt.y[seg] <- T
    }
    #Sys.sleep(0.8)
  }
  
  return(which(y == 0 & filt.y))
}

do.filt <- F
## Fit smoothing splines and sample at regular time-points

## Expression
sc.rna.spline.fits <- mclapply(1:length(comm.genes), function(i){
  tmp <- sc.rna.genes.expr.pt.df.filt %>% dplyr::filter(GeneID == comm.genes[i]) %>%
    transmute(GeneID = GeneID, x = pt.shifted.scaled, y = log2.expr)
  
  pt.phase <- sc.rna.genes.expr.pt.df.filt %>% transmute(pt = pt.shifted.scaled, phase = phase)
  
  if(do.filt){
    ind <- smoothFilter(tmp$y)
    tmp2 <- tmp[-ind,]
    y <- tmp2$y
    t <- tmp2$x
  }else{
    y <- tmp$y
    t <- tmp$x
    
  }
  sc.rna.sp <- ss(t, y, periodic = T, lambda = 0.0001)
  #sc.rna.sp <- smooth.spline(t, y, lambda = 0.1)
  sc.rna.sp <- predict(sc.rna.sp, seq(0, 6, by = 1/3)) 
  #sc.rna.pp <- fitPsplines(t, y)
  #plot(tmp$x, tmp$y)
  #points(t, y, col = 'red')
  #points(sc.rna.sp$x, sc.rna.sp$y, type = 'l', col = 'red')
  #points(sc.rna.sp$x, sc.rna.sp$y, type = 'l', col = 'blue')
  #points(sc.rna.pp$x, sc.rna.sp$y, type = 'l', col = 'green')
  mu <- data.frame(x = sc.rna.sp$x, y = sc.rna.sp$y) 
  mu <- data.frame(GeneID = rep(tmp$GeneID[1], length(mu[,1])), x = mu[,1], y = mu[,2])
  return(mu)
}, mc.cores = num.cores)

sc.rna.spline.fits <- bind_rows(sc.rna.spline.fits)

saveRDS(sc.rna.spline.fits, '../Input/toxo_cdc/rds/sc_rna_ss_spline_fits.rds')

## Access
sc.atac.spline.fits <- mclapply(1:length(comm.genes), function(i){
  tmp <- sc.atac.genes.expr.pt.df.filt %>% dplyr::filter(GeneID == comm.genes[i]) %>%
    transmute(GeneID = GeneID, x = pt.shifted.scaled, y = log2.expr)
  
  pt.phase <- sc.atac.genes.expr.pt.df.filt %>% transmute(pt = pt.shifted.scaled, phase = phase)
  
  if(do.filt){
    ind <- smoothFilter(tmp$y)
    tmp2 <- tmp[-ind,]
    y <- tmp2$y
    t <- tmp2$x
  }else{
    y <- tmp$y
    t <- tmp$x
    
  }
  sc.atac.sp <- ss(t, y, periodic = T, lambda = 0.0001)
  #sc.atac.sp <- smooth.spline(t, y, lambda = 0.1)
  sc.atac.sp <- predict(sc.atac.sp, seq(0, 6, by = 1/3)) 
 
  mu <- data.frame(x = sc.atac.sp$x, y = sc.atac.sp$y)
  mu <- data.frame(GeneID = rep(tmp$GeneID[1], length(mu[,1])), x = mu[,1], y = mu[,2])
  return(mu)
}, mc.cores = num.cores)

sc.atac.spline.fits <- bind_rows(sc.atac.spline.fits)

saveRDS(sc.atac.spline.fits, '../Input/toxo_cdc/rds/sc_atac_ss_spline_fits.rds')

