library(tidyverse)
library(openxlsx)
library(doParallel)
library(splines)

## Read in the data.
sc.rna.cell.cycle.genes.df.adj  <- readRDS('../Input/toxo_cdc/rds/rna_seq_gene_expr_pseudo_time.rds')
sc.atac.cell.cycle.genes.df.adj <- readRDS('../Input/toxo_cdc/rds/atac_seq_gene_access_pseudo_time.rds')

## Common genes between the data sets
comm.genes <- intersect(unique(sc.rna.cell.cycle.genes.df.adj$GeneID), unique(sc.atac.cell.cycle.genes.df.adj$GeneID))

getCurvePeakLoc <- function(t, y, prob = 0.8){
  
  ## Fitting the estimated kernel with smooth splines
  spline.fit <- smooth.spline(x = t, y = y)
  
  
  ## Compute the derivatives of the fitted splines
  s.0 <- predict(spline.fit, spline.fit$x, deriv=0)
  s.1 <- predict(spline.fit, spline.fit$x, deriv=1)
  s.derv <- data.frame(s0=s.0$y, s1=s.1$y)
  
  ## Get the location of the extrema
  locs <- rle(den.sign <- sign(s.derv$s1))
  
  
  ## Maxima
  inc.ind <- which(locs$values == 1)
  if(length(inc.ind) > 1){
    maxima.ind = {}
    for(i in inc.ind){
      maxima.ind = c(maxima.ind, sum(locs$lengths[1:i]))
    }
    ## Interpolate a point between the location where derivative changes sign
    maxima = (spline.fit$x[maxima.ind] + spline.fit$x[(maxima.ind + 1)]) / 2
    maxima = maxima[!is.na(maxima)]
    ## Get the maximum values
    maxval = predict(spline.fit, maxima)
    
    ## Get the outliers
    maxima.outliers = which(maxval$y >= quantile(maxval$y, prob = prob))
    
    ## Peaks for entities of interest
    entity.x = maxval$x[maxima.outliers]
    entity.y = maxval$y[maxima.outliers]
  }else{
    entity.x <- spline.fit$x[which.max(spline.fit$y)]
    entity.y <- spline.fit$y[which.max(spline.fit$y)]
  }
  
  return(entity.x)
}


fitPsplines <- function(t, y){
  # B-spline matrix
  ## Number of internal knots for B-spline basis 
  niknots <- 3
  
  ## Smoothing parameters for spline smoothing
  sparx <- c(5, 1, 1)
  
  knots <- c(rep(0,3), seq(0, 6, len = niknots), rep(6, 3)) 
  gridout <- seq(from = 0, to = 6, by = 1/3)
  gridsize <- length(gridout)
  Bout <- splineDesign(knots, gridout)
  nknots <- ncol(Bout)
  
  ## Penalty matrices
  v1 <- Bout[gridsize,] - Bout[1,]
  v1 <- v1 / sqrt(sum(v1^2))
  DB <- splineDesign(knots, c(0, 6), derivs = 1L)
  v2 <- DB[2,] - DB[1,]
  v2 <- v2 / sqrt(sum(v2^2))
  D2B <- splineDesign(knots, gridout, derivs = 2L)
  G <- crossprod(D2B) 
  G <- G / norm(G, "2")
  PENx <- sparx[1] * tcrossprod(v1) + 
    sparx[2] * tcrossprod(v2) + sparx[3] * G

  ## P-spline
  Bx <- splineDesign(knots, t)
  BtB <- crossprod(Bx)
  nrmBtB <- norm(BtB, "2")
  ## Smoothing
  dim(y) <- c(length(y), 1)
  xcoef <- solve(BtB + nrmBtB * PENx, t(Bx)) %*% y
  
  xs <- Bout %*% xcoef
  
  # Calculate first and second derivatives
  DB <- splineDesign(knots, gridout, derivs = 1L)
  D2B <- splineDesign(knots, gridout, derivs = 2L)
  Dxs <- DB %*% xcoef
  D2xs <- D2B %*% xcoef
  
  s.out <- data.frame(x = gridout, y = xs, yp = Dxs, ypp = D2xs)
  
  return(s.out)
}


smoothFilter <- function(y){
  window.length <- floor(length(y) / 30)
  filt.y <- rep(F, length(y))
  for(i in 1:30){
    # d <- density(y[((i - 1) * window.length + 1):(i * window.length)])
    # plot(d)
    # 
    seg <- ((i - 1) * window.length + 1):(i * window.length)
    prop.non.zero <- sum(y[seg] > 0) / 
      length(seg)
    
    print(prop.non.zero)
    if(prop.non.zero > 0.1){
      filt.y[seg] <- T
    }
    #Sys.sleep(0.8)
  }
  
 return(which(y == 0 & filt.y))
}

## Fit smoothing splines and sample at regular time-points

## Expression
sc.rna.spline.fits <- mclapply(1:length(comm.genes), function(i){
  tmp <- sc.rna.cell.cycle.genes.df.adj %>% dplyr::filter(GeneID == comm.genes[i]) %>%
    transmute(GeneID = GeneID, x = adj.pt, y = log2.expr)
  sc.rna.sp <- fitPsplines(tmp$x, tmp$y)
  mu <- data.frame(x = sc.rna.sp$x, y = sc.rna.sp$y) 
  mu <- data.frame(GeneID = rep(tmp$GeneID[1], length(mu[,1])), x = mu[,1], y = mu[,2])
  return(mu)
}, mc.cores = num.cores)

sc.rna.spline.fits <- bind_rows(sc.rna.spline.fits)

saveRDS(sc.rna.spline.fits, '../Input/toxo_cdc/rds/sc_rna_p_spline_fits.rds')

## Access
sc.atac.spline.fits <- mclapply(1:length(comm.genes), function(i){
  tmp <- sc.atac.cell.cycle.genes.df.adj %>% dplyr::filter(GeneID == comm.genes[i]) %>%
    transmute(GeneID = GeneID, x = adj.time, y = log2.expr)
  sc.atac.sp <- fitPsplines(tmp$x, tmp$y)
  mu <- data.frame(x = sc.atac.sp$x, y = sc.atac.sp$y)
  mu <- data.frame(GeneID = rep(tmp$GeneID[1], length(mu[,1])), x = mu[,1], y = mu[,2])
  return(mu)
}, mc.cores = num.cores)

sc.atac.spline.fits <- bind_rows(sc.atac.spline.fits)

saveRDS(sc.atac.spline.fits, '../Input/toxo_cdc/rds/sc_atac_p_spline_fits.rds')

