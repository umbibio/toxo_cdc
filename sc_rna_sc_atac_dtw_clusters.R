library(dtwclust)
library(doParallel)
library(bigmemory)
library(tidyverse)
library(openxlsx)


## For parallel calculations
num.cores <- detectCores(all.tests = FALSE, logical = TRUE)

getCurvePeakLoc <- function(t, y){
  
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
    maxima.outliers = which(maxval$y >= quantile(maxval$y, prob = 0.8))
    
    ## Peaks for entities of interest
    entity.x = maxval$x[maxima.outliers]
    entity.y = maxval$y[maxima.outliers]
  }else{
    entity.x <- spline.fit$x[which.max(spline.fit$y)]
    entity.y <- spline.fit$y[which.max(spline.fit$y)]
  }
  
  return(entity.x)
}


sc.rna.spline.fits <- readRDS('../Input/toxo_cdc/rds/sc_rna_spline_fits_cell_cycle_genes.rds')
sc.atac.spline.fits <- readRDS('../Input/toxo_cdc/rds/sc_atac_spline_fits_cell_cycle_genes.rds')


## Clustering using Dynamic Time Warping (independently done for Expression and Access)
dtwClustCurves <- function(tc.mus, nclust = 6L){
  ## Calculate clusters in parallel
  num.cores <- detectCores(all.tests = FALSE, logical = TRUE)
  workers <- makeCluster(num.cores)
  invisible(clusterEvalQ(workers, library("dtwclust")))
  registerDoParallel(workers)
  tc.mus <- lapply(1:ncol(tc.mus), function(i) c(as.numeric(tc.mus[,i])))
  hc_dtw <- tsclust(tc.mus, 
                    type = "h", 
                    k = nclust, 
                    distance = "dtw", 
                    control = hierarchical_control(method = "complete"),
                    centroid = shape_extraction, 
                    #preproc = NULL, 
                    preproc = zscore,
                    trace = T,
                    args = tsclust_args(dist = list(window.size = 4L))
  )
  
  stopCluster(workers)
  registerDoSEQ()
  
  return(hc_dtw)
}


## Turn the data into wide format (time by gene) and center & scale each gene
sc.rna.dtw.wide <- sc.rna.spline.fits %>% 
  pivot_wider(names_from = 'GeneID', values_from = 'y') %>% 
  mutate_at(vars(matches("TGME")), scale) %>%
  as.data.frame()

sc.atac.dtw.wide <- sc.atac.spline.fits %>% 
  pivot_wider(names_from = 'GeneID', values_from = 'y') %>% 
  mutate_at(vars(matches("TGME")), scale) %>%
  as.data.frame()


## Set cluster numbers
num.clust <- 8L

sc.rna.hc_dtws <- dtwClustCurves(sc.rna.dtw.wide[2:ncol(sc.rna.dtw.wide)], nclust = num.clust)

plot(sc.rna.hc_dtws, type = 'sc')
plot(sc.rna.hc_dtws, type = "series", clus = 1L)
plot(sc.rna.hc_dtws, type = "centroids", clus = 1L)


sc.atac.hc_dtws <- dtwClustCurves(sc.atac.dtw.wide[2:ncol(sc.atac.dtw.wide)], nclust = num.clust)

plot(sc.atac.hc_dtws, type = 'sc')
plot(sc.atac.hc_dtws, type = "series", clus = 1L)
plot(sc.atac.hc_dtws, type = "centroids", clus = 1L)



## Get the cluster info
rna.clust.info <- data.frame(GeneID = colnames(sc.rna.dtw.wide)[2:ncol(sc.rna.dtw.wide)], 
                             cluster = cutree(sc.rna.hc_dtws, k = num.clust))

atac.clust.info <- data.frame(GeneID = colnames(sc.atac.dtw.wide)[2:ncol(sc.atac.dtw.wide)], 
                             cluster = cutree(sc.atac.hc_dtws, k = num.clust))


## Join the two:
rna.atac.clust.info <- inner_join(rna.clust.info, atac.clust.info, by = 'GeneID')
colnames(rna.atac.clust.info) <- c('GeneID', 'cluster.rna', 'cluster.atac')

table(rna.atac.clust.info$cluster.rna, rna.atac.clust.info$cluster.atac)
