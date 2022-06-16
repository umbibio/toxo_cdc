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
library(dtwclust)
library(doParallel)
library(tidyverse)
library(tidytext)
library(tvReg)



## No time bin projections
sc.rna.cell.cycle.genes.df.adj  <- readRDS('../Output/scClockOut/forDavid/rna_seq_gene_expr_pseudo_time.RData')
sc.atac.cell.cycle.genes.df.adj <- readRDS('../Output/scClockOut/forDavid/atac_seq_gene_access_pseudo_time.RData')

comm.genes <- intersect(unique(sc.rna.cell.cycle.genes.df.adj$GeneID), unique(sc.atac.cell.cycle.genes.df.adj$GeneID))

sc.rna.spline.fits <- lapply(1:length(comm.genes), function(i){
  tmp <- sc.rna.cell.cycle.genes.df.adj %>% dplyr::filter(GeneID == comm.genes[i]) %>%
    transmute(GeneID = GeneID, x = adj.time, y = log2.expr)
  sc.rna.sp <- smooth.spline(tmp$x, tmp$y)
  sc.rna.sp <- predict(sc.rna.sp, seq(0, 6, by = 1/3)) 
  mu <- data.frame(x = sc.rna.sp$x, y = scale(sc.rna.sp$y)) ## For regression do not scale
  mu <- data.frame(GeneID = rep(tmp$GeneID[1], length(mu[,1])), x = mu[,1], y = mu[,2])
  return(mu)
})

sc.rna.spline.fits <- bind_rows(sc.rna.spline.fits)

sc.atac.spline.fits <- lapply(1:length(comm.genes), function(i){
  tmp <- sc.atac.cell.cycle.genes.df.adj %>% dplyr::filter(GeneID == comm.genes[i]) %>%
    transmute(GeneID = GeneID, x = adj.time, y = log2.expr)
  sc.atac.sp <- smooth.spline(tmp$x, tmp$y)
  sc.atac.sp <- predict(sc.atac.sp, seq(0, 6, by = 1/3)) 
  mu <- data.frame(x = sc.atac.sp$x, y = sc.atac.sp$y) %>% scale()
  mu <- data.frame(GeneID = rep(tmp$GeneID[1], length(mu[,1])), x = mu[,1], y = mu[,2])
  return(mu)
})

sc.atac.spline.fits <- bind_rows(sc.atac.spline.fits)

sc.rna.atac.spline.fits <- bind_cols(sc.rna.spline.fits, sc.atac.spline.fits[,2:3])
colnames(sc.rna.atac.spline.fits) <- c('GeneID', 'x1', 'y1', 'x2', 'y2')

sc.rna.atac.corr <- sc.rna.atac.spline.fits %>% group_by(GeneID) %>% summarise(m.ex = mean(y1), m.ax = mean(y2), co = cor(y1, y2))

par(mfrow = c(1,1))
for(i in 1:10){
  tmp <- sc.rna.atac.spline.fits %>% dplyr::filter(GeneID == unique(sc.rna.atac.spline.fits$GeneID)[i])
  y.lim = c(min(c(tmp$y1, tmp$y1)), max(c(tmp$y2, tmp$y2)))
  plot(tmp$x1, tmp$y1, type = 'l', col = 'red', main = unique(sc.rna.atac.spline.fits$GeneID)[i], ylim = y.lim)
  points(tmp$x2, tmp$y2, type = 'l', col = 'blue')
  
  # plot(tmp$x1, tmp$y1, type = 'l', col = 'red', main = paste('EXPR:', unique(sc.rna.atac.spline.fits$GeneID)[i]))
  # plot(tmp$x2, tmp$y2, type = 'l', col = 'blue', main = paste('ATAC:', unique(sc.rna.atac.spline.fits$GeneID)[i]))
  # plot(tmp$x2, tmp$y1 - tmp$y2, type = 'l', col = 'green', main = paste('Delay:', unique(sc.rna.atac.spline.fits$GeneID)[i]))
  Sys.sleep(0.9)
}


####### FUNCTIONAL REGRESSION: DO NOT DO
## Do a regression (from David)
tgrid <- seq(0, 6, by = 1/3)
gridsize <- length(tgrid) 

feat <- array(dim = c(gridsize, 4, 50))
dimnames(feat) <- list(time = 1:gridsize, 
                       coef = c("intercept", "slope", "y.expr", "y.fitted"), 
                       gene = unique(sc.rna.atac.spline.fits$GeneID)[1:50])

one <- rep(1, gridsize) # constant function
bw <- 2 # smoothing bandwidth, hand-selected
for (i in 1:50) {
  tmp <- sc.rna.atac.spline.fits %>% dplyr::filter(GeneID == unique(sc.rna.atac.spline.fits$GeneID)[i])
  fit <- tvOLS(x = cbind(one, tmp$y2), y = tmp$y1, 
               z = tgrid, bw = bw, est = "ll")	
  y.fitted <- apply(cbind(one, tmp$y2)*fit$coefficients, 1, sum)
  feat[,,i] <- c(fit$coefficients, tmp$y1, y.fitted)

} 


par(mfrow = c(3,1))
for(i in 1:50){
  tmp <- feat[,,i]
  plot(tgrid,tmp[,3], type = 'l', col = 'red', main = paste('EXPR:', dimnames(feat)$gene[i]))
  plot(tgrid,tmp[,4], type = 'l', col = 'blue', main = paste('Fitted:', dimnames(feat)$gene[i]))
  plot(tgrid,tmp[,2], type = 'l', col = 'green', main = paste('slope:', dimnames(feat)$gene[i]))
  Sys.sleep(0.9)
}
####### FUNCTIONAL REGRESSION: DO NOT DO

## Make sure expression is scaled
sc.rna.dtw.wide <- sc.rna.spline.fits %>% 
  pivot_wider(names_from = 'GeneID', values_from = 'y') %>%
  as.data.frame()

sc.atac.dtw.wide <- sc.atac.spline.fits %>% 
  pivot_wider(names_from = 'GeneID', values_from = 'y') %>%
  as.data.frame()



## Generate the clusters
num.clust <- 4L

sc.rna.hc_dtws <- dtwClustCurves(sc.rna.dtw.wide[2:ncol(sc.rna.dtw.wide)], nclust = num.clust)
plot(sc.rna.hc_dtws, type = 'sc')
plot(sc.rna.hc_dtws, type = "series", clus = 1L)
plot(sc.rna.hc_dtws, type = "centroids", clus = 1L)


sc.atac.hc_dtws <- dtwClustCurves(sc.atac.dtw.wide[2:ncol(sc.atac.dtw.wide)], nclust = num.clust)
plot(sc.atac.hc_dtws, type = 'sc')
plot(sc.atac.hc_dtws, type = "series", clus = 1L)
plot(sc.atac.hc_dtws, type = "centroids", clus = 1L)


TT <- table(sc.atac.hc_dtws@cluster, sc.rna.hc_dtws@cluster)
TT
sum(diag(TT))/sum(TT)



### Regression model: From David:
##################################
# Functon to fit regression model 
# y(t) = a + b * x(t + c) + e(t) 
# by ordinary least squares
# with x and y periodic functions
##################################


align.xy <- function(x, y)
{
  n <- length(x)
  xbar <- mean(x)
  ybar <- mean(y)
  
  ## Center x and y
  x <- x - xbar
  y <- y - ybar
  
  ## Squared norms
  sqnrmx <- sum(x^2)
  sqnrmy <- sum(y^2)
  if (sqnrmx == 0 || sqnrmy == 0) 
    return(list(a = ybar, b = 0, c = 0, xshifted = x))
  
  ## Create shifted versions of x
  shift <- 1:n
  idx <- sapply(1:n, seq.int, length.out = n)
  xx <- c(x,x)
  xshift <- matrix(xx[idx], nrow = n)
  
  ## Estimate shift
  cp <- crossprod(y, xshift) 
  idx <- which.max(cp)
  
  ## Estimate intercept and scale
  b <- as.numeric(cp[idx]) / sqnrmx
  a <- ybar - b * xbar
  
  ## Goodness of fit
  Rsquared <- b^2 * sqnrmx / sqnrmy 
  
  list(a = a, b = b, c = idx - 1, xshifted = xshift[,idx] + xbar, 
       Rsquared = Rsquared)
}

### Test the function
# x <- sin(seq(0, pi, len=101))
# y <- 5 + 2 * x
# test <- estimate.lin.transfo(x, y)
# print(test)

## xs and ys are accessibility, expression matricies (fitted splies, with same points)
n <- ncol(xs) # number of genes
T <- nrow(xs) # number of time points
# i <- sample.int(n,1)
# matplot(tgrid, cbind(xs[,i], ys[,i]), type = "l")
out <- lapply(1:n, function(i) align.xy(xs[,i], ys[,i]))
xshift <- sapply(out, "[[", "xshifted") 
a <- sapply(out, "[[", "a")
b <- sapply(out, "[[", "b")
c <- sapply(out, "[[", "c")
# Reformulate shift in [-T/2,T/2] instead of [0,T] 
c[c > (T/2)] <- c[c > (T/2)] - T
Rsquared <- sapply(out, "[[", "Rsquared")

hist(c*6/T, main = "Distribution of shifts (Period = 6)", xlab = "Time shift")
summary(Rsquared)

