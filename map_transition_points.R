library(tidyverse)
library(openxlsx)
library(Seurat)
library(ggplot2)
library(ggfortify)
library(jcolors)
require(gridExtra)
library(grid)
library(tidytext)




## For parallel calculations
num.cores <- detectCores(all.tests = FALSE, logical = TRUE)

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


sc.rna.spline.fits <- readRDS('../Input/toxo_cdc/rds/sc_rna_spline_fits.rds')
sc.atac.spline.fits <- readRDS('../Input/toxo_cdc/rds/sc_atac_spline_fits.rds')

## Turn the data into wide format (time by gene) and center & scale each gene
sc.rna.dtw.wide <- sc.rna.spline.fits %>% 
  pivot_wider(names_from = 'GeneID', values_from = 'y') %>% 
  mutate_at(vars(matches("TGME")), scale) %>%
  as.data.frame()

sc.atac.dtw.wide <- sc.atac.spline.fits %>% 
  pivot_wider(names_from = 'GeneID', values_from = 'y') %>% 
  mutate_at(vars(matches("TGME")), scale) %>%
  as.data.frame()


sc.rna.mu.scale <- sc.rna.dtw.wide %>% 
  pivot_longer(-x, names_to = 'GeneID', values_to = 'expr')

sc.atac.mu.scale <- sc.atac.dtw.wide %>% 
  pivot_longer(-x, names_to = 'GeneID', values_to = 'expr')


sc.rna.peak.order <- sc.rna.mu.scale %>% group_by(GeneID) %>% summarise(peak.ord = getCurvePeakLoc(x, expr))
sc.rna.mu.scale <- left_join(sc.rna.mu.scale, sc.rna.peak.order, by = 'GeneID')

sc.rna.mu.scale$GeneID <- factor(sc.rna.mu.scale$GeneID, 
                                 levels = unique(sc.rna.mu.scale$GeneID[order(-sc.rna.mu.scale$peak.ord)]))


sc.atac.peak.order <- sc.atac.mu.scale %>% group_by(GeneID) %>% summarise(peak.ord = getCurvePeakLoc(x, expr))
sc.atac.mu.scale <- left_join(sc.atac.mu.scale, sc.atac.peak.order, by = 'GeneID')

## Order by RNA peak expression
sc.atac.mu.scale$GeneID <- factor(sc.atac.mu.scale$GeneID, 
                                 levels = unique(sc.atac.mu.scale$GeneID[order(-sc.rna.mu.scale$peak.ord)]))

sc.atac.mu.scale$GeneID <- factor(sc.atac.mu.scale$GeneID,
                                  levels = unique(sc.atac.mu.scale$GeneID[order(-sc.atac.mu.scale$peak.ord)]))




p1 <- ggplot(sc.rna.mu.scale, aes(x = x, y = GeneID, fill = expr)) + 
  geom_tile() + 
  #scale_x_discrete(expand=c(0,0)) +
  ylab("Genes") + xlab("ps time") + ggtitle("scRNA") + 
  #scale_fill_gradientn(colours = hm.palette(10)) +
  scale_fill_gradientn(colours = viridis::inferno(10)) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks = element_blank(),
    axis.text.y  = element_blank(),
    plot.title = element_text(size=14, face="bold"),
    axis.title.x = element_text(size=14, face="bold"),
    axis.title.y = element_text(size=14, face="bold"),
    legend.position = "none")


plot(p1)  


ggsave(filename="../Output/toxo_cdc/figures/scRNA_intra_heatmap_no_clusters.png",
       plot=p1,
       width = 6, height = 7,
       units = "in", # other options are "in", "cm", "mm"
       dpi = 300
)


p2 <-  ggplot(sc.atac.mu.scale, aes(x = x, y = GeneID, fill = expr)) + 
  geom_tile() + 
  #scale_x_discrete(expand=c(0,0)) +
  ylab("Genes") + xlab("ps time") + ggtitle("scATAC") + 
  #scale_fill_gradientn(colours = hm.palette(10)) +
  scale_fill_gradientn(colours = viridis::inferno(10)) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks = element_blank(),
    axis.text.y  = element_blank(),
    plot.title = element_text(size=14, face="bold"),
    axis.title.x = element_text(size=14, face="bold"),
    axis.title.y = element_text(size=14, face="bold"),
    legend.position = "none")


plot(p2)  

ggsave(filename="../Output/toxo_cdc/figures/scATAC_heatmap_no_clusters.png",
       plot=p2,
       width = 6, height = 7,
       units = "in", # other options are "in", "cm", "mm"
       dpi = 300
)

tans.points <- function(mu.scale, lam = 0.05, prob = 0.01){
  peak.locs <- mu.scale %>% dplyr::select(GeneID, peak.ord) %>% distinct()
  spline.fit.peaks <- smooth.spline(x = 1:nrow(peak.locs), y = sort(peak.locs$peak.ord), lambda = lam)
  
  s.0 <- predict(spline.fit.peaks, seq(1, nrow(peak.locs), by = 0.1), deriv=0)
  s.1 <- predict(spline.fit.peaks, seq(1, nrow(peak.locs), by = 0.1), deriv=1)
  s.2 <- predict(spline.fit.peaks, seq(1, nrow(peak.locs), by = 0.1), deriv=2)
  
  spline.fit.peaks.smooth <- data.frame(x = seq(1, nrow(peak.locs), by = 0.1), 
                                        s0=s.0$y, s1=s.1$y, s2 = s.2$y)
  
  ## Get locations of peaks and vallies
  max.loc <- getCurvePeakLoc(spline.fit.peaks.smooth$x, spline.fit.peaks.smooth$s2, prob = prob)
  min.loc <- getCurvePeakLoc(spline.fit.peaks.smooth$x, -spline.fit.peaks.smooth$s2, prob = prob)
  
  transition.points <- predict(spline.fit.peaks, sort(c(max.loc, min.loc)))
  
  L <- list(peak.locs = peak.locs, spline.fit.peaks.smooth = spline.fit.peaks.smooth, 
            transition.points = transition.points)
  
  return(L)
  
}

L.trans.rna <- tans.points(sc.rna.mu.scale, lam = 1, prob = 0.0000)

saveRDS(L.trans.rna$transition.points, '../Input/toxo_cdc/rds/transition_points_rna.rds')


pdf(file="../Output/toxo_cdc/figures/peak_time_phase_rna.pdf",
    width=8, height=8)
par(mfrow = c(2,1))
plot(L.trans.rna$spline.fit.peaks.smooth$x, L.trans.rna$spline.fit.peaks.smooth$s0, 
     type = 'l', col = 'red', xlab = "genes", ylab = 'peak time', lwd = 2, main = 'sorted peak time')
points(sort(L.trans.rna$peak.locs$peak.ord), pch = 20,  cex = 0.1, col = 'blue')
abline(v = L.trans.rna$transition.points$x, col = "orange", lty=2, lwd = 2)

plot(L.trans.rna$spline.fit.peaks.smooth$x, L.trans.rna$spline.fit.peaks.smooth$s2, 
     type = 'l', col = 'green', xlab = "genes", ylab = 'peak time', lwd = 2, main = 'second derivative')
abline(v = L.trans.rna$transition.points$x, col = "orange", lty=2, lwd = 2)
dev.off()



L.trans.atac <- tans.points(sc.atac.mu.scale, lam = 1, prob = 0.0)

saveRDS(L.trans.atac$transition.points, '../Input/toxo_cdc/rds/transition_points_atac.rds')


pdf(file="../Output/toxo_cdc/figures/peak_time_phase_atac.pdf",
    width=8, height=8)
par(mfrow = c(2,1))
plot(L.trans.atac$spline.fit.peaks.smooth$x, L.trans.atac$spline.fit.peaks.smooth$s0, 
     type = 'l', col = 'red', xlab = "genes", ylab = 'peak time', lwd = 2, main = 'sorted peak time')
points(sort(L.trans.atac$peak.locs$peak.ord), pch = 20,  cex = 0.1, col = 'blue')
abline(v = L.trans.atac$transition.points$x, col = "orange", lty=2, lwd = 2)

plot(L.trans.atac$spline.fit.peaks.smooth$x, L.trans.atac$spline.fit.peaks.smooth$s2, 
     type = 'l', col = 'green', xlab = "genes", ylab = 'peak time', lwd = 2, main = 'second derivative')
abline(v = L.trans.atac$transition.points$x, col = "orange", lty=2, lwd = 2)
dev.off()



sc.rna.mu.scale <- sc.rna.mu.scale %>% 
  mutate(cluster.rna = case_when(peak.ord <  L.trans.rna$transition.points$y[1] ~ "RC1",
                             peak.ord >= L.trans.rna$transition.points$y[1] & 
                               peak.ord < L.trans.rna$transition.points$y[2] ~ "RC2",
                             peak.ord >= L.trans.rna$transition.points$y[2] ~ "RC3"),
         cluster.atac = case_when(peak.ord <  L.trans.atac$transition.points$y[1] ~ "AC1",
                                 peak.ord >= L.trans.atac$transition.points$y[1] & 
                                   peak.ord < L.trans.atac$transition.points$y[2] ~ "AC2",
                                 peak.ord >= L.trans.atac$transition.points$y[2] ~ "AC3"))

sc.rna.mu.scale$cluster.rna <- factor(sc.rna.mu.scale$cluster.rna, levels = sort(unique(sc.rna.mu.scale$cluster.rna)))
sc.rna.mu.scale$cluster.atac <- factor(sc.rna.mu.scale$cluster.atac, levels = sort(unique(sc.rna.mu.scale$cluster.atac)))

saveRDS(sc.rna.mu.scale, '../Input/toxo_cdc/rds/sc_rna_spline_mu_scale_transition_clusters.rds')
write.xlsx(sc.rna.mu.scale, '../Output/toxo_cdc/tabels/')

p1 <- ggplot(sc.rna.mu.scale, aes(x = x, y = GeneID, fill = expr)) + 
  geom_tile() + 
  facet_grid(cluster.rna~., scales = "free",  space='free',
             labeller=label_wrap_gen(multi_line = TRUE))+
  ylab("Genes") + xlab("time/cells") + ggtitle('scRNA') + 
  #scale_fill_gradientn(colours = hm.palette(10)) +
  scale_fill_gradientn(colours = viridis::inferno(10)) +
  theme(panel.spacing = unit(0.1, "lines")) + 
  theme(
    axis.text.x = element_blank(),
    axis.ticks = element_blank(),
    axis.text.y  = element_blank(),
    plot.title = element_text(size=14, face="bold"),
    axis.title.x = element_text(size=14, face="bold"),
    axis.title.y = element_text(size=14, face="bold"),
    legend.position = "none")#+
  #theme(plot.margin=unit(c(2.5,2.5,2.5,2.5),"cm"))

plot(p1)

ggsave(filename="../Output/toxo_cdc/figures/scRNA_intra_heatmap_rna_clust.png",
       plot=p1,
       width = 6, height = 7,
       units = "in", # other options are "in", "cm", "mm"
       dpi = 300
)


p2 <- ggplot(sc.rna.mu.scale, aes(x = x, y = GeneID, fill = expr)) + 
  geom_tile() + 
  facet_grid(cluster.atac~., scales = "free",  space='free',
             labeller=label_wrap_gen(multi_line = TRUE))+
  ylab("Genes") + xlab("time/cells") + ggtitle('scRNA') + 
  #scale_fill_gradientn(colours = hm.palette(10)) +
  scale_fill_gradientn(colours = viridis::inferno(10)) +
  theme(panel.spacing = unit(0.1, "lines")) + 
  theme(
    axis.text.x = element_blank(),
    axis.ticks = element_blank(),
    axis.text.y  = element_blank(),
    plot.title = element_text(size=14, face="bold"),
    axis.title.x = element_text(size=14, face="bold"),
    axis.title.y = element_text(size=14, face="bold"),
    legend.position = "none")#+
#theme(plot.margin=unit(c(2.5,2.5,2.5,2.5),"cm"))

plot(p2)

 
ggsave(filename="../Output/toxo_cdc/figures/scRNA_intra_heatmap_atac_clust.png",
         plot=p2,
         width = 6, height = 7,
         units = "in", # other options are "in", "cm", "mm"
         dpi = 300
)
  




sc.atac.mu.scale <- sc.atac.mu.scale %>% 
  mutate(cluster.rna = case_when(peak.ord <  L.trans.rna$transition.points$y[1] ~ "RC1",
                                 peak.ord >= L.trans.rna$transition.points$y[1] & 
                                   peak.ord < L.trans.rna$transition.points$y[2] ~ "RC2",
                                 peak.ord >= L.trans.rna$transition.points$y[2] & 
                                   peak.ord < L.trans.rna$transition.points$y[3] ~ "RC3",
                                 peak.ord >= L.trans.rna$transition.points$y[3] & 
                                   peak.ord < L.trans.rna$transition.points$y[4] ~ "RC4",
                                 peak.ord >= L.trans.rna$transition.points$y[4] & 
                                   peak.ord < L.trans.rna$transition.points$y[5] ~ "RC5",
                                 peak.ord >= L.trans.rna$transition.points$y[5] & 
                                   peak.ord < L.trans.rna$transition.points$y[6] ~ "RC6",
                                 peak.ord >= L.trans.rna$transition.points$y[6] & 
                                   peak.ord < L.trans.rna$transition.points$y[7] ~ "RC7",
                                 peak.ord >= L.trans.rna$transition.points$y[7] & 
                                   peak.ord < L.trans.rna$transition.points$y[8] ~ "RC8",
                                 peak.ord >= L.trans.rna$transition.points$y[8] ~ "RC9"),
         cluster.atac = case_when(peak.ord <  L.trans.atac$transition.points$y[1] ~ "AC1",
                                  peak.ord >= L.trans.atac$transition.points$y[1] & 
                                    peak.ord < L.trans.atac$transition.points$y[2] ~ "AC2",
                                  peak.ord >= L.trans.atac$transition.points$y[2] & 
                                    peak.ord < L.trans.atac$transition.points$y[3] ~ "AC3",
                                  peak.ord >= L.trans.atac$transition.points$y[3] & 
                                    peak.ord < L.trans.atac$transition.points$y[4] ~ "AC4",
                                  peak.ord >= L.trans.atac$transition.points$y[4] ~ "AC5"))

sc.atac.mu.scale$cluster.rna <- factor(sc.atac.mu.scale$cluster.rna, levels = sort(unique(sc.atac.mu.scale$cluster.rna)))
sc.atac.mu.scale$cluster.atac <- factor(sc.atac.mu.scale$cluster.atac, levels = sort(unique(sc.atac.mu.scale$cluster.atac)))

saveRDS(sc.atac.mu.scale, '../Input/toxo_cdc/rds/sc_atac_spline_mu_scale_transition_clusters.rds')


p1 <- ggplot(sc.atac.mu.scale, aes(x = x, y = GeneID, fill = expr)) + 
  geom_tile() + 
  facet_grid(cluster.rna~., scales = "free",  space='free',
             labeller=label_wrap_gen(multi_line = TRUE))+
  ylab("Genes") + xlab("time/cells") + ggtitle('sATAC') + 
  #scale_fill_gradientn(colours = hm.palette(10)) +
  scale_fill_gradientn(colours = viridis::inferno(10)) +
  theme(panel.spacing = unit(0.1, "lines")) + 
  theme(
    axis.text.x = element_blank(),
    axis.ticks = element_blank(),
    axis.text.y  = element_blank(),
    plot.title = element_text(size=14, face="bold"),
    axis.title.x = element_text(size=14, face="bold"),
    axis.title.y = element_text(size=14, face="bold"),
    legend.position = "none")#+
#theme(plot.margin=unit(c(2.5,2.5,2.5,2.5),"cm"))

plot(p1)


ggsave(filename="../Output/toxo_cdc/figures/scRNA_atac_heatmap_rna_clust.png",
       plot=p1,
       width = 6, height = 7,
       units = "in", # other options are "in", "cm", "mm"
       dpi = 300
)


p2 <- ggplot(sc.atac.mu.scale, aes(x = x, y = GeneID, fill = expr)) + 
  geom_tile() + 
  facet_grid(cluster.atac~., scales = "free",  space='free',
             labeller=label_wrap_gen(multi_line = TRUE))+
  ylab("Genes") + xlab("time/cells") + ggtitle('sATAC') + 
  #scale_fill_gradientn(colours = hm.palette(10)) +
  scale_fill_gradientn(colours = viridis::inferno(10)) +
  theme(panel.spacing = unit(0.1, "lines")) + 
  theme(
    axis.text.x = element_blank(),
    axis.ticks = element_blank(),
    axis.text.y  = element_blank(),
    plot.title = element_text(size=14, face="bold"),
    axis.title.x = element_text(size=14, face="bold"),
    axis.title.y = element_text(size=14, face="bold"),
    legend.position = "none")#+
#theme(plot.margin=unit(c(2.5,2.5,2.5,2.5),"cm"))

plot(p2)


ggsave(filename="../Output/toxo_cdc/figures/scRNA_atac_heatmap_atac_clust.png",
       plot=p2,
       width = 6, height = 7,
       units = "in", # other options are "in", "cm", "mm"
       dpi = 300
)


## Map back to PCA


sc.rna.cell.cycle.genes.df.adj <- sc.rna.cell.cycle.genes.df.adj %>% 
  mutate(cluster.rna = case_when(adj.time <  L.trans.rna$transition.points$y[1] ~ "RC1",
                                 adj.time >= L.trans.rna$transition.points$y[1] & 
                                   adj.time < L.trans.rna$transition.points$y[2] ~ "RC2",
                                 adj.time >= L.trans.rna$transition.points$y[2] & 
                                   adj.time < L.trans.rna$transition.points$y[3] ~ "RC3",
                                 adj.time >= L.trans.rna$transition.points$y[3] & 
                                   adj.time < L.trans.rna$transition.points$y[4] ~ "RC4",
                                 adj.time >= L.trans.rna$transition.points$y[4] & 
                                   adj.time < L.trans.rna$transition.points$y[5] ~ "RC5",
                                 adj.time >= L.trans.rna$transition.points$y[5] & 
                                   adj.time < L.trans.rna$transition.points$y[6] ~ "RC6",
                                 adj.time >= L.trans.rna$transition.points$y[6] & 
                                   adj.time < L.trans.rna$transition.points$y[7] ~ "RC7",
                                 adj.time >= L.trans.rna$transition.points$y[7] & 
                                   adj.time < L.trans.rna$transition.points$y[8] ~ "RC8",
                                 adj.time >= L.trans.rna$transition.points$y[8] ~ "RC9"),
         cluster.atac = case_when(adj.time <  L.trans.atac$transition.points$y[1] ~ "AC1",
                                  adj.time >= L.trans.atac$transition.points$y[1] & 
                                    adj.time < L.trans.atac$transition.points$y[2] ~ "AC2",
                                  adj.time >= L.trans.atac$transition.points$y[2] & 
                                    adj.time < L.trans.atac$transition.points$y[3] ~ "AC3",
                                  adj.time >= L.trans.atac$transition.points$y[3] & 
                                    adj.time < L.trans.atac$transition.points$y[4] ~ "AC4",
                                  adj.time >= L.trans.atac$transition.points$y[4] ~ "AC5"))

sc.rna.cell.cycle.genes.df.adj$cluster.rna <- factor(sc.rna.cell.cycle.genes.df.adj$cluster.rna, 
                                                     levels = sort(unique(sc.rna.cell.cycle.genes.df.adj$cluster.rna)))
sc.rna.cell.cycle.genes.df.adj$cluster.atac <- factor(sc.rna.cell.cycle.genes.df.adj$cluster.atac, 
                                                      levels = sort(unique(sc.rna.cell.cycle.genes.df.adj$cluster.atac)))



sc.atac.cell.cycle.genes.df.adj <- sc.atac.cell.cycle.genes.df.adj %>% 
  mutate(cluster.rna = case_when(adj.time <  L.trans.rna$transition.points$y[1] ~ "RC1",
                                 adj.time >= L.trans.rna$transition.points$y[1] & 
                                   adj.time < L.trans.rna$transition.points$y[2] ~ "RC2",
                                 adj.time >= L.trans.rna$transition.points$y[2] & 
                                   adj.time < L.trans.rna$transition.points$y[3] ~ "RC3",
                                 adj.time >= L.trans.rna$transition.points$y[3] & 
                                   adj.time < L.trans.rna$transition.points$y[4] ~ "RC4",
                                 adj.time >= L.trans.rna$transition.points$y[4] & 
                                   adj.time < L.trans.rna$transition.points$y[5] ~ "RC5",
                                 adj.time >= L.trans.rna$transition.points$y[5] & 
                                   adj.time < L.trans.rna$transition.points$y[6] ~ "RC6",
                                 adj.time >= L.trans.rna$transition.points$y[6] & 
                                   adj.time < L.trans.rna$transition.points$y[7] ~ "RC7",
                                 adj.time >= L.trans.rna$transition.points$y[7] & 
                                   adj.time < L.trans.rna$transition.points$y[8] ~ "RC8",
                                 adj.time >= L.trans.rna$transition.points$y[8] ~ "RC9"),
         cluster.atac = case_when(adj.time <  L.trans.atac$transition.points$y[1] ~ "AC1",
                                  adj.time >= L.trans.atac$transition.points$y[1] & 
                                    adj.time < L.trans.atac$transition.points$y[2] ~ "AC2",
                                  adj.time >= L.trans.atac$transition.points$y[2] & 
                                    adj.time < L.trans.atac$transition.points$y[3] ~ "AC3",
                                  adj.time >= L.trans.atac$transition.points$y[3] & 
                                    adj.time < L.trans.atac$transition.points$y[4] ~ "AC4",
                                  adj.time >= L.trans.atac$transition.points$y[4] ~ "AC5"))

sc.atac.cell.cycle.genes.df.adj$cluster.rna <- factor(sc.atac.cell.cycle.genes.df.adj$cluster.rna, 
                                                     levels = sort(unique(sc.atac.cell.cycle.genes.df.adj$cluster.rna)))
sc.atac.cell.cycle.genes.df.adj$cluster.atac <- factor(sc.atac.cell.cycle.genes.df.adj$cluster.atac, 
                                                      levels = sort(unique(sc.atac.cell.cycle.genes.df.adj$cluster.atac)))



saveRDS(sc.rna.cell.cycle.genes.df.adj, '../Input/toxo_cdc/rds/rna_seq_gene_expr_pseudo_time_transition_clusters.rds')
saveRDS(sc.atac.cell.cycle.genes.df.adj, '../Input/toxo_cdc/rds/atac_seq_gene_access_pseudo_time_transition_clusters.rds')


S.O.integrated <- readRDS('../Input/toxo_cdc/rds/S.O.intra_atac_integrated.rds')
Idents(S.O.integrated) <- 'orig.ident'

atac_sub<- subset(S.O.integrated, ident = 'scATAC')
rna_sub <- subset(S.O.integrated, ident = 'scRNA')


rna.trans.d.meta <- sc.rna.cell.cycle.genes.df.adj %>% 
  dplyr::select(Sample, adj.time, adj.time.idx, cluster.rna, cluster.atac) %>% distinct() %>%
  data.frame()
rownames(rna.trans.d.meta) <- rna.trans.d.meta$Sample

atac.trans.d.meta <- sc.atac.cell.cycle.genes.df.adj %>% 
  dplyr::select(Sample, adj.time, adj.time.idx, cluster.rna, cluster.atac) %>% distinct() %>%
  data.frame()
rownames(atac.trans.d.meta) <- atac.trans.d.meta$Sample

rna_sub <- AddMetaData(rna_sub, rna.trans.d.meta)
atac_sub <- AddMetaData(atac_sub, atac.trans.d.meta)

saveRDS(rna_sub, '../Input/toxo_cdc/rds/S.O.intra_rna_from_integrated_trnasition_clusters.rds')
saveRDS(atac_sub, '../Input/toxo_cdc/rds/S.O.atac_from_integrated_trnasition_clusters.rds')




## Scale the pt based on know biology: Radke et. al 2000
scalePT <- function(pt.boundaries){
  G <- c(0, 3) # 3h
  S <- c(3, 4.7) # 1.7h
  M <- c(4.7, 5) # ~20 min
  C <- c(5, 6) # 1h
  
  G.adj.pt <- c(min(pt.phase$pt[pt.phase$phase %in% c('G1.a', 'G1.b')]),
                max(pt.phase$pt[pt.phase$phase %in% c('G1.a', 'G1.b')]))
  
}

## Proportions
tmp <- rna_sub@meta.data

stats <- tmp %>% group_by(phase) %>% mutate(total.cells = n()) %>% 
  summarise(perc = total.cells[1]/nrow(tmp)) 


p <- ggplot(data=stats, aes(x=phase, y=perc, fill = phase)) +
  geom_bar(stat="identity", position=position_dodge(width = 1.0), width=0.8) +
  geom_text(aes(label=round(perc, 2)), vjust=2,  color="black", 
            size=6, fontface="bold", position=position_dodge(1.0), angle=0)+
  #scale_fill_manual(values = c("G1" = "#ff6c67","S/M" ='#a2a700', 'C' = '#00c377')) +
  theme_minimal() + 
  theme(panel.spacing = unit(0.5, "lines")) + 
  theme(axis.text.x = element_text(face="bold", size=16, angle=0)) +
  theme(axis.text.y = element_text(face="bold", size=16, angle=0)) +
  theme(
    axis.title.x = element_text(size=18, face="bold"),
    axis.title.y = element_text(size=18, face="bold")
  ) #+  coord_flip()

plot(p)
