library(plotly)
library(Seurat)
library(tidyverse)
library(SeuratWrappers)
library(openxlsx)
library(slingshot)


source('./util_funcs.R')

getPcaUmapMetaData <- function(S.O){
  pc <- S.O@reductions$pca@cell.embeddings
  pc <- data.frame(pc) %>% dplyr::mutate(Sample = rownames(pc)) %>% 
    transmute(Sample = Sample, PC_1 = PC_1, PC_2 = PC_2, PC_3 = PC_3)
  umap <- S.O@reductions$umap@cell.embeddings
  umap <- data.frame(umap) %>% dplyr::mutate(Sample = rownames(umap)) %>% 
    transmute(Sample = Sample, UMAP_1 = UMAP_1, UMAP_2 = UMAP_2, UMAP_3 = UMAP_3)
  
  meta.data <- data.frame(Sample = rownames(S.O@meta.data), phase = S.O@meta.data$phase, 
                          spp = S.O@meta.data$spp, clusters = S.O@meta.data$seurat_clusters)
  meta.data <- left_join(meta.data,
                         pc, by = 'Sample')
  meta.data <- left_join(meta.data, umap, by = 'Sample')
  return(meta.data)  
}


S.O.intra <- readRDS('../Input/toxo_cdc/rds/S.O.intra_lables.rds')
S.O.intra <- RunUMAP(S.O.intra, dims = 1:13, n.components = 3)

pca.tg <- getPcaUmapMetaData(S.O.intra)
pca.tg$phase <- factor(pca.tg$phase, 
                       levels = c('G1.a', 'G1.b', 'S', 'M', 'C'))




# sds <- slingshot(Embeddings(S.O.intra, 'pca'), 
#                  clusterLabels = S.O.intra$phase, 
#                  start.clus = 'G1.a', end.clus = 'C', stretch = 2)
# 
# pt <- slingPseudotime(sds)
# 
# s.data <- list(pt = pt, curves = sds@metadata$curves)
# trace1 <- data.frame(s.data$curves$Lineage1$s)
# trace1$type = 'curve'
# trace1$color = 'trace1'
# 
# trace2 <- data.frame(s.data$curves$Lineage2$s)
# trace2$type = 'curve'
# trace2$color = 'trace2'
# 
# traces <- bind_rows(trace1, trace2) %>% dplyr::select('PC_1', 'PC_2', 'PC_3', 'type', 'color')




x <- pca.tg %>% select(contains('PC')) %>% as.matrix()
rownames(x) <- pca.tg$Sample
#fit <- principal_curve(x, start = as.matrix(coord))
fit <- principal_curve(x, smoother = 'periodic_lowess')
pt <- fit$lambda

## reversing the order of time
pt <- max(pt) - pt
pt <- (pt - min(pt))/(max(pt) - min(pt))
cell.ord <- fit$ord[seq(length(fit$ord), 1, by = -1)]


sds.data <- data.frame(Sample = rownames(fit$s), 
                     cell.ord = cell.ord,
                     pt = pt,
                     PC_1 = fit$s[,1], 
                     PC_2 = fit$s[,2],
                     PC_3 = fit$s[,3])


fig <- plot_ly(pca.tg, x = ~PC_1, y = ~PC_2, z = ~PC_3, color = ~phase, 
               colors = c('#f8766d', '#a2a400', '#00bf7d', '#04b0f6', '#e76bf3')) 
fig <- fig %>% add_markers(size = 2)
fig <- fig %>% layout(scene = list(xaxis = list(title = 'PC_1'),
                                   yaxis = list(title = 'PC_2'),
                                   zaxis = list(title = 'PC_3')))

fig


PCAs <- pca.tg %>% select(contains('PC'))
PCAs$type <- 'cells'
PCAs$color <- pca.tg$phase

traces <- sds.data %>% dplyr::select(contains('PC'))
traces$type <- 'curve'
traces$color <- 'trace'

df <- bind_rows(PCAs, traces)
df$color <- factor(df$color, levels = c('G1.a', 'G1.b', 'S', 'M', 'C', 'trace'))
df$type <- factor(df$type, c('cells', 'curve'))

fig <- plot_ly(df, x = ~PC_1, y = ~PC_2, z = ~PC_3, color = ~color, 
               colors = c('#f8766d', '#a2a400', '#00bf7d', '#04b0f6', '#e76bf3', 'black')) 
fig <- fig %>% add_markers(size = c(4,4)[df$type])
fig <- fig %>% layout(scene = list(xaxis = list(title = 'PC_1'),
                                   yaxis = list(title = 'PC_2'),
                                   zaxis = list(title = 'PC_3')))

fig

plot(sds.data$PC_1[sds.data$cell.ord], sds.data$PC_2[sds.data$cell.ord], type = 'l')
