library(Seurat)
library(openxlsx)
library(ggplot2)
library(ggfortify)
library(jcolors)
require(gridExtra)
library(grid)
library(matrixStats)
library(tidyverse)
library(RColorBrewer)
library(parallel)
library(openxlsx)
library(plotly)


source('./util_funcs.R')

## Extracting the data for control over PCA directions
getPcaMetaData <- function(S.O){
  pc <- S.O@reductions$pca@cell.embeddings
  pc <- data.frame(pc) %>% dplyr::mutate(Sample = rownames(pc)) %>% 
    transmute(Sample = Sample, PC_1 = PC_1, PC_2 = PC_2, PC_3 = PC_3)
  umap <- S.O@reductions$umap@cell.embeddings
  umap <- data.frame(umap) %>% dplyr::mutate(Sample = rownames(umap)) %>% 
    transmute(Sample = Sample, UMAP_1 = UMAP_1, UMAP_2 = UMAP_2)
  
  meta.data <- data.frame(Sample = rownames(S.O@meta.data), 
                          spp = S.O@meta.data$spp, phase = S.O@meta.data$phase, transition = S.O@meta.data$transition.atac)
  meta.data <- left_join(meta.data,
                         pc, by = 'Sample')
  meta.data <- left_join(meta.data, umap, by = 'Sample')
  return(meta.data)  
}



# rna_sub <- readRDS('../Input/toxo_cdc/rds/S.O_intra_lables_pt.rds')
# atac_sub <- readRDS('../Input/toxo_cdc/rds/S.O_intra_atac_lables_pt.rds')


## With Transition maps
rna_sub  <- readRDS('../Input/toxo_cdc/rds/S.O.intra_rna_atac_trnasition.rds')
atac_sub <- readRDS('../Input/toxo_cdc/rds/S.O.intra_atac_atac_trnasition.rds')

rna_sub.pca <- getPcaMetaData(rna_sub)


p1  <- ggplot(rna_sub.pca, aes(x= PC_1,y=PC_2)) +
  geom_point(aes(#fill = lable.prob,
    fill = phase,
    color = phase, 
  ), #color = 'blue', 
  alpha = 0.4,
  shape=21, size = 1)+ 
  # scale_color_manual(values = c("BBig" = "firebrick","BBov" ="darkorchid3", 'BDiv_Cow' = 'darkslateblue', 
  #                               'BDiv_Human' = 'darkolivegreen4')) +
  # scale_fill_manual(values = c("BBig" = "firebrick","BBov" ="darkorchid3", 'BDiv_Cow' = 'darkslateblue', 
  #                              'BDiv_Human' = 'darkolivegreen4')) +
  
  theme_bw(base_size = 14) +
  theme(legend.position = "right") +
  #scale_fill_gradientn(colours = viridis::inferno(10)) +
  #scale_fill_gradientn(colours = col_range(10)) +
  #scale_fill_gradient(low = "gray66", high = "blue", midpoint = mid_point) + 
  #scale_fill_brewer(palette = "BuPu") +
  ylab('PC1') + xlab('PC2') +
  theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 12, face="bold")) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 12, face="bold")) +
  theme(strip.background = element_rect(colour="black", fill="white",
                                        size=0.5, linetype="solid")) +
  #facet_wrap(spp~.) + 
  #ggtitle(titles[i]) +
  theme(
    plot.title = element_text(size=14, face = "bold.italic", color = 'red'),
    axis.title.x = element_text(size=14, face="bold", hjust = 1),
    axis.title.y = element_text(size=14, face="bold")
  ) + 
  theme(#legend.position = c(0.15, 0.85),
    legend.position = 'none',
    legend.title = element_text(colour="black", size=12, 
                                face="bold"),
    legend.text = element_text(colour="black", size=12, 
                               face="bold"))


plot(p1)

ggsave(filename="../Output/toxo_cdc/figures/pca_rna.pdf",
       plot=p1,
       width = 6, height = 6,
       units = "in"
)



atac_sub.pca <- getPcaMetaData(atac_sub)


p2  <- ggplot(atac_sub.pca, aes(x= PC_1,y=PC_2)) +
  geom_point(aes(#fill = lable.prob,
    fill = phase,
    color = phase, 
  ), #color = 'blue', 
  alpha = 0.4,
  shape=21, size = 1)+ 
  # scale_color_manual(values = c("BBig" = "firebrick","BBov" ="darkorchid3", 'BDiv_Cow' = 'darkslateblue', 
  #                               'BDiv_Human' = 'darkolivegreen4')) +
  # scale_fill_manual(values = c("BBig" = "firebrick","BBov" ="darkorchid3", 'BDiv_Cow' = 'darkslateblue', 
  #                              'BDiv_Human' = 'darkolivegreen4')) +
  
  theme_bw(base_size = 14) +
  theme(legend.position = "right") +
  #scale_fill_gradientn(colours = viridis::inferno(10)) +
  #scale_fill_gradientn(colours = col_range(10)) +
  #scale_fill_gradient(low = "gray66", high = "blue", midpoint = mid_point) + 
  #scale_fill_brewer(palette = "BuPu") +
  ylab('') + xlab('PC2') +
  theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 12, face="bold")) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 12, face="bold")) +
  theme(strip.background = element_rect(colour="black", fill="white",
                                        size=0.5, linetype="solid")) +
  #facet_wrap(spp~.) + 
  #ggtitle(titles[i]) +
  theme(
    plot.title = element_text(size=14, face = "bold.italic", color = 'red'),
    axis.title.x = element_text(size=14, face="bold", hjust = 1),
    axis.title.y = element_text(size=14, face="bold")
  ) + 
  theme(#legend.position = c(0.15, 0.85),
    legend.position = 'none',
    legend.title = element_text(colour="black", size=12, 
                                face="bold"),
    legend.text = element_text(colour="black", size=12, 
                               face="bold"))


plot(p2)

ggsave(filename="../Output/toxo_cdc/figures/pca_atac.pdf",
       plot=p2,
       width = 6, height = 6,
       units = "in"
)



p3  <- ggplot(rna_sub.pca, aes(x= PC_1,y=PC_2)) +
  geom_point(aes(#fill = lable.prob,
    fill = transition,
    color = transition, 
  ), #color = 'blue', 
  alpha = 0.4,
  shape=21, size = 1)+ 
  # scale_color_manual(values = c("BBig" = "firebrick","BBov" ="darkorchid3", 'BDiv_Cow' = 'darkslateblue', 
  #                               'BDiv_Human' = 'darkolivegreen4')) +
  # scale_fill_manual(values = c("BBig" = "firebrick","BBov" ="darkorchid3", 'BDiv_Cow' = 'darkslateblue', 
  #                              'BDiv_Human' = 'darkolivegreen4')) +
  
  theme_bw(base_size = 14) +
  theme(legend.position = "right") +
  #scale_fill_gradientn(colours = viridis::inferno(10)) +
  #scale_fill_gradientn(colours = col_range(10)) +
  #scale_fill_gradient(low = "gray66", high = "blue", midpoint = mid_point) + 
  #scale_fill_brewer(palette = "BuPu") +
  ylab('PC1') + xlab('PC2') +
  theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 12, face="bold")) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 12, face="bold")) +
  theme(strip.background = element_rect(colour="black", fill="white",
                                        size=0.5, linetype="solid")) +
  #facet_wrap(spp~.) + 
  #ggtitle(titles[i]) +
  theme(
    plot.title = element_text(size=14, face = "bold.italic", color = 'red'),
    axis.title.x = element_text(size=14, face="bold", hjust = 1),
    axis.title.y = element_text(size=14, face="bold")
  ) + 
  theme(#legend.position = c(0.15, 0.85),
    legend.position = 'none',
    legend.title = element_text(colour="black", size=12, 
                                face="bold"),
    legend.text = element_text(colour="black", size=12, 
                               face="bold"))


plot(p3)

ggsave(filename="../Output/toxo_cdc/figures/pca_rna_transition.pdf",
       plot=p3,
       width = 6, height = 6,
       units = "in"
)




p4  <- ggplot(atac_sub.pca, aes(x= PC_1,y=PC_2)) +
  geom_point(aes(#fill = lable.prob,
    fill = transition,
    color = transition, 
  ), #color = 'blue', 
  alpha = 0.4,
  shape=21, size = 1)+ 
  # scale_color_manual(values = c("BBig" = "firebrick","BBov" ="darkorchid3", 'BDiv_Cow' = 'darkslateblue', 
  #                               'BDiv_Human' = 'darkolivegreen4')) +
  # scale_fill_manual(values = c("BBig" = "firebrick","BBov" ="darkorchid3", 'BDiv_Cow' = 'darkslateblue', 
  #                              'BDiv_Human' = 'darkolivegreen4')) +
  
  theme_bw(base_size = 14) +
  theme(legend.position = "right") +
  #scale_fill_gradientn(colours = viridis::inferno(10)) +
  #scale_fill_gradientn(colours = col_range(10)) +
  #scale_fill_gradient(low = "gray66", high = "blue", midpoint = mid_point) + 
  #scale_fill_brewer(palette = "BuPu") +
  ylab('') + xlab('PC2') +
  theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 12, face="bold")) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 12, face="bold")) +
  theme(strip.background = element_rect(colour="black", fill="white",
                                        size=0.5, linetype="solid")) +
  #facet_wrap(spp~.) + 
  #ggtitle(titles[i]) +
  theme(
    plot.title = element_text(size=14, face = "bold.italic", color = 'red'),
    axis.title.x = element_text(size=14, face="bold", hjust = 1),
    axis.title.y = element_text(size=14, face="bold")
  ) + 
  theme(#legend.position = c(0.15, 0.85),
    legend.position = 'none',
    legend.title = element_text(colour="black", size=12, 
                                face="bold"),
    legend.text = element_text(colour="black", size=12, 
                               face="bold"))


plot(p4)

ggsave(filename="../Output/toxo_cdc/figures/pca_atac_transiiton.pdf",
       plot=p4,
       width = 6, height = 6,
       units = "in"
)

## Cell cycle Markers
Intra.markers.sig <- readRDS('../Input/toxo_cdc/rds/Intra_markers_sig.rds')

ss <- Intra.markers.sig %>% group_by(cluster) %>% summarise(num.DEG = n())
ss$cluster <- factor(ss$cluster, levels = c('G1.a', 'G1.b', 'S', 'M', 'C'))
p <- ggplot(data=ss, aes(x=cluster, y=num.DEG, fill=cluster)) +
  geom_bar(stat="identity")+
  geom_text(aes(label=num.DEG), vjust=1.6, color="black", size=5, fontface="bold")+
  theme_minimal() + 
  ylab('') + xlab('') +
  theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 16, face="bold")) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 16, face="bold")) +
  theme(strip.background = element_rect(colour="black", fill="white",
                                        size=0.5, linetype="solid")) +
  #facet_wrap(spp~.) + 
  #ggtitle(titles[i]) +
  theme(
    plot.title = element_text(size=18, face = "bold.italic", color = 'red'),
    axis.title.x = element_text(size=18, face="bold", hjust = 1),
    axis.title.y = element_text(size=18, face="bold")
  ) + 
  theme(#legend.position = c(0.15, 0.85),
    legend.position = 'none',
    legend.title = element_text(colour="black", size=12, 
                                face="bold"),
    legend.text = element_text(colour="black", size=12, 
                               face="bold"))


plot(p)

ggsave(filename="../Output/toxo_cdc/figures/tg_Intra_deg_numbers.pdf",
       plot=p,
       width = 6, height = 6,
       units = "in", # other options are "in", "cm", "mm"
       dpi = 300
)


## Enrichment
e.a.atac.trans <- read.xlsx('../Output/toxo_cdc/tabels/GO_toxo_db_atact_transition_points_sig_KZ.xlsx')

e.a.atac.trans <- e.a.atac.trans %>% dplyr::filter(Benjamini <= 0.1 & rank < 30) %>% arrange(Benjamini)

# e.a.sig$GO.name2 <- gsub("transport", "trans.",
#                          gsub("obsolete oxidation-reduction", "obs. ox-reduc.",
#                               gsub("symbiont-containing", "symb.-cont.", 
#                                    gsub("membrane", "mem.",
#                                         gsub("biosynthetic", "biosynth.", 
#                                              gsub("phosphorylation", "phosph.",
#                                                   gsub("activity", "act.",
#                                                        gsub("aminoacylation for prot. translation", "aa prot. trans.",
#                                                             gsub("synthesis coupled", "synth. coup.", 
#                                                                  gsub("heterodimerization", "hetero.dim", 
#                                                                       gsub("protein", "prot.", 
#                                                                            gsub("process", "proc.", 
#                                                                                 gsub("cellular", "cell.", 
#                                                                                      gsub("aminoacyl-tRNA ", "aa-tRNA ",
#                                                                                           gsub("complex", "comp.", 
#                                                                                                gsub("threonine-type ", "", 
#                                                                                                     gsub("proteolysis involved in cellular protein catabolic process", "proteolysis/catabolic proc.",
#                                                                                                          gsub(", alpha-subunit complex", "", e.a.sig$GO.name))))))))))))))))))
# 
# e.a.sig <- e.a.sig %>% mutate(GO.t = as.factor(GO.name2),
#                               name = reorder_within(GO.name2, by = qval, within = category))

e.a.atac.trans$Name <- factor(e.a.atac.trans$Name, levels = unique(e.a.atac.trans$Name[sort(e.a.atac.trans$Benjamini, index.return = T)$ix]))
p1 <- ggplot(e.a.atac.trans, aes(x = cluster, y = Name)) + 
  geom_point(aes(color = Benjamini, size = -log(Benjamini))) +
  theme_bw(base_size = 12) +
  scale_colour_gradient(limits=c(0, 0.01), low="red") +
  ylab(NULL) + xlab(NULL) +
  #facet_grid(category~., scales='free') + 
  theme(strip.text = element_text(size = 12, face="bold", angle = 90)) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8, face="bold")) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 8, face="bold")) +
  theme(strip.background = element_rect(colour="black", fill="white",
                                        size=0.5, linetype="solid")) +
  #facet_wrap(spp~.) + 
  #ggtitle(titles[i]) +
  theme(
    plot.title = element_text(size=10, face = "bold.italic", color = 'red'),
    axis.title.x = element_text(size=10, face="bold", hjust = 1),
    axis.title.y = element_text(size=10, face="bold")
  ) + 
  theme(#legend.position = c(0.15, 0.85),
    #legend.position = 'none',
    legend.title = element_text(colour="black", size=12, 
                                face="bold"),
    legend.text = element_text(colour="black", size=12, 
                               face="bold"))



plot(p1)

ggsave(filename="../Output/toxo_cdc/figures/enrichment_atac_transition.pdf",
       plot=p1,
       width = 6, height = 10,
       units = "in", # other options are "in", "cm", "mm"
       dpi = 300
)


e.a.phase <- read.xlsx('../Output/toxo_cdc/tabels/GO_toxo_db_phases_sig_KZ.xlsx')

e.a.phase$cluster <- factor(e.a.phase$cluster, levels = c('G1.a', 'G1.b', 'S', 'M', 'C'))
e.a.phase$Name <- factor(e.a.phase$Name, levels = unique(e.a.phase$Name[sort(e.a.phase$Benjamini, index.return = T)$ix]))
p2 <- ggplot(e.a.phase, aes(x = cluster, y = Name)) + 
  geom_point(aes(color = Benjamini, size = -log(Benjamini))) +
  theme_bw(base_size = 12) +
  scale_colour_gradient(limits=c(0, 0.01), low="red") +
  ylab(NULL) + xlab(NULL) +
  #facet_grid(category~., scales='free') + 
  theme(strip.text = element_text(size = 12, face="bold", angle = 90)) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8, face="bold")) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 8, face="bold")) +
  theme(strip.background = element_rect(colour="black", fill="white",
                                        size=0.5, linetype="solid")) +
  #facet_wrap(spp~.) + 
  #ggtitle(titles[i]) +
  theme(
    plot.title = element_text(size=10, face = "bold.italic", color = 'red'),
    axis.title.x = element_text(size=10, face="bold", hjust = 1),
    axis.title.y = element_text(size=10, face="bold")
  ) + 
  theme(#legend.position = c(0.15, 0.85),
    #legend.position = 'none',
    legend.title = element_text(colour="black", size=12, 
                                face="bold"),
    legend.text = element_text(colour="black", size=12, 
                               face="bold"))



plot(p2)

ggsave(filename="../Output/toxo_cdc/figures/enrichment_cell_cycle_markers.pdf",
       plot=p2,
       width = 6, height = 10,
       units = "in", # other options are "in", "cm", "mm"
       dpi = 300
)


## Pseudo time
sc.rna.genes.expr.pt <- readRDS('../Input/toxo_cdc/rds/sc_rna_genes_expr_pt.rds')
sc.atac.genes.expr.pt <- readRDS('../Input/toxo_cdc/rds/sc_atac_genes_expr_pt.rds')


S.O.rna  <- readRDS('../Input/toxo_cdc/rds/S.O_intra_lables_pt.rds')
S.O.atac <- readRDS('../Input/toxo_cdc/rds/S.O_intra_atac_lables_pt.rds')

rna.sds.data <- readRDS('../Input/toxo_cdc/rds/sc_rna_sds_data.rds')
atac.sds.data <- readRDS('../Input/toxo_cdc/rds/sc_atac_sds_data.rds')


L <- wiskerPlot(S.O.rna)

pdf(file="../Output/toxo_cdc/figures/whiskers_rna.pdf",
    width=6, height=6)


par(mar = c(5, 5, 4, 4) + 0.1)
plot(x = -35:10, y = -25:20, type = 'n', xlab = '', ylab = '',  xaxt = "n", yaxt = "n", axes=FALSE,
     lwd = 2, cex.lab = 1.5, cex.main = 2, cex.axis = 1.5)
whiskers(as.matrix(L$pc[,c(1,2)]), L$fit$s, col = "gray")
color = rep(NA, length=length(rna.sds.data$phase))
color[which(rna.sds.data$phase=="G1.a")] = "#f8766d"
color[which(rna.sds.data$phase=="G1.b")] = "#a3a400"
color[which(rna.sds.data$phase=="S")] = "#00bf7d"
color[which(rna.sds.data$phase=="M")] = "#04b0f6"
color[which(rna.sds.data$phase=="C")] = "#e76bf3"
points(rna.sds.data$PC_1, rna.sds.data$PC_2, cex = 0.5, col = color, pch = 20)
points(rna.sds.data$sc1[rna.sds.data$cell.ord],rna.sds.data$sc2[rna.sds.data$cell.ord], cex = 0.2, col = 'red')
# grid(nx = NULL, ny = NULL,
#      lty = 1, col = "gray", lwd = 1)
# 

dev.off()


L <- wiskerPlot(S.O.atac)

pdf(file="../Output/toxo_cdc/figures/whiskers_atac.pdf",
    width=6, height=6)


par(mar = c(5, 5, 4, 4) + 0.1)
plot(x = -35:10, y = -25:20, type = 'n',  xlab = '', ylab = '',  xaxt = "n", yaxt = "n", axes=FALSE,  
     lwd = 2, cex.lab = 1.5, cex.main = 2, cex.axis = 1.5)
whiskers(as.matrix(L$pc[,c(1,2)]), L$fit$s, col = "gray")
color = rep(NA, length=length(atac.sds.data$phase))
color[which(atac.sds.data$phase=="G1.a")] = "#f8766d"
color[which(atac.sds.data$phase=="G1.b")] = "#a3a400"
color[which(atac.sds.data$phase=="S")] = "#00bf7d"
color[which(atac.sds.data$phase=="M")] = "#04b0f6"
color[which(atac.sds.data$phase=="C")] = "#e76bf3"
points(atac.sds.data$PC_1, atac.sds.data$PC_2, cex = 0.5, col = color, pch = 20)
points(atac.sds.data$sc1[atac.sds.data$cell.ord],atac.sds.data$sc2[atac.sds.data$cell.ord], cex = 0.2, col = 'red')

dev.off()


picewise_scale <- function(sds.data){
  G <- c(0, 3) # 3h
  S <- c(3, 4.7) # 1.7h
  M <- c(4.7, 5) # ~20 min
  C <- c(5, 6) # 1h
  
  
  ind.G1.a <- which(sds.data$phase == 'G1.a')
  ind.G1.b <- which(sds.data$phase == 'G1.b')
  ind.S <- which(sds.data$phase == 'S')
  ind.M <- which(sds.data$phase == 'M')
  ind.C <- which(sds.data$phase == 'C')
  
  t <- sds.data$pt
  
  t0 <- 0
  t1 <- quantile(t[ind.G1.b], prob=0.25) ## Start of G1.b
  t2 <- (quantile(t[ind.G1.b], prob=0.75) + 
           quantile(t[ind.S], prob=0.25)) / 2 ## End of G1.b, Start of S
  t3 <- (quantile(t[ind.S], prob=0.75) + 
           quantile(t[ind.M], prob=0.25)) / 2## End of S, Start of M
  t4 <- (quantile(t[ind.M], prob=0.75) + 
           quantile(t[ind.C], prob=0.25)) / 2## End of M, Start of C
  t5 <- 6
  
  slp.g <- (G[2] - G[1]) / (t2 - t0)
  inc.g  <- c(t0, G[1])
  
  slp.s <- (S[2] - S[1]) / (t3 - t2)
  inc.s  <- c(t2, S[1])
  
  slp.m <- (M[2] - M[1]) / (t4 - t3)
  inc.m  <- c(t3, M[1])
  
  slp.c <- (C[2] - C[1]) / (t5 - t4)
  inc.c  <- c(t4, C[1])
  
  s.t <- data.frame(pt = seq(0, 6, 0.01))  
  s.t <- s.t %>% mutate(phase = case_when(pt >= t0 & pt < t1 ~ 'G1.a',
                                          pt >= t1 & pt < t2 ~ 'G1.b',
                                          pt >= t2 & pt < t3 ~ 'S',
                                          pt >= t3 & pt < t4 ~ 'M',
                                          pt >= t4 ~ 'C'), 
                        t = case_when(pt >= t0 & pt < t2 ~ inc.g[2] + slp.g * (pt - inc.g[1]),
                                      pt >= t2 & pt < t3 ~ inc.s[2] + slp.s * (pt - inc.s[1]),
                                      pt >= t3 & pt < t4 ~ inc.m[2] + slp.m * (pt - inc.m[1]),
                                      pt >= t4 ~ inc.c[2] + slp.c * (pt - inc.c[1])))
  
  return(s.t)
}

rna.s.t <- picewise_scale(rna.sds.data)

rna.s.t$phase <- factor(rna.s.t$phase, levels = c('G1.a', 'G1.b', 'S', 'M', 'C'))
p1  <- ggplot(rna.s.t, aes(x= pt,y=t)) +
  geom_line(aes(#fill = lable.prob,
    color = phase, 
  ), #color = 'blue', 
  alpha = 0.9,
  size = 1) + 
  geom_segment(aes(x = pt[phase == 'G1.b'][1], y = 0, 
                   xend = pt[phase == 'G1.b'][1], 
                   yend = t[phase == 'G1.b'][1]), linetype=2, color = 'black') +
  geom_segment(aes(x = 0, y = t[phase == 'G1.b'][1], 
                   xend = pt[phase == 'G1.b'][1], 
                   yend = t[phase == 'G1.b'][1]), linetype=2, color = 'black') +
  geom_segment(aes(x = pt[phase == 'S'][1], y = 0, 
                   xend = pt[phase == 'S'][1], 
                   yend = t[phase == 'S'][1]), linetype=2, color = 'black') +
  geom_segment(aes(x = 0, y = t[phase == 'S'][1], 
                   xend = pt[phase == 'S'][1], 
                   yend = t[phase == 'S'][1]), linetype=2, color = 'black') +
  geom_segment(aes(x = pt[phase == 'M'][1], y = 0, 
                   xend = pt[phase == 'M'][1], 
                   yend = t[phase == 'M'][1]), linetype=2, color = 'black') +
  geom_segment(aes(x = 0, y = t[phase == 'M'][1], 
                   xend = pt[phase == 'M'][1], 
                   yend = t[phase == 'M'][1]), linetype=2, color = 'black') +
  geom_segment(aes(x = pt[phase == 'C'][1], y = 0, 
                   xend = pt[phase == 'C'][1], 
                   yend = t[phase == 'C'][1]), linetype=2, color = 'black') +
  geom_segment(aes(x = 0, y = t[phase == 'C'][1], 
                   xend = pt[phase == 'C'][1], 
                   yend = t[phase == 'C'][1]), linetype=2, color = 'black') +
  geom_segment(aes(x = 6, y = 0, 
                   xend = 6, 
                   yend = 6), linetype=2, color = 'black') +
  geom_segment(aes(x = 0, y = 6, 
                   xend = 6, 
                   yend = 6), linetype=2, color = 'black') +
  
  
  theme_bw(base_size = 14) +
  theme(legend.position = "right") +
  #scale_fill_gradientn(colours = viridis::inferno(10)) +
  #scale_fill_gradientn(colours = col_range(10)) +
  #scale_fill_gradient(low = "gray66", high = "blue", midpoint = mid_point) + 
  #scale_fill_brewer(palette = "BuPu") +
  ylab('Scaled time') + xlab('Pseudo time') +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 16, face="bold")) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 16, face="bold")) +
  theme(strip.background = element_rect(colour="black", fill="white",
                                        size=0.5, linetype="solid")) +
  
  scale_y_continuous(breaks=c(0, 1, 2, 
                              round(rna.s.t$t[rna.s.t$phase == 'G1.b'][1], 1), 
                              3, 4,
                              round(rna.s.t$t[rna.s.t$phase == 'M'][1], 1), 5, 6)) + 
  scale_x_continuous(breaks=c(0, 1, 2, round(rna.s.t$pt[rna.s.t$phase == 'G1.b'][1], 1), 
                              round(rna.s.t$pt[rna.s.t$phase == 'S'][1], 1), 4, 
                              round(rna.s.t$pt[rna.s.t$phase == 'M'][1], 1),
                              round(rna.s.t$pt[rna.s.t$phase == 'C'][1], 1), 5, 6)) + 
  
  #facet_wrap(spp~.) + 
  #ggtitle(titles[i]) +
  theme(
    plot.title = element_text(size=14, face = "bold.italic", color = 'red'),
    axis.title.x = element_text(size=18, face="bold", hjust = 1),
    axis.title.y = element_text(size=18, face="bold")
  ) + 
  theme(#legend.position = c(0.9, 0.3),
        legend.position = 'none',
        legend.title = element_text(colour="black", size=12, 
                                    face="bold"),
        legend.text = element_text(colour="black", size=12, 
                                   face="bold"))


plot(p1)

ggsave(filename="../Output/toxo_cdc/figures/ps_time_piecewise_scale.pdf",
       plot=p1,
       width = 6, height = 6,
       units = "in"
)

## Proportions
ss <- rna_sub@meta.data %>% group_by(phase) %>% summarise(total = n()) %>% 
  mutate(prop = total/nrow(rna_sub@meta.data))

ss$phase <- factor(ss$phase, levels = c('G1.a', 'G1.b', 'S', 'M', 'C'))
ss$phase0 <- gsub('\\..*', '', ss$phase)
ss$phase0 <- factor(ss$phase0, levels = c('G1', 'S', 'M', 'C'))
ss$label <- ifelse(ss$phase0 == 'G1', paste(ss$phase, round(ss$prop, 2), sep = ' : '), round(ss$prop, 2))
p <- ggplot(data=ss, aes(x=phase0, y=prop, fill=phase)) +
  geom_bar(stat ="identity",
           position = "stack")+
  geom_text(aes(label=label), vjust=1.6, color="black", size=5, fontface="bold")+
  theme_minimal() + 
  ylab('') + xlab('') +
  theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 16, face="bold")) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 16, face="bold")) +
  theme(strip.background = element_rect(colour="black", fill="white",
                                        size=0.5, linetype="solid")) +
  #facet_wrap(spp~.) + 
  #ggtitle(titles[i]) +
  theme(
    plot.title = element_text(size=18, face = "bold.italic", color = 'red'),
    axis.title.x = element_text(size=18, face="bold", hjust = 1),
    axis.title.y = element_text(size=18, face="bold")
  ) + 
  theme(#legend.position = c(0.15, 0.85),
    legend.position = 'none',
    legend.title = element_text(colour="black", size=12, 
                                face="bold"),
    legend.text = element_text(colour="black", size=12, 
                               face="bold"))



plot(p)

ggsave(filename="../Output/toxo_cdc/figures/phase_proportions.pdf",
       plot=p,
       width = 6, height = 6,
       units = "in", # other options are "in", "cm", "mm"
       dpi = 300
)


## Heatmaps
sc.rna.mu.scale <- readRDS('../Input/toxo_cdc/rds/sc_rna_spline_mu_scale_phase.rds')
sc.atac.mu.scale <- readRDS('../Input/toxo_cdc/rds/sc_atac_spline_mu_scale_phase.rds')

## Filter for Marker genes
sc.rna.mu.scale <- sc.rna.mu.scale[sc.rna.mu.scale$GeneID %in% marker.genes$gene,]
sc.atac.mu.scale <- sc.atac.mu.scale[sc.atac.mu.scale$GeneID %in% marker.genes$gene,]


p1 <- ggplot(sc.rna.mu.scale, aes(x = x, y = GeneID, fill = expr)) + 
  geom_tile() + 
  facet_grid(phase~., scales = "free",  space='free',
             labeller=label_wrap_gen(multi_line = TRUE))+
  theme(strip.text = element_text(size = 14, face="bold", angle = 90)) + 
  theme(strip.background = element_rect(colour="black", fill="white",
                                        size=0.2, linetype="solid")) +
  
  ylab("") + xlab("") + ggtitle('') + 
  scale_x_continuous(expand = expansion(mult = c(0, 0))) + 
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


ggsave(filename="../Output/toxo_cdc/figures/scRNA_heatmap_phase.png",
       plot=p1,
       width = 6, height = 7,
       units = "in", # other options are "in", "cm", "mm"
       dpi = 300
)


p2 <- ggplot(sc.atac.mu.scale, aes(x = x, y = GeneID, fill = expr)) + 
  geom_tile() + 
  facet_grid(phase~., scales = "free",  space='free',
             labeller=label_wrap_gen(multi_line = TRUE))+
  theme(strip.text = element_text(size = 14, face="bold", angle = 90)) + 
  theme(strip.background = element_rect(colour="black", fill="white",
                                        size=0.2, linetype="solid")) +
  
  ylab("") + xlab("") + ggtitle('') + 
  scale_x_continuous(expand = expansion(mult = c(0, 0))) + 
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


ggsave(filename="../Output/toxo_cdc/figures/scATAC_heatmap_phase.png",
       plot=p2,
       width = 6, height = 7,
       units = "in", # other options are "in", "cm", "mm"
       dpi = 300
)



## Cross correlation
cc.dat <- readRDS('../Input/toxo_cdc/rds/sc_rna_sc_atac_cross_cor_lag.rds')



p2 <- ggplot(cc.dat, aes(x = ccs)) + 
  geom_histogram(aes(y = ..density..),
                 colour = 1, fill = "steelblue") + 
  geom_density(lwd = 1.2,
               linetype = 1,
               colour = 2) + 
  
  theme_bw(base_size = 14) +
  theme(legend.position = "right") +
  #scale_fill_gradientn(colours = viridis::inferno(10)) +
  #scale_fill_gradientn(colours = col_range(10)) +
  #scale_fill_gradient(low = "gray66", high = "blue", midpoint = mid_point) + 
  #scale_fill_brewer(palette = "BuPu") +
  xlab('scRNA/scATAC cross corr.') + ylab('density') +
  theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 16, face="bold")) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 16, face="bold")) +
  theme(strip.background = element_rect(colour="black", fill="white",
                                        size=0.5, linetype="solid")) +
  #facet_wrap(spp~.) + 
  #ggtitle(titles[i]) +
  theme(
    plot.title = element_text(size=18, face = "bold.italic", color = 'red'),
    axis.title.x = element_text(size=18, face="bold", hjust = 1),
    axis.title.y = element_text(size=18, face="bold")
  ) + 
  theme(#legend.position = c(0.15, 0.85),
    legend.position = 'none',
    legend.title = element_text(colour="black", size=12, 
                                face="bold"),
    legend.text = element_text(colour="black", size=12, 
                               face="bold"))


plot(p2)



ggsave(filename="../Output/toxo_cdc/figures/sc_rna_sc_atac_cross_corr.pdf",
       plot=p2,
       width = 6, height = 6,
       units = "in"
)

## AP2 clusters
sc.rna.sc.atac.joint.long <- readRDS('../Input//toxo_cdc/rds/AP2_sc_rna_sc_atac_joint_dtw_clust.rds')
p1  <- ggplot(sc.rna.sc.atac.joint.long, aes(x= time,y=normExpr)) +
  geom_path(aes(color = Name,),alpha = 0.8, size = 1)+ 
  theme_bw(base_size = 14) +
  theme(legend.position = "right") +
  ylab('normExpr') + xlab('Time') +
  theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 12, face="bold")) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 12, face="bold")) +
  theme(strip.background = element_rect(colour="black", fill="white",
                                        size=0.5, linetype="solid")) +
  theme(strip.text = element_text(size = 14, face="bold", angle = 0)) + 
  
  #coord_cartesian(xlim = c(0,6.5)) + 
  geom_text_repel(aes(label = label), size = 2.5, fontface = "bold",
                  box.padding = unit(0.6, "lines"),
                  max.overlaps = 300,
                  #segment.angle = 180,
                  nudge_x = 0.25, 
                  nudge_y = 0.25,
                  hjust=0.25,
                  #nudge_x=0.25, 
                  segment.size = 0.1,
                  na.rm = TRUE)+ 
  facet_grid(cluster.RNA~data, scales = 'free', space = 'free') +
  
  
  #ggtitle(titles[i]) +
  theme(
    plot.title = element_text(size=14, face = "bold.italic", color = 'red'),
    axis.title.x = element_text(size=14, face="bold", hjust = 1),
    axis.title.y = element_text(size=14, face="bold")
  ) + 
  theme(#legend.position = c(0.15, 0.85),
    legend.position = 'none',
    legend.title = element_text(colour="black", size=12, 
                                face="bold"),
    legend.text = element_text(colour="black", size=12, 
                               face="bold")) 


plot(p1)

ggsave(filename="../Output/toxo_cdc/figures/AP2_rna_atac_clusters.pdf",
       plot=p1,
       width = 12, height = 8,
       units = "in"
)


### ATAC coverage plot
Tg_ATAC <- readRDS('../Input/toxo_cdc/rds/S.O_ATAC_peak.rds')

region <- StringToGRanges(regions = gsub('TGME49-', 'TGME49_', top.da$gene[4]),  sep = c("-", "-"))
#region <- StringToGRanges(regions = gsub('TGME49-', 'TGME49_', rownames(da_peaks)[5]),  sep = c("-", "-"))

xx <- findOverlaps(region, tx_trans)
tx_trans[xx@to]$gene_id

my.gene <- tx_trans[xx@to]$gene_id[1]

Idents(Tg_ATAC) <- 'phase'
Tg_ATAC@active.ident <- factor(Tg_ATAC@active.ident, levels = c('G1.a', 'G1.b', 'S', 'M', 'C'))
DefaultAssay(Tg_ATAC) <- 'peaks'
p2 <- CoveragePlot(
  object = Tg_ATAC,
  sep = c("-", "-"),
  #region = gsub('TGME49-', 'TGME49_', rownames(Tg_ATAC)[1:3]),
  region = region,
  #region = gsub('TGME49-', 'TGME49_', VariableFeatures(Tg_ATAC)[20]),
  extend.upstream = 4000,
  extend.downstream = 2000
)

plot(p2)
ggsave(filename="../Output/scClockFigs/RON2_track_pileup_by_phase.pdf",
       plot=p2,
       width = 5, height = 6,
       units = "in", # other options are "in", "cm", "mm"
       dpi = 300
)


