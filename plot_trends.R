library(Seurat)
library(ggplot2)
library(ggfortify)
library(jcolors)
require(gridExtra)
library(grid)
library(gam)
library(princurve)
library(parallel)
library(tidyverse)
library(MyEllipsefit)
library(sctransform)
library(openxlsx)
library(doParallel)
library(tidytext)
library(ggrepel)
library(dtwclust)
library(geomtextpath)


source('./util_funcs.R')

## IDs
prod.desc  <- read.xlsx('../Input/toxo_genomics/genes/ProductDescription_GT1.xlsx')
TGGT1_ME49 <- read.xlsx('../Input/toxo_genomics/Orthologs/TGGT1_ME49 Orthologs.xlsx')

## scDATA

rna_sub <- readRDS('../Input/toxo_cdc/rds/S.O_intra_lables_pt.rds')
atac_sub <- readRDS('../Input/toxo_cdc/rds/S.O_intra_atac_lables_pt.rds')


## Splines
sc.rna.spline.fits <- readRDS('../Input/toxo_cdc/rds/sc_rna_spline_fits_all_genes.rds')
sc.atac.spline.fits <- readRDS('../Input/toxo_cdc/rds/sc_atac_spline_fits_all_genes.rds')

## Turn the data into wide format (time by gene) and center & scale each gene
sc.rna.dtw.wide <- sc.rna.spline.fits %>% 
  pivot_wider(names_from = 'GeneID', values_from = 'y') %>% 
  mutate_at(vars(matches("TGME")), ~scale(., center = T, scale = T)) %>%
  as.data.frame()

sc.atac.dtw.wide <- sc.atac.spline.fits %>% 
  pivot_wider(names_from = 'GeneID', values_from = 'y') %>% 
  mutate_at(vars(matches("TGME")), ~scale(., center = T, scale = T)) %>%
  as.data.frame()


sc.rna.mu.scale <- sc.rna.dtw.wide %>% 
  pivot_longer(-x, names_to = 'GeneID', values_to = 'expr')

sc.atac.mu.scale <- sc.atac.dtw.wide %>% 
  pivot_longer(-x, names_to = 'GeneID', values_to = 'expr')




plot_trends <- function(my.GeneID){
  
  my.rna <- sc.rna.mu.scale %>% dplyr::filter(GeneID == my.GeneID)
  my.atac <- sc.atac.mu.scale %>% dplyr::filter(GeneID == my.GeneID)
  
  p1  <- ggplot(my.rna , aes(x= x,y=expr)) +
    geom_line(color = 'blue',alpha = 0.8, size = 0.8)+ 
    theme_bw(base_size = 14) +
    theme(legend.position = "right") +
    ylab('normExpr') + xlab('Time') +
    theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 12, face="bold")) +
    theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 12, face="bold")) +
    theme(strip.background = element_rect(colour="black", fill="white",
                                          size=0.5, linetype="solid")) +
    ggtitle(paste('rna', my.GeneID)) + 
    theme(strip.text = element_text(size = 14, face="bold", angle = 0)) + 
    theme(
      plot.title = element_text(size=14, face = "bold.italic", color = 'black'),
      axis.title.x = element_text(size=14, face="bold", hjust = 1),
      axis.title.y = element_text(size=14, face="bold")
    ) 
  
  
  p2  <- ggplot(my.atac , aes(x= x,y=expr)) +
    geom_line(color = 'red',alpha = 0.8, size = 0.8)+ 
    theme_bw(base_size = 14) +
    theme(legend.position = "right") +
    ylab('normExpr') + xlab('Time') +
    theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 12, face="bold")) +
    theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 12, face="bold")) +
    theme(strip.background = element_rect(colour="black", fill="white",
                                          size=0.5, linetype="solid")) +
    ggtitle(paste('atac', my.GeneID)) + 
    theme(strip.text = element_text(size = 14, face="bold", angle = 0)) + 
    theme(
      plot.title = element_text(size=14, face = "bold.italic", color = 'black'),
      axis.title.x = element_text(size=14, face="bold", hjust = 1),
      axis.title.y = element_text(size=14, face="bold")
    ) 
  
  p <- grid.arrange(p1, p2)
  
  return(p)
}

conserved_TFs <- c('TGME49-286710',
'TGME49-206650',
'TGME49-236840',
'TGME49-252310',
'TGME49-203380',
'TGME49-201220',
'TGME49-242320')


for(i in 1:length(conserved_TFs)){
  
  p <- plot_trends(conserved_TFs[i])
  f.n <- paste("../Output/toxo_cdc/figures/conserved_TFs/", conserved_TFs[i], '.pdf', sep = '')
  ggsave(filename=f.n,
         plot=p,
         width = 6, height = 6,
         units = "in", # other options are "in", "cm", "mm"
         dpi = 300
  )
  
  
}

p <- plot_trends('TGME49-305330')

plot(p)

FeaturePlot(rna_sub, 'TGME49-305330', reduction = 'pca')



Crk2 <- 'TGGT1-218220'
Ark3 <- 'TGGT1-203010'
CycP1 <- 'TGGT1-267580'
CycP2 <- 'TGGT1-305330'
Cyc5 <- 'TGGT1-293280'
my.features <- c(Crk2, Ark3, CycP1, CycP2, Cyc5)

for(i in 1:length(my.features)){
  
  p <- plot_trends(conserved_TFs[i])
  f.n <- paste("../Output/toxo_cdc/figures/cyclins/", my.features[i], '.pdf', sep = '')
  ggsave(filename=f.n,
         plot=p,
         width = 6, height = 6,
         units = "in", # other options are "in", "cm", "mm"
         dpi = 300
  )
  
  
}

