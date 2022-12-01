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


## Raw Count files sparse (cewll x gene matrix)
## Intra cellular Tachyzoiutes
intra.file.csv <- "../Input/toxo_scRNA_MJ/RH.intra.expr.csv"

## IDs (Converting from GT1 to ME49)
prod.desc  <- read.xlsx('../Input/toxo_genomics/genes/ProductDescription_GT1.xlsx')
TGGT1_ME49 <- read.xlsx('../Input/toxo_genomics/Orthologs/TGGT1_ME49 Orthologs.xlsx')


## Get the expression matrix and convert the ids to ME49
getExpr <- function(in.file, TGGT1_ME49){
  file.counts <- read.csv(in.file)
  genes <- file.counts$X
  ind <- which(genes %in% TGGT1_ME49$TGME49)
  file.counts <- file.counts[ind, ]
  genes <- genes[ind]
  genes <- TGGT1_ME49$TGGT1[match(genes, TGGT1_ME49$TGME49)]
  expr <- file.counts[,-1]
  rownames(expr) <- genes
  
  return(expr)
}

intra.counts <- getExpr(intra.file.csv, TGGT1_ME49)


## Create the Seurat Object and filter cells/genes 
feats <- c("nFeature_RNA","nCount_RNA")
S.O <- CreateSeuratObject(counts = intra.counts)
S.O$orig.ident <- 'intra'
VlnPlot(S.O, features = feats, pt.size = 0.1,ncol = 2) + NoLegend()
FeatureScatter(S.O, "nCount_RNA", "nFeature_RNA",  pt.size = 0.5)
selected_c <- WhichCells(S.O, expression = nFeature_RNA > 100 & nFeature_RNA < 950)
selected_f <- rownames(S.O)[ Matrix::rowSums(S.O) > 5]
S.O  <- subset(S.O, features=selected_f, cells=selected_c)


## Downsample to 8000 cells for processing speed
S.O <- subset(x = S.O, downsample = 8000)


## This function processes the Seurat object to generate clusters (graph-based k-means)
## generate PCA and UMAP
prep_S.O <- function(S.O, res = 0.1){
  S.O <- NormalizeData(S.O, normalization.method = "LogNormalize", scale.factor = 10000)
  S.O <- FindVariableFeatures(S.O, selection.method = "vst", nfeatures = 6000)

  all.genes <- rownames(S.O)
  S.O <- ScaleData(S.O, features = all.genes)
  S.O <- RunPCA(S.O, features = VariableFeatures(object = S.O))
  S.O <- FindNeighbors(S.O, dims = 1:10, reduction = 'pca')
  S.O <- FindClusters(S.O, resolution = res)
  S.O <- RunUMAP(S.O, dims = 1:13, n.components = 3L)
  return(S.O)
}



# transfer cell cycle phase label Boothroyed data
S.O.tg.boothroyd <- readRDS('../Input/boothroyd_sc_all_data/rds/S.O.tg_RH_boothroyd.rds')
S.O <- prep_S.O(S.O, res = 0.2)

anchors <- FindTransferAnchors(reference = S.O.tg.boothroyd, query = S.O, dims = 1:30)
predictions <- TransferData(anchorset = anchors, refdata = S.O.tg.boothroyd@meta.data$phase,dims = 1:30)
predictions$phase <- predictions$predicted.id
#predictions$phase[which(predictions$prediction.score.max < 0.7)] <- 'NA'
S.O <- AddMetaData(object = S.O, metadata = predictions)



## Test plots
Idents(S.O) <- 'phase'
p <- DimPlot(S.O, reduction = "pca", 
             #group.by = "cell", 
             #split.by = 'spp',
             pt.size = 1,
             #shape.by='spp',
             label = F, label.size = 4) + #NoLegend() + 
  theme(panel.spacing = unit(0.5, "lines")) + 
  theme(axis.text.x = element_text(face="bold", size=12, angle=0)) +
  theme(axis.text.y = element_text(face="bold", size=12, angle=0)) +
  theme(
    axis.title.x = element_text(size=14, face="bold"),
    axis.title.y = element_text(size=14, face="bold")
  )


plot(p)

## To get the meta data:
S.O@meta.data

##Change the identity of cells to seurate generated clusters
Idents(S.O) <- 'seurat_clusters'
p <- DimPlot(S.O, reduction = "pca", 
             #group.by = "cell", 
             #split.by = 'spp',
             pt.size = 1,
             #shape.by='spp',
             label = F, label.size = 4) + #NoLegend() + 
  theme(panel.spacing = unit(0.5, "lines")) + 
  theme(axis.text.x = element_text(face="bold", size=12, angle=0)) +
  theme(axis.text.y = element_text(face="bold", size=12, angle=0)) +
  theme(
    axis.title.x = element_text(size=14, face="bold"),
    axis.title.y = element_text(size=14, face="bold")
  )


plot(p)

## Get the PCA
pc <- S.O[['pca']]@cell.embeddings
pc <- data.frame(pc) %>% dplyr::mutate(Sample = rownames(pc)) 
pc$cluster <- S.O$seurat_clusters


## Save the data
saveRDS(S.O, '../Input/toxo_cdc/rds/S.O.intra_lables.rds')
