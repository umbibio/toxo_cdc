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
library(Matrix)
library(graphlayouts)
library(Signac)
library(Seurat)
library(patchwork)
library(hdf5r)
library(GenomeInfoDbData)
library(GenomicRanges)
library(GenomicAlignments)
library(Biostrings)
library(rtracklayer)
library(GenomicFeatures)


source('./util_funcs.R')



source('./util_funcs.R')


## Generating count data
input.dir <-  "../Input/toxo_scRNA_MJ_ME49_59/27-30-33-hpi_with_10M_RH-YFP2/"
barcode.path <- paste(input.dir, "filtered_feature_bc_matrix/barcodes.tsv.gz", sep = "")
features.path <- paste(input.dir, "filtered_feature_bc_matrix/features.tsv.gz", sep = "")
matrix.path <- paste(input.dir, "filtered_feature_bc_matrix/matrix.mtx.gz", sep = "")

mat <- readMM(file = matrix.path)

feature.names = read.delim(features.path,
                           header = FALSE,
                           stringsAsFactors = FALSE)

barcode.names = read.delim(barcode.path,
                           header = FALSE,
                           stringsAsFactors = FALSE)

colnames(mat) = barcode.names$V1
rownames(mat) = feature.names$V1

expr <- as.data.frame(as.matrix(mat))

write.csv(expr, "../Input/toxo_scRNA_MJ_ME49_59/RH.intra.expr.csv")




num.cores <- detectCores(all.tests = FALSE, logical = TRUE)

## Tw of the datasets are mutant lines.
## ID of the KO genes
Crk2 <- 'TGGT1-218220'
Ark3 <- 'TGGT1-203010'

## Count files
intra.file.csv <- "../Input/toxo_scRNA_MJ_ME49_59/RH.intra.expr.csv" ## New version of the Genome

## IDs
prod.desc  <- read.xlsx('../Input/toxo_genomics/genes/ProductDescription_GT1.xlsx')
TGGT1_ME49 <- read.xlsx('../Input/toxo_genomics/Orthologs/TGGT1_ME49 Orthologs.xlsx')


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

## individual Seurat objects
feats <- c("nFeature_RNA","nCount_RNA")

# Intra
S.O.intra <- CreateSeuratObject(counts = intra.counts)
S.O.intra$orig.ident <- 'intra'
VlnPlot(S.O.intra, features = feats, pt.size = 0.1,ncol = 2) + NoLegend()
FeatureScatter(S.O.intra, "nCount_RNA", "nFeature_RNA",  pt.size = 0.5)
selected_c <- WhichCells(S.O.intra, expression = nFeature_RNA > 100 & nFeature_RNA < 1200)
selected_f <- rownames(S.O.intra)[ Matrix::rowSums(S.O.intra) > 5]
selected_f <- unique(c(selected_f, Crk2, Ark3)) ## Make sure Crk2 and Ark 3 are selected
S.O.intra  <- subset(S.O.intra, features=selected_f, cells=selected_c)

saveRDS(S.O.intra, '../Input/toxo_cdc/rds/ME49_59/S.O_intra_not_down_sample.rds')

set.seed(100)
S.O <- subset(x = S.O.intra, downsample = 8000)


## Individually process the data and transfer labels
## Boothroyed data
S.O.tg.boothroyd <- readRDS('../Input/boothroyd_sc_all_data/rds/S.O.tg_RH_boothroyd.rds')

S.O <- prep_S.O(S.O, res = 0.4)
anchors <- FindTransferAnchors(reference = S.O.tg.boothroyd, query = S.O, dims = 1:30)
predictions <- TransferData(anchorset = anchors, refdata = S.O.tg.boothroyd@meta.data$phase,dims = 1:30)
predictions$phase <- predictions$predicted.id
#predictions$phase[which(predictions$prediction.score.max < 0.7)] <- 'NA'
S.O <- AddMetaData(object = S.O, metadata = predictions)

saveRDS(S.O, '../Input/toxo_cdc/rds/ME49_59/S.O.intra_lables.rds')


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

########## ATAC-Seq


## Read scRAN-Seq data
S.O <- readRDS('../Input/toxo_cdc/rds/ME49_59/S.O.intra_lables.rds')

## IDs
prod.desc  <- read.xlsx('../Input/toxo_genomics/genes/ProductDescription_GT1.xlsx')
TGGT1_ME49 <- read.xlsx('../Input/toxo_genomics/Orthologs/TGGT1_ME49 Orthologs.xlsx')

## Map to ME49
counts = S.O@assays$RNA@counts
rownames(counts) <- gsub('_', '-', TGGT1_ME49$TGME49[match(gsub('-', '_', rownames(counts)), TGGT1_ME49$TGGT1)])

S.O.ME49 <- CreateSeuratObject(counts = counts)
S.O.ME49$orig.ident <- 'scRNA'
S.O.ME49 <- AddMetaData(S.O.ME49, S.O@meta.data)
Idents(S.O.ME49) <- 'phase'

S.O.ME49@meta.data$phase <- factor(S.O.ME49@meta.data$phase, levels = c('G1.a', 'G1.b', 'S', 'M', 'C'))
S.O.ME49 <- prep_S.O(S.O.ME49)
Idents(S.O.ME49) <- 'phase'
DimPlot(S.O.ME49, reduction = 'pca')


## Now read scATAC data
ME49.fasta <- readDNAStringSet("../Input/toxo_genomics/genome/ToxoDB-59_TgondiiME49_Genome.fasta")
chrs <- names(ME49.fasta)[grep("TGME49_chr", names(ME49.fasta))]

chr.len <- data.frame(chr = gsub(" ", "", unlist(lapply(strsplit(chrs, split = '\\|'), `[[`, 1))),
                      len = as.numeric(gsub('length=', '', unlist(lapply(strsplit(chrs, split = '\\|'), `[[`, 4)))))

txdb <- makeTxDbFromGFF(file="../Input/toxo_genomics/genome/ToxoDB-59_TgondiiME49_filter.gtf",
                        dataSource="Toxodb",
                        organism="Toxoplasma")

trans_biotypes <- select(txdb, keys=keys(txdb, "TXID"), 
                         columns = "TXTYPE", keytype =  "TXID")


genome(txdb) <- 'ME49'
#tx_genes <- genes(txdb)

#tx_trans <- unlist(transcriptsBy(txdb, by = c("gene", "exon", "cds")))

tx_trans <- exonsBy(txdb, by = "tx", use.names = TRUE)
tx_names <- names(tx_trans)
num.exons <- lapply(tx_trans, function(x) length(x))
tx_names <- rep(tx_names, unlist(num.exons))
tx_trans <- unlist(tx_trans)
tx_trans$tx_id <- tx_names
tx_trans$gene_id <- gsub('-t.*', '', tx_trans$tx_id)
tx_trans$gene_name <- tx_trans$gene_id
tx_trans$type <- 'exon'
tx_trans$gene_biotype <- 'protein_coding'
tx_trans$exon_name <- tx_trans$exon_rank

tmp <- chr.len$len[match(names(seqlengths(txdb)), chr.len$chr)]
names(tmp) <- names(seqlengths(txdb))
#seqlengths(tx_genes) <- tmp
seqlengths(tx_trans) <- tmp
#seqlevels(tx_genes)
seqlevels(tx_trans)
#inds <- c(c(5,6,1,2,3,7,8,10,11,9,4,12,13,14) + 100, seq(1, 100, by = 1))
inds <- c(5,6,1,2,3,7,8,10,11,9,4,12,13,14) 

#seqlevels(tx_genes) <- gsub('TGME49_', '', seqlevels(tx_genes)[inds])
#seqlevels(tx_genes) <- seqlevels(tx_genes)[inds]
#isCircular(tx_genes) <- rep(F, length(isCircular(tx_genes)))

seqlevels(tx_trans) <- seqlevels(tx_trans)[inds]
isCircular(tx_trans) <- rep(F, length(isCircular(tx_trans)))

#seqinfo(tx_genes)
seqinfo(tx_trans)

counts <- Read10X_h5(filename = "../Input/toxo_scATAC_MJ_ME49_59/filtered_peak_bc_matrix.h5")
metadata <- read.csv(
  file = "../Input/toxo_scATAC_MJ_ME49_59/singlecell.csv",
  header = TRUE,
  row.names = 1
)

metadata.filt <- metadata
metadata.filt$Sample <- rownames(metadata.filt)
metadata.filt <- metadata.filt[metadata.filt$Sample %in% colnames(counts), ]
peak_anno <- read_tsv("../Input/toxo_scATAC_MJ_ME49_59/peaks.bed", col_names = c('Chr', 'strt', 'stp'))

#counts <- Read10X_h5(filename = "../Input/scATAC/ME49_cell_ranger/filtered_peak_bc_matrix.h5")
#peak_anno <- read_tsv("../Input/scATAC/ME49_cell_ranger/peak_annotation.tsv")

chrom_assay <- CreateChromatinAssay(
  counts = counts,
  sep = c(":", "-"),
  genome = seqinfo(tx_trans),
  fragments = '../Input/toxo_scATAC_MJ_ME49_59/fragments.tsv.gz',
  min.cells = 5,
  min.features = 100
)

Tg_ATAC <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "peaks",
  meta.data = metadata.filt
)


Tg_ATAC[['peaks']]
granges(Tg_ATAC)
#annotations <- tx_genes
annotations <- tx_trans

#annotations$gene_biotype <- 'protein_coding'
#annotations$gene_name <- annotations$gene_id
#annotations$gene_name <- gsub('-t.*', '', annotations$tx_name)
#annotations$gene_id <- annotations$gene_name
#annotations$type <- 'exon'

Annotation(Tg_ATAC) <- annotations
Tg_ATAC <- NucleosomeSignal(object = Tg_ATAC)
Tg_ATAC <- TSSEnrichment(object = Tg_ATAC, fast = FALSE)

peaks <- CallPeaks(
  object = Tg_ATAC,
  macs2.path = "/Users/kouroshz/miniconda3/envs/macs2/bin/macs2",
  extsize = 100,
  additional.args = "--nomodel -B --SPMR"
)

## Save peaks
peaks.dat <- as.data.frame(peaks)
peaks.dat <- peaks.dat[grep('TGME49', peaks.dat$seqnames), ] %>% 
  dplyr::filter(10^(-neg_log10qvalue_summit) < 0.1) %>% 
  dplyr::select(seqnames, start, end, width, strand,  score, fold_change, contains('log')) 
colnames(peaks.dat)[1] <- 'chr'

saveRDS(peaks.dat, '../Input/toxo_cdc/rds/ME49_59/sc_atac_peaks_macs2.rds')
##

# fragments <- CreateFragmentObject(
#   path = "../Input/scATAC/ME49_cell_ranger/fragments.tsv.gz",
#   cells = colnames(Tg_ATAC),
#   validate.fragments = FALSE
# )
# #> Computing hash
# Fragments(Tg_ATAC) <- fragments





# add blacklist ratio and fraction of reads in peaks
Tg_ATAC$pct_reads_in_peaks <- Tg_ATAC$peak_region_fragments / Tg_ATAC$passed_filters * 100
Tg_ATAC$blacklist_ratio <- Tg_ATAC$blacklist_region_fragments / Tg_ATAC$peak_region_fragments

Tg_ATAC$high.tss <- ifelse(Tg_ATAC$TSS.enrichment > 3, 'High', 'Low')
TSSPlot(Tg_ATAC, group.by = 'high.tss') + NoLegend()

Tg_ATAC$nucleosome_group <- ifelse(Tg_ATAC$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
FragmentHistogram(object = Tg_ATAC, group.by = 'nucleosome_group')

VlnPlot(
  object = Tg_ATAC,
  features = c('pct_reads_in_peaks', 'peak_region_fragments',
               'TSS.enrichment', 'blacklist_ratio', 'nucleosome_signal'),
  pt.size = 0.1,
  ncol = 5
)


Tg_ATAC <- subset(
  x = Tg_ATAC,
  subset = peak_region_fragments > 200 &
    peak_region_fragments < 6000 &
    pct_reads_in_peaks > 40 &
    #blacklist_ratio < 0.05 &
    nucleosome_signal < 4 &
    TSS.enrichment > 2
)
Tg_ATAC



Tg_ATAC <- RunTFIDF(Tg_ATAC)
Tg_ATAC <- FindTopFeatures(Tg_ATAC, min.cutoff = 'q0')
Tg_ATAC <- RunSVD(Tg_ATAC)

DepthCor(Tg_ATAC)

## Must remove highly correlating components
Tg_ATAC <- RunUMAP(object = Tg_ATAC, reduction = 'lsi', dims = seq(1:30)[-c(1,3)])
Tg_ATAC <- FindNeighbors(object = Tg_ATAC, reduction = 'lsi', dims = seq(1:30)[-c(1,3)])
Tg_ATAC <- FindClusters(object = Tg_ATAC, verbose = FALSE, algorithm = 3)


DimPlot(object = Tg_ATAC, label = TRUE, reduction = 'umap') + NoLegend()

DefaultAssay(Tg_ATAC) <- "peaks"
gene.activities <- GeneActivity(Tg_ATAC, extend.upstream = 600,
                                extend.downstream = 200)


##### Merging GeneActivity with scRNA

## gene.activity matrix created with other approaches can be passed here
S.O.ATAC <- CreateSeuratObject(counts = gene.activities)
S.O.ATAC$orig.ident <- 'scATAC'

S.O.ATAC <- NormalizeData(
  object = S.O.ATAC,
  assay = 'RNA',
  normalization.method = 'LogNormalize',
  scale.factor = median(S.O.ATAC$nCount_RNA)
)

saveRDS(S.O.ATAC, '../Input/toxo_cdc/rds/ME49_59/S.O_ATAC_not_integrated_not_down_samples.rds')

DefaultAssay(S.O.ATAC) <- 'RNA'
S.O.ATAC <- FindVariableFeatures(S.O.ATAC, selection.method = "vst", nfeatures = 6000)
S.O.list <- list(RNA = S.O.ME49, ATAC = S.O.ATAC)
features <- SelectIntegrationFeatures(object.list = S.O.list, nfeatures = 6000)
reference_dataset <- 1
anchors <- FindIntegrationAnchors(object.list = S.O.list, 
                                  anchor.features = features, reference = reference_dataset)
S.O.integrated <- IntegrateData(anchorset = anchors)
# switch to integrated assay. Make sure to set to RNA for Differential Expression
DefaultAssay(S.O.integrated) <- "integrated"
S.O.integrated <- ScaleData(object = S.O.integrated, verbose = FALSE)
S.O.integrated <- RunPCA(S.O.integrated, features = VariableFeatures(object = S.O.integrated))
S.O.integrated <- FindVariableFeatures(S.O.integrated, nfeatures = 6000)
S.O.integrated <- FindNeighbors(S.O.integrated, dims = 1:10, reduction = 'pca')
S.O.integrated <- FindClusters(S.O.integrated, resolution = 0.2)
S.O.integrated <- RunUMAP(S.O.integrated, dims = 1:13)



## Transfer labels to scATAC
Idents(S.O.integrated) <- 'orig.ident'

atac_sub <- subset(S.O.integrated, ident = 'scATAC')
rna_sub <- subset(S.O.integrated, ident = 'intra')

anchors <- FindTransferAnchors(reference = rna_sub, query = atac_sub, dims = 1:30)
predictions <- TransferData(anchorset = anchors, refdata = rna_sub@meta.data$phase,dims = 1:30)
atac_sub <- AddMetaData(object = atac_sub, metadata = predictions)
atac_sub@meta.data$phase <- atac_sub@meta.data$predicted.id
Idents(atac_sub) <- 'phase'
DimPlot(atac_sub, reduction = 'pca')

ind1 <- S.O.integrated@meta.data$orig.ident == 'scATAC'
ind2 <- match(rownames(S.O.integrated@meta.data)[ind1], rownames(atac_sub@meta.data))
S.O.integrated@meta.data$phase[ind1] <- atac_sub@meta.data$phase[ind2]

ind <- S.O.integrated$orig.ident == 'intra'
S.O.integrated$orig.ident[ind] <- 'scRNA'
Idents(S.O.integrated) <- 'phase'
S.O.integrated$phase <- factor(S.O.integrated$phase, levels = c('G1.a', 'G1.b', 'S', 'M', 'C'))
S.O.integrated@meta.data$phase <- factor(S.O.integrated@meta.data$phase, levels = c('G1.a', 'G1.b', 'S', 'M', 'C'))
p <- DimPlot(S.O.integrated, reduction = "pca", 
             #group.by = "cell", 
             split.by = 'orig.ident',
             pt.size = 1,
             #shape.by='spp',
             label = T, label.size = 5) + NoLegend() + 
  theme(panel.spacing = unit(0.5, "lines")) + 
  theme(axis.text.x = element_text(face="bold", size=12, angle=0)) +
  theme(axis.text.y = element_text(face="bold", size=12, angle=0)) +
  theme(
    axis.title.x = element_text(size=14, face="bold"),
    axis.title.y = element_text(size=14, face="bold")
  )


plot(p)



saveRDS(S.O.integrated, '../Input/toxo_cdc/rds/ME49_59/S.O.intra_atac_integrated.rds')



### Differential peak expression
Tg_ATAC[['RNA']] <- CreateAssayObject(counts = gene.activities)

Tg_ATAC <- NormalizeData(
  object = Tg_ATAC,
  assay = 'RNA',
  normalization.method = 'LogNormalize',
  scale.factor = median(Tg_ATAC$nCount_RNA)
)

DefaultAssay(Tg_ATAC) <- 'RNA'

## For PCA
DefaultAssay(Tg_ATAC) <- 'RNA'
Tg_ATAC <- FindVariableFeatures(Tg_ATAC, selection.method = "vst", nfeatures = 2000)

all.genes <- rownames(Tg_ATAC)
Tg_ATAC <- ScaleData(Tg_ATAC, features = all.genes)
Tg_ATAC <- RunPCA(Tg_ATAC, features = VariableFeatures(object = Tg_ATAC))
Tg_ATAC <- FindNeighbors(Tg_ATAC, dims = 1:10, reduction = 'pca')
Tg_ATAC <- FindClusters(Tg_ATAC, resolution = 0.2)
Tg_ATAC <- RunTSNE(object = Tg_ATAC,features = VariableFeatures(object = Tg_ATAC) )
DimPlot(object = Tg_ATAC, reduction = "tsne", label = TRUE) + NoLegend()


## Coverage Browser

Tg_ATAC <- AddMetaData(Tg_ATAC, atac_sub@meta.data)
#levels(Tg_ATAC) <- c("CD4 Naive","CD4 Memory","CD8 Naive","CD8 Effector","DN T","NK CD56bright","NK CD56Dim","pre-B",'pro-B',"pDC","DC","CD14 Mono",'CD16 Mono')
Idents(Tg_ATAC) <- 'phase'

##Find Markers
DefaultAssay(Tg_ATAC) <- 'peaks'
da_peaks <- FindAllMarkers(
  object = Tg_ATAC,
  only.pos = T,
  min.pct = 0.2,
  test.use = 'LR',
  latent.vars = 'peak_region_fragments'
)

head(da_peaks)

plot1 <- VlnPlot(
  object = Tg_ATAC,
  features = rownames(da_peaks)[1],
  pt.size = 0.4,
  idents = c("G1.a","G1.b", 'S', 'M', 'C')
)

plot2 <- FeaturePlot(
  object = Tg_ATAC,
  features = rownames(da_peaks)[1],
  reduction = 'pca',
  pt.size = 0.4
)

plot1 | plot2

head(da_peaks)

top.da <- da_peaks %>% group_by(cluster) %>% slice_max(n = 1, order_by = avg_log2FC)
region <- StringToGRanges(regions = gsub('TGME49-', 'TGME49_', top.da$gene[4]),  sep = c("-", "-"))
#region <- StringToGRanges(regions = gsub('TGME49-', 'TGME49_', rownames(da_peaks)[5]),  sep = c("-", "-"))

xx <- findOverlaps(region, tx_trans)
tx_trans[xx@to]$gene_id

my.gene <- tx_trans[xx@to]$gene_id[1]

DefaultAssay(Tg_ATAC) <- 'RNA'

DefaultAssay(atac_sub) <- "RNA"
p1 <- FeaturePlot(
  object = atac_sub,
  features = gsub('_', '-', my.gene),
  pt.size = 0.4,
  max.cutoff = 'q0',
  ncol = 1,
  reduction = 'pca'
)

plot(p1)


saveRDS(Tg_ATAC, '../Input/toxo_cdc/rds/ME49_59/S.O_ATAC_peak.rds')

#### region_gene assignments
## Generating region/gene peak  counts
##Find Markers

## OR Extract the region turn into a bed file only and Use Bulk Peak-Gene Assignment here on region_gene_assignment

peak.regions <- rownames(Tg_ATAC@assays$peaks@data)

regions <- lapply(peak.regions, function(rr){
  region <- StringToGRanges(regions = gsub('TGME49-', 'TGME49_', rr),  sep = c("-", "-"))
  #xx <- findOverlaps(region, tx_trans, maxgap = 100)
  xx <- findOverlaps(region, tx_trans, ignore.strand=TRUE)
  tx_trans[xx@to]$gene_id
  L = list(region = region, gene_id = tx_trans[xx@to]$gene_id)
})

regions <- lapply(regions, function(rr){
  if(length(rr$gene_id) == 0){
    rr$gene_id <- NA
  }
  data.frame(rr$gene_id, as.data.frame(rr$region))
})

regions <- bind_rows(regions)
colnames(regions) <- c('GeneID', 'chr', 'strt', 'stp', 'width', 'strand')

region_gene_assignment <- regions %>% 
  transmute(chr = chr, strt = strt, stp = stp, width = width, strand = strand, GeneID = GeneID)

region_gene_assignment <- region_gene_assignment %>% na.omit()

gtf.file <- "../Input/toxo_genomics/genome/ToxoDB-59_TgondiiME49_filter.gtf"
gtf <- read.table(gtf.file, header = F, sep = '\t', quote = NULL)
gtf.filt <- gtf %>% dplyr::filter(grepl('TGME49*', V1))
## filter gtf for transcripts only 
gtf.filt.trn <- gtf.filt %>% filter(V3 == "transcript")
gtf.filt.trn$gene_name <- gsub("\\;.*", "", gsub("transcript_id ", "", gsub("-t.*", "", gtf.filt.trn$V9)))
gtf.filt.trn$gene_name <- gsub("\"", "", gtf.filt.trn$gene_name)

region_gene_assignment$strand <- gtf.filt.trn$V7[match(region_gene_assignment$GeneID, gtf.filt.trn$gene_name)]
saveRDS(region_gene_assignment, '../Input/toxo_cdc/rds/ME49_59/sc_atac_regions_gene_assignment.rds')


region_gene_assignment.bed <- region_gene_assignment %>% 
  transmute(chr = chr, strt = strt, stp = stp, GeneID = GeneID, width = width, strand = strand)


write.table(region_gene_assignment.bed, '../Input/toxo_cdc/rds/ME49_59/sc_atac_regions_gene_assignment_bed.bed', 
            quote = F, row.names = F, col.names = F, sep = '\t')



#### Pseudo-time
library(slingshot)
library(gam)
library(princurve)
library(parallel)
library(dtwclust)
library(doParallel)
library(tidyverse)
library(tidytext)
library(fda)
library(MyEllipsefit)
library(openxlsx)


num.cores <- detectCores(all.tests = FALSE, logical = TRUE)

## Fit a pseudo-time curve and align using sync data
S.O.integrated <- readRDS('../Input/toxo_cdc/rds/ME49_59/S.O.intra_atac_integrated.rds')
S.O.integrated@meta.data$Sample <- rownames(S.O.integrated@meta.data)
Idents(S.O.integrated) <- 'orig.ident'

atac_sub <- subset(S.O.integrated, ident = 'scATAC')
rna_sub <- subset(S.O.integrated, ident = 'scRNA')


## Pseudo-time analysis with SLingshot

fitTime <- function(S.O, method = 'rna', reverse.t = F){
  pc.tg <- getPCA(S.O)
  sds.data <- getPrinCurve(pc.tg)
  ind <- match(sds.data$Sample, pc.tg$Sample)
  sds.data$PC_1 <- pc.tg$PC_1[ind]
  sds.data$PC_2 <- pc.tg$PC_2[ind]
  sds.data$phase <- S.O@meta.data$phase[match(sds.data$Sample, rownames(S.O@meta.data))]
  sds.data$phase <- factor(sds.data$phase, levels = c('G1.a', 'G1.b', 'S', 'M', 'C'))
  
  Idents(S.O) <- 'phase'
  DimPlot(S.O, reduction = 'pca', label = T)
  
  
  pt <- sds.data$pt
  
  sds.data$pt <- 6 * ((as.numeric(pt) - min(as.numeric(pt)))/(max(as.numeric(pt)) - min(as.numeric(pt))))
  
  plot(sds.data$phase, sds.data$pt)
  
  if(reverse.t){
    sds.data$pt <- 6 - sds.data$pt
  }
  
  # Shift the time to start at G1.a
  tmp <- sds.data %>% dplyr::filter(phase == 'G1.a') %>% arrange(pt)
  lag.time <- tmp$pt[which.max((tmp$pt[2:length(tmp$pt)] - tmp$pt[1:(length(tmp$pt) - 1)])) + 1]
  lag.time <- quantile(tmp$pt, p = 0.205) ## excluse the ones overlapping with G1.b
  
  sds.data$pt <- (sds.data$pt - lag.time + 6) %% 6
  
  
  plot(sds.data$phase, sds.data$pt)
  
  
  ### NEW
  plot(sds.data$phase, sds.data$pt)
  outlier.c <- sds.data$Sample[(sds.data$phase == "C" & sds.data$pt < 1)]
  sds.data$phase[(sds.data$phase == "C" & sds.data$pt < 1)] <- 'G1.a'
  outlier.G1a <- sds.data$Sample[(sds.data$phase == "G1.a" & sds.data$pt > 4)]
  sds.data$phase[(sds.data$phase == "G1.a" & sds.data$pt > 4)] <- 'C'
  
  ####
  
  
  ind.G1.a <- which(sds.data$phase == 'G1.a')
  ind.G1.b <- which(sds.data$phase == 'G1.b')
  ind.S <- which(sds.data$phase == 'S')
  ind.M <- which(sds.data$phase == 'M')
  ind.C <- which(sds.data$phase == 'C')
  
  
  L <- wiskerPlot(S.O)
  
  par(mar = c(5, 5, 4, 4) + 0.1)
  plot(x = -35:10, y = -25:20, type = 'n', xlab = 'PC1', ylab = 'PC2',  lwd = 2, cex.lab = 1.5, cex.main = 2, cex.axis = 1.5)
  whiskers(as.matrix(L$pc[,c(1,2)]), L$fit$s, col = "gray")
  points(sds.data$PC_1, sds.data$PC_2, cex = 0.6, col = sds.data$phase, pch = 20)
  points(sds.data$sc1[sds.data$cell.ord],sds.data$sc2[sds.data$cell.ord], cex = 0.2, col = 'red')
  
  
  # Scale the pt based on know biology: Radke et. al 2000
  G <- c(0, 3) # 3h
  S <- c(3, 4.7) # 1.7h
  M <- c(4.7, 5) # ~20 min
  C <- c(5, 6) # 1h
  
  t1 <- (quantile(sds.data$pt[ind.G1.a], prob=0.75) + 
           quantile(sds.data$pt[ind.G1.b], prob=0.25))/2 ## Start of G1.b
  t2 <- (quantile(sds.data$pt[ind.G1.b], prob=0.75) + 
           quantile(sds.data$pt[ind.S], prob=0.25)) / 2 ## End of G1.b, Start of S
  t3 <- (quantile(sds.data$pt[ind.S], prob=0.75) + 
           quantile(sds.data$pt[ind.M], prob=0.25)) / 2## End of S, Start of M
  t4 <- (quantile(sds.data$pt[ind.M], prob=0.75) + 
           quantile(sds.data$pt[ind.C], prob=0.25)) / 2## End of M, Start of C
  
  t0 <- 0
  t5 <- 6
  
  # sds.data$pt.shift <- (sds.data$pt + t.shift) %% 6
  # sds.data$pt.shift <- 6 * ((sds.data$pt.shift - min(sds.data$pt.shift)) / 
  #                             (max(sds.data$pt.shift) - min(sds.data$pt.shift)))
  # plot(sds.data$phase, sds.data$pt.shift)
  
  slp.g <- (G[2] - G[1]) / (t2 - t0)
  inc.g  <- c(t0, G[1])
  
  slp.s <- (S[2] - S[1]) / (t3 - t2)
  inc.s  <- c(t2, S[1])
  
  slp.m <- (M[2] - M[1]) / (t4 - t3)
  inc.m  <- c(t3, M[1])
  
  slp.c <- (C[2] - C[1]) / (t5 - t4)
  inc.c  <- c(t4, C[1])
  
  sds.data <- sds.data %>%  
    mutate(pt.shifted.scaled = case_when(phase %in% c('G1.a', 'G1.b') ~ inc.g[2] + slp.g * (pt - inc.g[1]),
                                         phase == 'S' ~ inc.s[2] + slp.s * (pt - inc.s[1]),
                                         phase == 'M' ~ inc.m[2] + slp.m * (pt - inc.m[1]),
                                         phase == 'C' ~ inc.c[2] + slp.c * (pt - inc.c[1])))
  
  plot(sds.data$phase, sds.data$pt.shifted.scaled)
  
  ## Exclude outlier samples
  #q.ex <- quantile(sds.data$pt.shifted.scaled, p = 0.998)
  q.ex <- 6.5
  sds.data <- sds.data %>% dplyr::filter(pt.shifted.scaled <= q.ex)
  ## Rescale to [0, 6]
  sds.data$pt.shifted.scaled <- 6 * ((sds.data$pt.shifted.scaled - min(sds.data$pt.shifted.scaled))/
                                       (max(sds.data$pt.shifted.scaled) - min(sds.data$pt.shifted.scaled)))
  plot(sds.data$phase, sds.data$pt.shifted.scaled)
  # plot(sort(sds.data$pt.shifted.scaled))
  # 
  
  S.O.filt <- S.O
  S.O.filt@meta.data$Sample <- rownames(S.O.filt@meta.data)
  Idents(S.O.filt) <- 'Sample'
  S.O.filt <- subset(S.O.filt, idents = sds.data$Sample)
  
  
  genes.expr <- as.matrix(S.O.filt@assays$RNA@data)
  #genes.expr <- as.matrix(S.O@assays$RNA@data)
  
  
  genes.df <- data.frame(GeneID = rownames(genes.expr),
                         genes.expr) %>%
    pivot_longer(-c(GeneID), names_to = 'Sample', values_to = 'log2.expr')
  
  if(method == 'atac'){
    genes.df$Sample <- gsub('\\.', '-', genes.df$Sample)
  }
  
  genes.df <- inner_join(genes.df, sds.data, by = 'Sample')
  
  sds.data <- as.data.frame(sds.data)
  rownames(sds.data) <- sds.data$Sample
  
  ## Add the new clusters as meta-data
  S.O.filt <- AddMetaData(S.O.filt, sds.data)
  #S.O <- AddMetaData(S.O, sds.data)
  
  L <- list(sds.data = sds.data, genes.df = genes.df,
            S.O = S.O.filt)
  
  return(L)
}


DefaultAssay(atac_sub) <- 'RNA'
L.atac <- fitTime(atac_sub, 'atac', reverse.t = T)

DefaultAssay(rna_sub) <- 'RNA'
L.rna <- fitTime(rna_sub, 'rna', reverse.t = T)


saveRDS(L.rna$genes.df, '../Input/toxo_cdc/rds/ME49_59/sc_rna_genes_expr_pt.rds')
saveRDS(L.atac$genes.df, '../Input/toxo_cdc/rds/ME49_59/sc_atac_genes_expr_pt.rds')

saveRDS(L.rna$sds.data, '../Input/toxo_cdc/rds/ME49_59/sc_rna_sds_data.rds')
saveRDS(L.atac$sds.data, '../Input/toxo_cdc/rds/ME49_59/sc_atac_sds_data.rds')

saveRDS(L.rna$S.O, '../Input/toxo_cdc/rds/ME49_59/S.O_intra_lables_pt.rds')
saveRDS(L.atac$S.O, '../Input/toxo_cdc/rds/ME49_59/S.O_intra_atac_lables_pt.rds')



## Marker Analysis

## IDs
prod.desc  <- read.xlsx('../Input/toxo_genomics/genes/ProductDescription_GT1.xlsx')
TGGT1_ME49 <- read.xlsx('../Input/toxo_genomics/Orthologs/TGGT1_ME49 Orthologs.xlsx')
prod.desc <- left_join(prod.desc, TGGT1_ME49, by = c('GeneID' = 'TGGT1'))

## AP2s
AP2s <- read.xlsx('../Input/compScBdTgPb/genes/TF_Info_Updated_kz.xlsx')
AP2s <- left_join(AP2s, TGGT1_ME49, by = c('GeneID' = 'TGGT1'))
AP2s$TGME49 <- gsub('_', '-', AP2s$TGME49)

AP2s <- AP2s[grep("AP2", AP2s$Ap2Name),]


atac_sub <- readRDS('../Input/toxo_cdc/rds/ME49_59/S.O_intra_atac_lables_pt.rds')
rna_sub <- readRDS('../Input/toxo_cdc/rds/ME49_59/S.O_intra_lables_pt.rds')

# ## Differential gene expression
Idents(rna_sub) <- 'phase'
Intra.markers <- FindAllMarkers(object = rna_sub, only.pos = T, min.pct = 0)

Intra.markers$GeneID <- gsub('-', '_', Intra.markers$gene)
Intra.markers.top <- Intra.markers %>% group_by(cluster) %>% top_n(2, avg_log2FC)
FeaturePlot(object = rna_sub, 
            features = Intra.markers.top$gene, 
            cols = c("grey", "blue"), reduction = "pca")

Intra.markers.sig <- Intra.markers %>% dplyr::filter(avg_log2FC > log2(1.5) & p_val_adj < 0.05)

saveRDS(Intra.markers.sig, '../Input/toxo_cdc/rds/ME49_59/Intra_markers_sig.rds')

write.xlsx(Intra.markers.sig, '../Output/toxo_cdc/tabels/Intra_markers_sig_ME49_59.xlsx')
ss <- Intra.markers.sig %>% group_by(cluster) %>% summarise(num.DEG = n())
ss$cluster <- factor(ss$cluster, levels = c('G1.a', 'G1.b', 'S', 'M', 'C'))
p <- ggplot(data=ss, aes(x=cluster, y=num.DEG)) +
  geom_bar(stat="identity", fill="steelblue")+
  geom_text(aes(label=num.DEG), vjust=1.6, color="black", size=4, fontface="bold")+
  theme_minimal() + 
  ylab('DEGs') + xlab('') +
  theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 12, face="bold")) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 12, face="bold")) +
  theme(strip.background = element_rect(colour="black", fill="white",
                                        size=0.5, linetype="solid")) +
  theme(
    plot.title = element_text(size=14, face = "bold.italic", color = 'red'),
    axis.title.x = element_text(size=14, face="bold", hjust = 1),
    axis.title.y = element_text(size=14, face="bold")
  ) 

plot(p)

ggsave(filename="../Output/toxo_cdc/figures/tg_Intra_deg_numbers.pdf",
       plot=p,
       width = 8, height = 6,
       units = "in", # other options are "in", "cm", "mm"
       dpi = 300
)

print(ss)



## Differential accessibility analysis: ATAC
Idents(atac_sub) <- 'phase'
ATAC.markers <- FindAllMarkers(object = atac_sub, only.pos = TRUE, min.pct = 0)

ATAC.markers$GeneID <- gsub('-', '_', ATAC.markers$gene)
ATAC.markers.top <- ATAC.markers %>% group_by(cluster) %>% top_n(2, avg_log2FC)
FeaturePlot(object = atac_sub, 
            features = ATAC.markers.top$gene, 
            cols = c("grey", "blue"), reduction = "pca")

ATAC.markers.sig <- ATAC.markers %>% dplyr::filter(avg_log2FC > log2(1.5) & p_val_adj < 0.05)

saveRDS(ATAC.markers.sig, '../Input/toxo_cdc/rds/ME49_59/ATAC_markers_sig.rds')

ss <- ATAC.markers.sig %>% group_by(cluster) %>% summarise(num.DEG = n())
ss$cluster <- factor(ss$cluster, levels = c('G1.a', 'G1.b', 'S', 'M', 'C'))
p <- ggplot(data=ss, aes(x=cluster, y=num.DEG)) +
  geom_bar(stat="identity", fill="steelblue")+
  geom_text(aes(label=num.DEG), vjust=1.6, color="black", size=4, fontface="bold")+
  theme_minimal() + 
  ylab('DEGs') + xlab('') +
  theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 12, face="bold")) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 12, face="bold")) +
  theme(strip.background = element_rect(colour="black", fill="white",
                                        size=0.5, linetype="solid")) +
  theme(
    plot.title = element_text(size=14, face = "bold.italic", color = 'red'),
    axis.title.x = element_text(size=14, face="bold", hjust = 1),
    axis.title.y = element_text(size=14, face="bold")
  ) 

plot(p)

ggsave(filename="../Output/toxo_cdc/figures/tg_atac_deg_numbers.pdf",
       plot=p,
       width = 8, height = 6,
       units = "in", # other options are "in", "cm", "mm"
       dpi = 300
)

print(ss)


## DEGs (both Up & Down to identify depleated genes)
Idents(rna_sub) <- 'phase'
DEGs <- FindAllMarkers(object = rna_sub, only.pos = F, min.pct = 0)

DEGs$GeneID <- gsub('-', '_', DEGs$gene)
DEGs.sig <- DEGs %>% dplyr::filter(abs(avg_log2FC) > 0.5 & p_val_adj < 0.05)

DEGs.sig <- left_join(DEGs.sig, prod.desc, by = c('GeneID' = 'TGME49'))
saveRDS(DEGs.sig, '../Input/toxo_cdc/rds/Intra_DEGs_up_and_down_sig.rds')
write.xlsx(DEGs.sig, '../Output/toxo_cdc/tabels/Intra_DEGs_up_and_down_sig.xlsx')


### ATAC region-gene assignment phases

## Identify phase of the atach regions assigned to genes
Intra.markers.sig <- readRDS('../Input/toxo_cdc/rds/ME49_59/Intra_markers_sig.rds')
region_gene_assignment <- readRDS('../Input/toxo_cdc/rds/ME49_59/sc_atac_regions_gene_assignment.rds')
cell.cycle.markers <- Intra.markers.sig %>% transmute(GeneID = GeneID, phase = cluster)
region_cell_clycle_gene_assignment <- inner_join(region_gene_assignment, cell.cycle.markers, by = 'GeneID')

## IDs
prod.desc  <- read.xlsx('../Input/toxo_genomics/genes/ProductDescription_GT1.xlsx')
TGGT1_ME49 <- read.xlsx('../Input/toxo_genomics/Orthologs/TGGT1_ME49 Orthologs.xlsx')
TGGT1_ME49 <- inner_join(TGGT1_ME49, prod.desc, by = c('TGGT1' = 'GeneID'))
region_cell_clycle_gene_assignment <- left_join(region_cell_clycle_gene_assignment, TGGT1_ME49, by = c('GeneID' = 'TGME49'))
saveRDS(region_cell_clycle_gene_assignment, '../Input/toxo_cdc/rds/ME49_59/sc_atac_regions_cell_cycle_gene_assignment.rds')

## Fit Smoothing Splines
library(npreg)


sc.rna.genes.expr.pt <- readRDS('../Input/toxo_cdc/rds/ME49_59/sc_rna_genes_expr_pt.rds')
sc.atac.genes.expr.pt <- readRDS('../Input/toxo_cdc/rds/ME49_59/sc_atac_genes_expr_pt.rds')

#Intra.markers.sig <- readRDS('../Input/toxo_cdc/rds/Intra_markers_sig.rds')
#DEG.sig <- readRDS('../Input/toxo_cdc/rds/Intra_DEGs_up_and_down_sig.rds') ## both up & down


## Common genes between the data sets
comm.genes <- intersect(unique(sc.rna.genes.expr.pt$GeneID), unique(sc.atac.genes.expr.pt$GeneID))


do.filt <- F
## Fit smoothing splines and sample at regular time-points

AP2XII_8 <- 'TGME49-250800'
## Expression
sc.rna.spline.fits <- mclapply(1:length(comm.genes), function(i){
  tmp <- sc.rna.genes.expr.pt %>% dplyr::filter(GeneID == comm.genes[i]) %>%
    transmute(GeneID = GeneID, x = pt.shifted.scaled, y = log2.expr)
  y <- tmp$y
  t <- tmp$x
  
  #sc.rna.sp <- ss(t, y, periodic = T, lambda = 1e-4)
  w <- rep(1, length(y))
  w[which(y == 0)] <- 1/3
  sc.rna.sp <- smooth.spline(t, y, lambda = 0.1, w = w)
  
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

#saveRDS(sc.rna.spline.fits, '../Input/toxo_cdc/rds/sc_rna_spline_fits_cell_cycle_genes.rds')
saveRDS(sc.rna.spline.fits, '../Input/toxo_cdc/rds/ME49_59/sc_rna_spline_fits_all_genes.rds')

## Access
sc.atac.spline.fits <- mclapply(1:length(comm.genes), function(i){
  tmp <- sc.atac.genes.expr.pt %>% dplyr::filter(GeneID == comm.genes[i]) %>%
    transmute(GeneID = GeneID, x = pt.shifted.scaled, y = log2.expr)
  
  if(do.filt){
    ind <- smoothFilter(tmp$y)
    tmp2 <- tmp[-ind,]
    y <- tmp2$y
    t <- tmp2$x
  }else{
    y <- tmp$y
    t <- tmp$x
    
  }
  #sc.atac.sp <- ss(t, y, periodic = T, lambda = 1e-4)
  w <- rep(1, length(y))
  w[which(y == 0)] <- 1/3
  sc.atac.sp <- smooth.spline(t, y, lambda = 0.1, w = w)
  sc.atac.sp <- predict(sc.atac.sp, seq(0, 6, by = 1/3)) 
  
  mu <- data.frame(x = sc.atac.sp$x, y = sc.atac.sp$y)
  mu <- data.frame(GeneID = rep(tmp$GeneID[1], length(mu[,1])), x = mu[,1], y = mu[,2])
  return(mu)
}, mc.cores = num.cores)

sc.atac.spline.fits <- bind_rows(sc.atac.spline.fits)

#saveRDS(sc.atac.spline.fits, '../Input/toxo_cdc/rds/sc_atac_spline_fits_cell_cycle_genes.rds')
saveRDS(sc.atac.spline.fits, '../Input/toxo_cdc/rds/ME49_59/sc_atac_spline_fits_all_genes.rds')


### Plot trends
## IDs
prod.desc  <- read.xlsx('../Input/toxo_genomics/genes/ProductDescription_GT1.xlsx')
TGGT1_ME49 <- read.xlsx('../Input/toxo_genomics/Orthologs/TGGT1_ME49 Orthologs.xlsx')

## scDATA

rna_sub <- readRDS('../Input/toxo_cdc/rds/ME49_59/S.O_intra_lables_pt.rds')
atac_sub <- readRDS('../Input/toxo_cdc/rds/ME49_59/S.O_intra_atac_lables_pt.rds')


## Splines
sc.rna.spline.fits <- readRDS('../Input/toxo_cdc/rds/ME49_59/sc_rna_spline_fits_all_genes.rds')
sc.atac.spline.fits <- readRDS('../Input/toxo_cdc/rds/ME49_59/sc_atac_spline_fits_all_genes.rds')

## Turn the data into wide format (time by gene) and center & scale each gene
sc.rna.dtw.wide <- sc.rna.spline.fits %>% 
  pivot_wider(names_from = 'GeneID', values_from = 'y') %>% 
  mutate_at(vars(matches("TGME")), ~scale(., center = F, scale = F)) %>%
  as.data.frame()

sc.atac.dtw.wide <- sc.atac.spline.fits %>% 
  pivot_wider(names_from = 'GeneID', values_from = 'y') %>% 
  mutate_at(vars(matches("TGME")), ~scale(., center = F, scale = F)) %>%
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

AP2_degree <- read.xlsx('../Input/toxo_cdc/AP2_list/AP2_TFs_Network_Degree.xlsx')

AP2_degree <- left_join(AP2_degree, TGGT1_ME49, by = c('GeneID' = 'TGGT1'))
AP2_degree$Name <- gsub('AP2 domain transcription factor ', '', AP2_degree$ProductDescription)


for(i in 1:nrow(AP2_degree)){
  
  p <- plot_trends(gsub('_', '-', AP2_degree$TGME49[i]))
  f.n <- paste("../Output/toxo_cdc/figures/ME49_59/AP2_degree/", AP2_degree$Name[i], '.pdf', sep = '')
  ggsave(filename=f.n,
         plot=p,
         width = 6, height = 6,
         units = "in", # other options are "in", "cm", "mm"
         dpi = 300
  )
  
  
}


sig.AP2s <- readRDS('../Input/toxo_cdc/rds/sig_AP2s.rds')
TGGT1_ME49 <- read.xlsx('../Input/toxo_genomics/Orthologs/TGGT1_ME49 Orthologs.xlsx')
sig.AP2s <- sig.AP2s %>% transmute(GeneID = GeneID, Name = Ap2Name) %>% distinct()
sig.AP2s <- left_join(sig.AP2s, TGGT1_ME49, by = c('GeneID' = 'TGGT1'))


for(i in 1:nrow(sig.AP2s)){
  
  p <- plot_trends(gsub('_', '-', sig.AP2s$TGME49[i]))
  f.n <- paste("../Output/toxo_cdc/figures/sig_AP2s/", sig.AP2s$Name[i], '.pdf', sep = '')
  ggsave(filename=f.n,
         plot=p,
         width = 6, height = 6,
         units = "in", # other options are "in", "cm", "mm"
         dpi = 300
  )
  
  
}


AP2VIIa6 <- 'TGME49-237800'
p <- plot_trends(AP2VIIa6)
Idents(atac_sub) <- 'phase'
FeaturePlot(atac_sub, AP2VIIa6, reduction = 'pca', label = T)


AP2XII_2 = 'TGME49-289710'
p <- plot_trends(AP2XII_2)
Idents(atac_sub) <- 'phase'
FeaturePlot(rna_sub, AP2XII_2, reduction = 'pca', label = T)

#### Map transition points


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

## Read in the data.
## IDs
prod.desc  <- read.xlsx('../Input/toxo_genomics/genes/ProductDescription_GT1.xlsx')
TGGT1_ME49 <- read.xlsx('../Input/toxo_genomics/Orthologs/TGGT1_ME49 Orthologs.xlsx')
prod.desc <- left_join(prod.desc, TGGT1_ME49, by = c('GeneID' = 'TGGT1'))

marker.genes <- readRDS('../Input/toxo_cdc/rds/ME49_59/Intra_markers_sig.rds')
marker.genes.phase <- marker.genes %>% transmute(GeneID = gene, phase = cluster) %>% distinct()


sc.rna.spline.fits <- readRDS('../Input/toxo_cdc/rds/ME49_59/sc_rna_spline_fits_all_genes.rds')
sc.atac.spline.fits <- readRDS('../Input/toxo_cdc/rds/ME49_59/sc_atac_spline_fits_all_genes.rds')



rna_sub <- readRDS('../Input/toxo_cdc/rds/ME49_59/S.O_intra_lables_pt.rds')
atac_sub <- readRDS('../Input/toxo_cdc/rds/ME49_59/S.O_intra_atac_lables_pt.rds')

sc.rna.spline.fits.markers <- sc.rna.spline.fits %>% dplyr::filter(GeneID %in% marker.genes.phase$GeneID)
sc.atac.spline.fits.markers <- sc.atac.spline.fits %>% dplyr::filter(GeneID %in% marker.genes.phase$GeneID)

## Turn the data into wide format (time by gene) and center & scale each gene
sc.rna.dtw.wide <- sc.rna.spline.fits.markers %>% 
  pivot_wider(names_from = 'GeneID', values_from = 'y') %>% 
  mutate_at(vars(matches("TGME")), scale) %>%
  as.data.frame()

sc.atac.dtw.wide <- sc.atac.spline.fits.markers %>% 
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

sc.atac.mu.scale$GeneID <- factor(sc.atac.mu.scale$GeneID, 
                                  levels = unique(sc.atac.mu.scale$GeneID[order(-sc.atac.mu.scale$peak.ord)]))

## Filter to include markers only
sc.rna.mu.scale <- sc.rna.mu.scale %>% dplyr::filter(GeneID %in% marker.genes$gene)
sc.atac.mu.scale <- sc.atac.mu.scale %>% dplyr::filter(GeneID %in% marker.genes$gene)

sc.rna.mu.scale <- sc.rna.mu.scale %>% mutate(t.t = case_when(peak.ord >= 0 & peak.ord < 3 ~ 'G1',
                                                              peak.ord >= 3 & peak.ord < 4.7 ~ 'S',
                                                              peak.ord >= 4.7 & peak.ord < 5 ~ 'M',
                                                              peak.ord >= 5 ~ 'C'))

sc.rna.mu.scale$t.t <- factor(sc.rna.mu.scale$t.t, levels = c('G1', 'S', 'M', 'C'))

sc.atac.mu.scale <- sc.atac.mu.scale %>% mutate(t.t = case_when(peak.ord >= 0 & peak.ord < 3 ~ 'G1',
                                                                peak.ord >= 3 & peak.ord < 4.7 ~ 'S',
                                                                peak.ord >= 4.7 & peak.ord < 5 ~ 'M',
                                                                peak.ord >= 5 ~ 'C'))

sc.atac.mu.scale$t.t <- factor(sc.atac.mu.scale$t.t, levels = c('G1', 'S', 'M', 'C'))


tans.points <- function(mu.scale, lam = 0.01, prob = 0){
  peak.locs <- mu.scale %>% dplyr::select(GeneID, peak.ord) %>% distinct()
  spline.fit.peaks <- smooth.spline(x = 1:nrow(peak.locs), y = sort(peak.locs$peak.ord), lambda = lam)
  
  s.0 <- predict(spline.fit.peaks, seq(1, nrow(peak.locs), by = 0.1), deriv=0)
  s.1 <- predict(spline.fit.peaks, seq(1, nrow(peak.locs), by = 0.1), deriv=1)
  s.2 <- predict(spline.fit.peaks, seq(1, nrow(peak.locs), by = 0.1), deriv=2)
  
  spline.fit.peaks.smooth <- data.frame(x = seq(1, nrow(peak.locs), by = 0.1), 
                                        s0=s.0$y, s1=s.1$y, s2 = s.2$y)
  
  ## Get locations of peaks and valies
  max.loc <- getCurvePeakLoc(spline.fit.peaks.smooth$x, spline.fit.peaks.smooth$s1, prob = prob)
  min.loc <- getCurvePeakLoc(spline.fit.peaks.smooth$x, -spline.fit.peaks.smooth$s1, prob = prob)
  
  transition.points <- predict(spline.fit.peaks, sort(c(max.loc, min.loc)))
  
  L <- list(peak.locs = peak.locs, spline.fit.peaks.smooth = spline.fit.peaks.smooth, 
            transition.points = transition.points)
  
  return(L)
  
}

## ATAC transition points
L.trans.atac <- tans.points(sc.atac.mu.scale, lam = 0.1, prob = 0)

### NEW
L.trans.atac <- tans.points(sc.atac.mu.scale, lam = 0.01, prob = 0)
L.trans.atac$transition.points$y[L.trans.atac$transition.points$y < 0] <- 0
L.trans.atac$transition.points$y[L.trans.atac$transition.points$y > 6] <- 6
###

saveRDS(L.trans.atac, '../Input/toxo_cdc/rds/atac_based_transition_points.rds')

atac_peaks.dat <- L.trans.atac$spline.fit.peaks.smooth %>% 
  transmute(g = x, y = s0, yp = s1) %>% pivot_longer(-g, names_to = 'drivs', values_to = 'value')

### NEW
atac_peaks.dat <- L.trans.atac$spline.fit.peaks.smooth %>%  transmute(g = x, y = s0, yp = s1) %>% pivot_longer(-g, names_to = 'drivs', values_to = 'value')
atac_peaks.dat$value[atac_peaks.dat$value < 0] <- 0
atac_peaks.dat$value[atac_peaks.dat$value > 6] <- 6
###

atac_peaks.dat$drivs <- factor(atac_peaks.dat$drivs, levels = c('y', 'yp'))
p1  <- ggplot(atac_peaks.dat, aes(x= g,y=value)) +
  geom_path(aes(color = drivs),alpha = 0.8, size = 1)+ 
  theme_bw(base_size = 14) +
  geom_vline(xintercept=L.trans.atac$transition.points$x, linetype=2, color = 'orange', size = 0.6) + 
  ylab('peak time') + xlab('genes') +
  theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 12, face="bold")) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 12, face="bold")) +
  theme(strip.background = element_rect(colour="black", fill="white",
                                        size=0.5, linetype="solid")) +
  theme(strip.text = element_text(size = 14, face="bold", angle = 0)) + 
  
  
  facet_grid(drivs~., scales = 'free') +
  
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


## RNA transition points
L.trans.rna <- tans.points(sc.rna.mu.scale, lam = 0.1, prob = 0)



###### NEW
L.trans.rna <- tans.points(sc.rna.mu.scale, lam = 0.01, prob = 0)
L.trans.rna$transition.points$y[L.trans.rna$transition.points$y < 0] <- 0
L.trans.rna$transition.points$y[L.trans.rna$transition.points$y > 6] <- 6
#####



#saveRDS(L.trans.rnac, '../Input/toxo_cdc/rds/rna_based_transition_points.rds')

rna_peaks.dat <- L.trans.rna$spline.fit.peaks.smooth %>% 
  transmute(g = x, y = s0, yp = s1) %>% pivot_longer(-g, names_to = 'drivs', values_to = 'value')

####NEW
rna_peaks.dat <- L.trans.rna$spline.fit.peaks.smooth %>% transmute(g = x, y = s0, yp = s1) %>% pivot_longer(-g, names_to = 'drivs', values_to = 'value')
rna_peaks.dat$value[rna_peaks.dat$value < 0] <- 0
rna_peaks.dat$value[rna_peaks.dat$value > 6] <- 6

####

rna_peaks.dat$drivs <- factor(rna_peaks.dat$drivs, levels = c('y', 'yp'))
p2  <- ggplot(rna_peaks.dat, aes(x= g,y=value)) +
  geom_path(aes(color = drivs),alpha = 0.8, size = 1)+ 
  theme_bw(base_size = 14) +
  geom_vline(xintercept=L.trans.rna$transition.points$x, linetype=2, color = 'orange', size = 0.6) + 
  ylab('peak time') + xlab('genes') +
  theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 12, face="bold")) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 12, face="bold")) +
  theme(strip.background = element_rect(colour="black", fill="white",
                                        size=0.5, linetype="solid")) +
  theme(strip.text = element_text(size = 14, face="bold", angle = 0)) + 
  
  
  facet_grid(drivs~., scales = 'free') +
  
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

## Wrapping the curve around a circle (ignore)
x <- L.trans.atac$spline.fit.peaks.smooth$x
y <- L.trans.atac$spline.fit.peaks.smooth$s1
r <- sqrt(x^2 + y^2)
tt <- seq(0, 2*pi, length.out = length(y))
par(mfrow = c(1,1))
plot(cos(tt), y * sin(tt), type = 'l')

ggsave(filename="../Output/toxo_cdc/figures/atac_trnsition_curve.pdf",
       plot=p1,
       width = 8, height = 6,
       units = "in"
)


mapTransToPCA <- function(x, y, p){
  z = rep('T0', length(x))
  z[x >= 0 & x < y[1]] <- 'T1'
  for(i in 1:(length(y) - 1)){
    z[x >= y[i] & x < y[(i+1)]] <- paste('T', (i+1), sep = '')
  }
  
  z[x >= y[length(y)]] <- paste('T', length(y), sep = '')
  #z[p == 'C' & z == 'T1'] <- paste('T', length(y), sep = '')
  
  return(z)
  
}

t.atac.over.rna <- mapTransToPCA(rna_sub@meta.data$pt.shifted.scaled, L.trans.atac$transition.points$y, rna_sub@meta.data$phase)
rna_sub@meta.data$transition.atac <- factor(t.atac.over.rna, levels = sort(unique(t.atac.over.rna)))

t.rna.over.rna <- mapTransToPCA(rna_sub@meta.data$pt.shifted.scaled, L.trans.rna$transition.points$y, rna_sub@meta.data$phase)
rna_sub@meta.data$transition.rna <- factor(t.rna.over.rna, levels = sort(unique(t.rna.over.rna)))

t.atac.over.atac <- mapTransToPCA(atac_sub@meta.data$pt.shifted.scaled, L.trans.atac$transition.points$y, atac_sub@meta.data$phase)
atac_sub@meta.data$transition.atac <- factor(t.atac.over.atac, levels = sort(unique(t.atac.over.atac)))



saveRDS(rna_sub, '../Input/toxo_cdc/rds/S.O.intra_rna_atac_trnasition.rds')
saveRDS(atac_sub, '../Input/toxo_cdc/rds/S.O.intra_atac_atac_trnasition.rds')

Idents(rna_sub) <- 'phase'
p1 <- DimPlot(rna_sub, reduction = 'pca', label = T) + NoLegend()


Idents(rna_sub) <- 'transition.atac'
p2 <- DimPlot(rna_sub, reduction = 'pca', label = T) + NoLegend()




Idents(rna_sub) <- 'transition.rna'
p3 <- DimPlot(rna_sub, reduction = 'pca', label = T,  dims = c(2,3)) + NoLegend()

p1|p2|p3

p3

## Do a marker analysis on newly mapped transition points

# ## Differential gene expression
Idents(rna_sub) <- 'transition.atac'
transition.markers <- FindAllMarkers(object = rna_sub, only.pos = T, min.pct = 0.3)

transition.markers$GeneID <- gsub('-', '_', transition.markers$gene)
transition.markers.top <- transition.markers %>% group_by(cluster) %>% top_n(2, avg_log2FC)
FeaturePlot(object = rna_sub, 
            features = transition.markers.top$gene, 
            cols = c("grey", "blue"), reduction = "pca", label = T)

transition.markers.sig <- transition.markers %>% dplyr::filter(avg_log2FC > log2(1.5) & p_val_adj < 0.05)

saveRDS(transition.markers.sig, '../Input/toxo_cdc/rds/rna_transition_markers_sig.rds')

write.xlsx(transition.markers.sig, '../Output/toxo_cdc/tabels/rna_transition_markers_sig.xlsx')


## ATAC
# ## Differential gene expression
Idents(atac_sub) <- 'transition.atac'
transition.markers <- FindAllMarkers(object = atac_sub, only.pos = T, min.pct = 0)

transition.markers$GeneID <- gsub('-', '_', transition.markers$gene)
transition.markers.top <- transition.markers %>% group_by(cluster) %>% top_n(2, avg_log2FC)
FeaturePlot(object = atac_sub, 
            features = transition.markers.top$gene, 
            cols = c("grey", "blue"), reduction = "pca")

transition.markers.sig <- transition.markers %>% dplyr::filter(avg_log2FC > log2(1.5) & p_val_adj < 0.05)

saveRDS(transition.markers.sig , '../Input/toxo_cdc/rds/atac_transition_markers_sig.rds')

write.xlsx(transition.markers.sig, '../Output/toxo_cdc/tabels/atac_transition_markers_sig.xlsx')


ind.T2 <- rna_sub@meta.data$transition.rna == 'T2'
ind.T4 <- rna_sub@meta.data$transition.rna == 'T4'

plot(sort(rna_sub@meta.data$pt.shifted.scaled[ind.T2]), ylim = c(3, 6))
points(sort(rna_sub@meta.data$pt.shifted.scaled[ind.T4]), col = 'red')

max(rna_sub@meta.data$pt.shifted.scaled[ind.T4])
