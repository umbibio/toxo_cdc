library(openxlsx)
library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(gridExtra)
library(fdrtool)


GO.terms <- read.xlsx('../Input/toxo_genomics/gene_function/ME49_toxo_db_go_terms.xlsx')



GO.terms.by.id <- GO.terms %>% group_by(GO.ID, GO.name, category, type)  %>% 
  summarise(genes = list(as.character(GeneID)), total = n())


gene.clusters <- read.xlsx('../Output/toxo_cdc/tabels/rna_gene_atac_clusters_phase.xlsx')
gene.clusters <- gene.clusters %>% dplyr::select(GeneID, transition.atac) %>%
  distinct() %>% group_by(transition.atac) %>%
  summarise(genes = list(gsub('-', '_', as.character(GeneID))), total = n())


marker.genes <- readRDS('../Input/toxo_cdc/rds/Intra_markers_sig.rds')
colnames(marker.genes)[7] <- 'genes'


marker.genes <- marker.genes  %>% group_by(cluster) %>%
  summarise(genes = list(gsub('-', '_', as.character(GeneID))), total = n())


## Total number of genes in toxo
all.genes.marker <- union(unique(unlist(GO.terms.by.id$genes)), unique(unlist(marker.genes$genes)))
nn <- length(all.genes.marker)

all.genes.cluster <- union(unique(unlist(GO.terms.by.id$genes)), unique(unlist(gene.clusters$genes)))
mm <- length(all.genes.cluster)



GO.terms.by.id$all <- nn
gene.clusters$all <- nn
marker.genes$all <-nn
XX <- full_join(GO.terms.by.id, gene.clusters, by = 'all')
XX <- XX %>% rowwise() %>% mutate(overlap = length(intersect(c(genes.x), c(genes.y))),
                                  overlap.genes = list(intersect(c(genes.x), c(genes.y))))
XX <- XX %>% rowwise() %>% 
  mutate(pvalue = fisher.test(matrix(c(overlap, total.x - overlap, total.y - overlap, 
                                       all - (total.x + total.y - overlap) ), byrow = T, ncol = 2, nrow = 2),
                              alternative = "greater")$p.value)

XX$qval <- fdrtool(XX$pvalue, statistic="pvalue",plot=F, color.figure=F, verbose=F, cutoff.method="fndr")$qval

saveRDS(XX, '../Input/toxo_cdc/rds/enrichment_cell_cycle_markers.rds')

XX.sig <- XX %>% dplyr::filter(pvalue <= 0.01) %>% arrange(pvalue)

p1 <- ggplot(XX.sig, aes(x = transition.atac, y = GO.name)) + 
  geom_point(aes(color = pvalue, size = -log(pvalue))) +
  theme_bw(base_size = 12) +
  scale_colour_gradient(limits=c(0, 0.01), low="red") +
  ylab(NULL) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + facet_grid(category~., scales='free') + 
  theme(strip.background = element_rect(colour="black", fill="white", 
                                        size=1, linetype="solid"))

plot(p1)

