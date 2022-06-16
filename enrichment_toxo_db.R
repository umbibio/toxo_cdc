library(openxlsx)
library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(gridExtra)
library(fdrtool)



## ATAC
in.dir <- '../Input/toxo_cdc/GO/atac_transition/'
all.files <- list.files(in.dir)

all.clust.items <- list()
for(f in all.files){
  nn <- gsub('\\.tsv', '', f)
  cluster <- strsplit(nn, split = '_')[[1]][1]
  GF <- strsplit(nn, split = '_')[[1]][4]
  tmp <- read_tsv(paste(in.dir, f, sep = ''))
  tmp$GF <- GF
  tmp$cluster <- cluster
  all.clust.items <- c(all.clust.items, list(tmp))
}


all.clust.items <- do.call(rbind, all.clust.items)

saveRDS(all.clust.items, '../Input/toxo_cdc/rds/GO_toxo_db_atact_transition_points.rds')

filtered.Go <- all.clust.items %>% arrange(cluster, Benjamini) %>% distinct() %>%
  group_by(cluster) %>% mutate(rank = row_number()) %>%
  dplyr::filter(Benjamini < 0.1 & rank < 30) %>% 
  arrange(cluster, Benjamini) 

# filtered.Go$phase <- gsub("Clust0" , "G1",
#                           gsub("Clust1" , "C", 
#                                gsub("Clust2", "M" ,
#                                     gsub("Clust3", "S", filtered.Go$phase))))
# 
filtered.Go$cluster <- factor(filtered.Go$cluster, levels = sort(unique(filtered.Go$cluster)))
filtered.Go$ID <- factor(filtered.Go$ID, level=unique(filtered.Go$ID))
filtered.Go$Name <- factor(filtered.Go$Name, level=unique(filtered.Go$Name))

write.xlsx(filtered.Go, '../Output/toxo_cdc/tabels/GO_toxo_db_atact_transition_points_sig.xlsx')

## Category of contrasts
p <- ggplot(filtered.Go, aes(x = cluster, y = Name)) + 
  geom_point(aes(colour = "red", size = -log(Benjamini))) +
  theme_bw(base_size = 14) +
  #scale_colour_gradient(limits=c(0, 0.01), low="red") +
  ylab(NULL) + xlab(NULL) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10, face="bold")) + 
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 10, face="bold")) +
  theme(legend.position="none") +
  theme(strip.background = element_rect(colour="black", fill="white", 
                                        size=0.5, linetype="solid")) + guides(color = FALSE)+
  ggtitle(strsplit(in.dir,split = "/")[[1]][5]) +
  theme(plot.title = element_text(size = 10))

plot(p)


ggsave(filename="../Output/compScBabsSpecies/figs/GO_Bdiv_human_TGGT1_Cell_Cycle_phases.png", 
       plot=p,
       width = 10, height = 14, 
       units = "in", # other options are "in", "cm", "mm" 
       dpi = 300
)


## Cell Cycle Phases
in.dir <- '../Input/toxo_cdc/GO/phases/'
all.files <- list.files(in.dir)

all.clust.items <- list()
for(f in all.files){
  nn <- gsub('\\.tsv', '', f)
  cluster <- strsplit(nn, split = '_')[[1]][1]
  GF <- strsplit(nn, split = '_')[[1]][3]
  tmp <- read_tsv(paste(in.dir, f, sep = ''))
  tmp$GF <- GF
  tmp$cluster <- cluster
  all.clust.items <- c(all.clust.items, list(tmp))
}


all.clust.items <- do.call(rbind, all.clust.items)

saveRDS(all.clust.items, '../Input/toxo_cdc/rds/GO_toxo_db_phases.rds')

filtered.Go <- all.clust.items %>% arrange(cluster, Benjamini) %>% distinct() %>%
  group_by(cluster) %>% mutate(rank = row_number()) %>%
  dplyr::filter(Benjamini < 0.1 & rank < 30) %>% 
  arrange(cluster, Benjamini) 

# filtered.Go$phase <- gsub("Clust0" , "G1",
#                           gsub("Clust1" , "C", 
#                                gsub("Clust2", "M" ,
#                                     gsub("Clust3", "S", filtered.Go$phase))))
# 
filtered.Go$cluster <- factor(filtered.Go$cluster, levels = c('G1.a', 'G1.b', 'S', 'M', 'C'))
filtered.Go$ID <- factor(filtered.Go$ID, level=unique(filtered.Go$ID))
filtered.Go$Name <- factor(filtered.Go$Name, level=unique(filtered.Go$Name))

write.xlsx(filtered.Go, '../Output/toxo_cdc/tabels/GO_toxo_db_phases_sig.xlsx')

## Category of contrasts
p <- ggplot(filtered.Go, aes(x = cluster, y = Name)) + 
  geom_point(aes(colour = "red", size = -log(Benjamini))) +
  theme_bw(base_size = 14) +
  #scale_colour_gradient(limits=c(0, 0.01), low="red") +
  ylab(NULL) + xlab(NULL) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10, face="bold")) + 
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 10, face="bold")) +
  theme(legend.position="none") +
  theme(strip.background = element_rect(colour="black", fill="white", 
                                        size=0.5, linetype="solid")) + guides(color = FALSE)+
  ggtitle(strsplit(in.dir,split = "/")[[1]][5]) +
  theme(plot.title = element_text(size = 10))

plot(p)


ggsave(filename="../Output/compScBabsSpecies/figs/GO_Bdiv_human_TGGT1_Cell_Cycle_phases.png", 
       plot=p,
       width = 10, height = 14, 
       units = "in", # other options are "in", "cm", "mm" 
       dpi = 300
)
