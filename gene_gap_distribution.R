library(tidyverse)
library(seqinr)
## Read the annotations
bbig.gtf <- read.table('../Input/compScBabesia/genes/PiroplasmaDB-54_BbigeminaBOND.gtf', header = F, sep = '\t', quote = "")
bbov.gtf <- read.table('../Input/compScBabesia/genes/PiroplasmaDB-54_BbovisT2Bo.gtf', header = F, sep = '\t', quote = "")
bdiv.gtf <- read.table('../Input/compScBabesia/genes/PiroplasmaDB-54_Bdivergens1802A.gtf', header = F, sep = '\t', quote = "")
Pf.gtf <- read.table('../Input/compScBabesia/genes/PlasmoDB-52_Pfalciparum3D7.gtf', header = F, sep = '\t', quote = "")
Tg.gtf <- read.table('../Input/compScBabesia/genes/ToxoDB-52_TgondiiME49.gtf', header = F, sep = '\t', quote = "")


all.gtf <- list(bbig = bbig.gtf, bbov = bbov.gtf, bdiv = bdiv.gtf, Pf = Pf.gtf, Tg = Tg.gtf)

## calculate the gap between genes
spps <- names(all.gtf)
all.gtf <- lapply(1:length(all.gtf), function(i){
  spp <- spps[i]
  X <- all.gtf[[i]]
  
  ## Look at the entire gene only
  X <- X %>% dplyr::filter(V3 == 'transcript') %>% 
    transmute(Chr = V1, strt = V4, stp = V5, strand = V7) %>% arrange(Chr, strt) %>% unique()
  
  ## Strand; Sense (+), Antisense (-)
  X$Strand <- ifelse(X$strand == '+', 'S', 'A')
  total.chrs <- unique(X$Chr)
  
  ## On each chromosome, calculate the gap between genes
  ## Moving 5' to 3' Classify neighboring genes to:
  ## SS (both genes are on S strand)
  ## SA (first gene is on S, and second on A)
  ## AS (first gene is on A, and second on S)
  ## SS (both genes are on A strand)
  gap.class <- lapply(total.chrs, function(x){
    tmp <- X %>% dplyr::filter(Chr == x)
    nn <- nrow(tmp)
    if(nn > 1){
      ## ignore the last gene (set to NA)
      class <- c(paste(tmp$Strand[1:(nn-1)], tmp$Strand[2:nn], sep = ''), 'NA')
      gap <- c(tmp$strt[2:nn] - tmp$stp[1:(nn-1)], 'NA')
    }else{ ## only one gene on the chromosome
      class <- 'NA'
      gap <- 'NA'
    }
    return(list(class = class, gap = gap))
  })
 
  ## Extract the class and gap per chromosome and concatinate
  X$class <- unlist(lapply(gap.class, `[[`, 1))
  X$gap <- unlist(lapply(gap.class, `[[`, 2))

  X <- X %>% dplyr::filter(gap != 'NA')
  X$gap <- as.numeric(X$gap)
  X$gap[X$gap < 0] <- 0 ## overlapping genes
  X$spp <- spp 
  
  return(X)
})

names(all.gtf) <- spps

xx <- all.gtf$Pf
head(xx)

all.gtf <- bind_rows(all.gtf)

all.gtf$spp <- factor(all.gtf$spp, levels = spps)
all.gtf$class <- factor(all.gtf$class, levels = c('SS', 'SA', 'AS', 'AA'))


p <- ggplot(all.gtf, aes(x=gap, color = spp)) +
  geom_density(size = 1) + theme_bw() + 
  #facet_grid(class~.) + 
  theme(
    plot.title = element_text(size=14, face = "bold.italic", color = 'red'),
    axis.title.x = element_text(size=14, face="bold"),
    axis.title.y = element_text(size=14, face="bold")
  )+ expand_limits(x=c(0,12))
plot(p)



gap.stats <- all.gtf %>% group_by(spp) %>% 
  summarise(mean.gap = mean(gap), median.gap = median(gap), promoter.length = 1 * median(gap)/2)

print(gap.stats)


gap.stats.no.pf <- gap.stats %>% dplyr::filter(spp != 'Pf')
p <- ggplot(gap.stats.no.pf, aes(x=spp, y = median.gap, fill = spp)) +
  geom_bar(stat="identity")+
  geom_text(aes(label=round(median.gap)), vjust=1.6, color="black", size=5, fontface="bold")+
  theme_minimal() + 
  ylab('Median gap') + xlab('') +
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

ggsave(filename="../Output/compScBabsSpecies/revisions/figs/median_gap.pdf",
       plot=p,
       width = 6, height = 6,
       units = "in", # other options are "in", "cm", "mm"
       dpi = 300
)


## strand specific
gap.stats.ss <- all.gtf %>% group_by(spp, class) %>% 
  summarise(mean.gap = mean(gap), median.gap = median(gap), promoter.length = median(gap)/2)

print(gap.stats.ss)

gene.class.totals <- all.gtf %>% group_by(spp, class) %>% 
  summarise(total.genes = n()) %>% group_by(spp) %>%
  mutate(percent.genes = total.genes/sum(total.genes))

print(gene.class.totals)

write.xlsx(gap.stats.ss, '../Output/compScBabsSpecies/tables/intergenic_stats.xlsx')
write.xlsx(gene.class.totals, '../Output/compScBabsSpecies/tables/gene_class_totals.xlsx')


## Gap distribution in B. Divergence
bdiv.gtf <- all.gtf %>% dplyr::filter(spp == 'bdiv')
p <- ggplot(bdiv.gtf, aes(x=gap)) + 
  #geom_density(size = 2, color = 'blue') +  
  #geom_histogram(aes(y = ..density..), alpha = 0.4, fill = "red", color = 'red', bins = 10) + theme_bw() + xlim(0, 100) + 
  geom_histogram(alpha = 0.4, fill = "red", color = 'red', bins = 30) + theme_bw() + xlim(0, 500) + 
  #facet_grid(class~.) + 
  theme(
    plot.title = element_text(size=14, face = "bold.italic", color = 'red'),
    axis.title.x = element_text(size=14, face="bold"),
    axis.title.y = element_text(size=14, face="bold")
  )+ expand_limits(x=c(0,12))
plot(p)

ggsave(filename="../Output/compScBabsSpecies/revisions/figs/bdiv_gap_dist_500.pdf",
       plot=p,
       width = 6, height = 6,
       units = "in", # other options are "in", "cm", "mm"
       dpi = 300
)

