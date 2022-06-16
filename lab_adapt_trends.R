library(ggplot2)
library(ggfortify)
library(jcolors)
library(gridExtra)
library(grid)
library(tidyverse)
library(openxlsx)
library(tidytext)
library(ggrepel)
library(geomtextpath)


## Read Mean expression and sd accross passages for extra and intracellular
toxo.tab.means.sd.qval <- read.xlsx('../Input/toxo_cdc/tables/lab_adapt_mean_sd.xlsx')

## Read Gene Info
prod.desc <- read.xlsx('../Input/toxo_genomics/genes/ProductDescription_GT1.xlsx')


## Get the AP2 trends
AP2s <- prod.desc %>% dplyr::filter(grepl('AP2', ProductDescription))

## Simplify the names
AP2s$ProductDescription <- gsub("PAP2 superfamily protein", "",gsub("AP2 domain-containing protein", "", 
                                                     gsub("AP2 domain transcription factor ", "", AP2s$ProductDescription)))
AP2s <- AP2s %>% dplyr::filter(AP2 != "")


## Specific AP2s
## AP2s associated with Bradyzoites
activators <- c('AP2XII-6', 'AP2XI-4','AP2IV-3', 'AP2Ib-1')
repressors <- c('AP2IX-9','AP2IX-4','AP2IV-4')
others <- c('AP2XI-5', 'AP2X-5', 'AP2X-7')
stress_induced <- c('AP2XI-5', 'AP2X-5','AP2IX-9', 'AP2IX-4', 'AP2IV-4', 'AP2IV-3', 'AP2XI-4')
my.AP2s <- c(others, activators, repressors)

my.AP2s <- AP2s %>% dplyr::filter(ProductDescription %in% my.AP2s)

## Get the mean and sd of AP2s from the expression table
Ap2.tab.means.sd.qval <- right_join(toxo.tab.means.sd.qval, my.AP2s, by = c('GeneName' = 'GeneID'))

## Turn into factors for ordered plots
Ap2.tab.means.sd.qval$GeneName <- factor(Ap2.tab.means.sd.qval$GeneName,
                                         levels = unique(Ap2.tab.means.sd.qval$GeneName))
Ap2.tab.means.sd.qval$Passage <- factor(Ap2.tab.means.sd.qval$Passage,
                                        levels = unique(Ap2.tab.means.sd.qval$Passage))
Ap2.tab.means.sd.qval$cond <- factor(Ap2.tab.means.sd.qval$cond,
                                     levels = unique(Ap2.tab.means.sd.qval$cond))
Ap2.tab.means.sd.qval$ProductDescription <- factor(Ap2.tab.means.sd.qval$ProductDescription,
                                    levels = unique(Ap2.tab.means.sd.qval$ProductDescription))

Ap2.tab.means.sd.qval <- Ap2.tab.means.sd.qval %>%
  mutate(sig =ifelse(qval < 0.0001, '***', ifelse(qval < 0.001, '**', ifelse(qval < 0.01, '*', ''))))
Ap2.tab.means.sd.qval$sig[Ap2.tab.means.sd.qval$cond == 'intra'] = ''
p1<- ggplot(Ap2.tab.means.sd.qval, aes(x = Passage, y = mean, group = cond, color = cond)) + 
  geom_line(size=1) + 
  geom_point(size = 1) + 
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.1) +  
  geom_text(aes(x = Passage, y = 140, label=sig), color = 'black') + 
  theme_bw() + 
  facet_wrap(ProductDescription~., scales='free') + 
  theme(strip.background = element_rect(colour="black", fill="white", 
                                        size=1, linetype="solid"))

plot(p1)
ggsave(filename="~/work/ToxoplasmaGondii/ToxoR/output/AP2_trends.pdf", plot=p1,
       width = 8, height = 8, 
       units = "in", # other options are "in", "cm", "mm" 
       dpi = 300
)

