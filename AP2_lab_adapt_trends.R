###
toxo.tab <- read.xlsx('~/work/ToxoplasmaGondii/Final.Modified.Summary_KZ.xlsx')
colnames(toxo.tab) <- gsub("fold.change", 'fold_change', gsub("qValue", 'q-value', colnames(toxo.tab)))

toxo.tab.B2 <- toxo.tab %>% dplyr::select(-contains("B4")) %>% dplyr::select(-contains("freshlyse"))  
colnames(toxo.tab.B2) <- gsub('B\\.', '',gsub('over', 'vs',gsub('6hes', 'extra', 
                                                                gsub('B2\\.', '',colnames(toxo.tab.B2)))))

toxo.tab.extra.expr.sd <- toxo.tab.B2 %>% 
  dplyr::select(matches(c('.*_FPKM|.*_sd'))) %>% 
  dplyr::select(-contains("intra")) %>% 
  dplyr::select(-contains("freshlyse"))  %>%  dplyr::select(-contains("_quant"))


toxo.tab.intra.expr.sd <- toxo.tab.B2 %>% 
  dplyr::select(matches(c('.*_FPKM|.*_sd'))) %>% 
  dplyr::select(-contains("extra")) %>% 
  dplyr::select(-contains("freshlyse"))  %>%  dplyr::select(-contains("_quant"))

toxo.tab.qval <- toxo.tab.B2 %>% 
  dplyr::select(matches(c('.*extra.*intra.*q-value|.*intra.*extra.*q-value' )))

## Replace RH with P1000 for sorting
colnames(toxo.tab.extra.expr.sd) <- gsub('RH', 'P1000', colnames(toxo.tab.extra.expr.sd))
colnames(toxo.tab.intra.expr.sd) <- gsub('RH', 'P1000', colnames(toxo.tab.intra.expr.sd))
colnames(toxo.tab.qval) <- gsub('RH', 'P1000', colnames(toxo.tab.qval))
colnames(toxo.tab.qval) <- gsub('\\.intra.*', '', gsub('\\.extra.*', '', colnames(toxo.tab.qval)))

nn.extra <- gsub("\\.extra_FPKM", "", colnames(toxo.tab.extra.expr.sd))
nn.extra.mean <- c(1:length(nn.extra))[-grep('_sd', nn.extra)]
nn.extra.sd <- c(1:length(nn.extra))[grep('_sd', nn.extra)]
sort.idx <- c(nn.extra.mean[sort(as.numeric(gsub("P", "", nn.extra[nn.extra.mean])), index.return = T)$ix], 
              nn.extra.sd[sort(as.numeric(gsub("_sd", "", gsub("P", "", nn.extra[nn.extra.sd]))), index.return = T)$ix])

toxo.tab.extra.expr.sd <- toxo.tab.extra.expr.sd[, sort.idx]
colnames(toxo.tab.extra.expr.sd) <- nn.extra[sort.idx]

nn.intra <- gsub("\\.intra_FPKM", "", colnames(toxo.tab.intra.expr.sd))
nn.intra.mean <- c(1:length(nn.intra))[-grep('_sd', nn.intra)]
nn.intra.sd <- c(1:length(nn.intra))[grep('_sd', nn.intra)]
sort.idx <- c(nn.intra.mean[sort(as.numeric(gsub("P", "", nn.intra[nn.intra.mean])), index.return = T)$ix], 
              nn.intra.sd[sort(as.numeric(gsub("_sd", "", gsub("P", "", nn.intra[nn.intra.sd]))), index.return = T)$ix])

toxo.tab.intra.expr.sd <- toxo.tab.intra.expr.sd[, sort.idx]
colnames(toxo.tab.intra.expr.sd) <- nn.intra[sort.idx]

nn.qval <- colnames(toxo.tab.qval)
sort.idx <- sort(as.numeric(gsub('P', '',nn.qval)), index.return = T)$ix
toxo.tab.qval <- toxo.tab.qval[, sort.idx]
colnames(toxo.tab.qval) <- nn.qval[sort.idx]

## There is a slight shift in time. Just take column names from extra
colnames(toxo.tab.extra.expr.sd) == colnames(toxo.tab.intra.expr.sd)

colnames(toxo.tab.intra.expr.sd) = colnames(toxo.tab.extra.expr.sd)

toxo.tab.extra.expr.sd <- data.frame(GeneName = as.character(toxo.tab$GeneName), toxo.tab.extra.expr.sd)
toxo.tab.intra.expr.sd <- data.frame(GeneName = as.character(toxo.tab$GeneName), toxo.tab.intra.expr.sd)
toxo.tab.qval <- data.frame(GeneName = as.character(toxo.tab$GeneName), toxo.tab.qval)

colnames(toxo.tab.extra.expr.sd) <- gsub('P1000', 'RH', colnames(toxo.tab.extra.expr.sd))
colnames(toxo.tab.intra.expr.sd) <- gsub('P1000', 'RH', colnames(toxo.tab.intra.expr.sd))
colnames(toxo.tab.qval) <- gsub('P1000', 'RH', colnames(toxo.tab.qval))

toxo.tab.extra.means <- toxo.tab.extra.expr.sd %>% dplyr::select(-contains("_sd")) %>% 
  gather(key = Passage, value = extra, P7:RH)
toxo.tab.extra.sds <- toxo.tab.extra.expr.sd %>% dplyr::select(matches("GeneName|.*_sd")) %>% 
  gather(key = Passage, value = extra, P7_sd:RH_sd)
toxo.tab.extra.sds$extra <- toxo.tab.extra.sds$extra * sqrt(3) / 1.96


toxo.tab.intra.means <- toxo.tab.intra.expr.sd %>% dplyr::select(-contains("_sd")) %>% 
  gather(key = Passage, value = intra, P7:RH)
toxo.tab.intra.sds <- toxo.tab.intra.expr.sd %>% dplyr::select(matches("GeneName|.*_sd")) %>% 
  gather(key = Passage, value = intra, P7_sd:RH_sd)
toxo.tab.intra.sds$intra <- toxo.tab.intra.sds$intra * sqrt(3) / 1.96
toxo.tab.qval <- toxo.tab.qval %>% 
  gather(key = Passage, value = qval, P7:RH)

toxo.tab.means <- full_join(toxo.tab.extra.means, toxo.tab.intra.means, by = c('GeneName', 'Passage'))
toxo.tab.sds <- full_join(toxo.tab.extra.sds, toxo.tab.intra.sds, by = c('GeneName', 'Passage'))
toxo.tab.sds$Passage <- gsub('_sd', '', toxo.tab.sds$Passage)

toxo.tab.means <- toxo.tab.means %>% gather(key = cond, value=mean,extra:intra)
toxo.tab.sds <- toxo.tab.sds %>% gather(key = cond, value=sd,extra:intra)
toxo.tab.means.sd <- full_join(toxo.tab.means, toxo.tab.sds, by = c('GeneName', 'Passage', 'cond'))
toxo.tab.means.sd.qval <- left_join(toxo.tab.means.sd, toxo.tab.qval, by = c('GeneName', 'Passage'))
toxo.tab.means.sd.qval$GeneName <- as.character(toxo.tab.means.sd.qval$GeneName)


AP2s.ind <- grep('AP2', toxo.tab.B2$Product.Description)
AP2s <- data.frame(GeneName = toxo.tab.B2$GeneName[AP2s.ind], AP2 = toxo.tab.B2$Product.Description[AP2s.ind],
                   stringsAsFactors = F)
AP2s$AP2 <- gsub("PAP2 superfamily protein", "",gsub("AP2 domain-containing protein", "", 
                                                     gsub("AP2 domain transcription factor ", "", AP2s$AP2)))

AP2s <- AP2s %>% dplyr::filter(AP2 != "")
my.AP2s <- c('AP2XI-5', 'AP2X-5','AP2IX-9', 'AP2IX-4', 'AP2IV-4', 'AP2IV-3', 'AP2XI-4')
my.AP2s <- c('AP2Ib-1', 'AP2IV-3', 'AP2VI-3', 'AP2VIIa-1', 'AP2VIII-4', 'AP2IX-9', 
             'AP2XI-5', 'AP2X-5', 'AP2IX-4','AP2IV-4', 'AP2XI-4')
my.AP2s <- c('AP2Ib-1', 'AP2IV-3', 'AP2IX-9', 
             'AP2XI-5', 'AP2X-5', 'AP2IX-4','AP2IV-4', 'AP2XI-4')

## AP2s associated with Bradyzoites
activators <- c('AP2XII-6', 'AP2XI-4','AP2IV-3', 'AP2Ib-1')
repressors <- c('AP2IX-9','AP2IX-4','AP2IV-4')
others <- c('AP2XI-5', 'AP2X-5', 'AP2X-7')
stress_induced <- c('AP2XI-5', 'AP2X-5','AP2IX-9', 'AP2IX-4', 'AP2IV-4', 'AP2IV-3', 'AP2XI-4')
#my.AP2s <- c(others, activators, repressors)

my.AP2s <- AP2s %>% dplyr::filter(AP2 %in% my.AP2s)

Ap2.tab.means.sd.qval <- right_join(toxo.tab.means.sd.qval, my.AP2s, by = 'GeneName')

Ap2.tab.means.sd.qval$GeneName <- factor(Ap2.tab.means.sd.qval$GeneName,
                                         levels = unique(Ap2.tab.means.sd.qval$GeneName))

Ap2.tab.means.sd.qval$Passage <- factor(Ap2.tab.means.sd.qval$Passage,
                                         levels = unique(Ap2.tab.means.sd.qval$Passage))

Ap2.tab.means.sd.qval$cond <- factor(Ap2.tab.means.sd.qval$cond,
                                        levels = unique(Ap2.tab.means.sd.qval$cond))

Ap2.tab.means.sd.qval$AP2 <- factor(Ap2.tab.means.sd.qval$AP2,
                                     levels = unique(Ap2.tab.means.sd.qval$AP2))
Ap2.tab.means.sd.qval <- Ap2.tab.means.sd.qval %>%
  mutate(sig =ifelse(qval < 0.0001, '***', ifelse(qval < 0.001, '**', ifelse(qval < 0.01, '*', ''))))
Ap2.tab.means.sd.qval$sig[Ap2.tab.means.sd.qval$cond == 'intra'] = ''
p1<- ggplot(xx, aes(x = Passage, y = mean, group = cond, color = cond)) + 
  geom_line(size=1) + 
  geom_point(size = 1) + 
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.1) +  
  geom_text(aes(x = Passage, y = 140, label=sig), color = 'black') + 
  theme_bw() + 
  facet_wrap(AP2~., scales='free') + 
  theme(strip.background = element_rect(colour="black", fill="white", 
                                        size=1, linetype="solid"))

plot(p1)
ggsave(filename="~/work/ToxoplasmaGondii/ToxoR/output/AP2_trends.pdf", plot=p1,
       width = 8, height = 8, 
       units = "in", # other options are "in", "cm", "mm" 
       dpi = 300
)


##### Clustering
library(TMixClust)
getCluster <- function(D, nb_clusters = 3){
  best_clust_obj <- analyse_stability(D, nb_clusters = nb_clusters, nb_clustering_runs = 20, nb_cores = 4)
  return(best_clust_obj)
}

toxo.ap2 <- toxo.tab %>% dplyr::filter(is.TF == 'yes') %>% 
  dplyr::filter(str_detect(Product.Description, 'AP2 domain transcription factor'))
toxo.ap2$Product.Description <- gsub("AP2 domain transcription factor ", "", toxo.ap2$Product.Description)
toxo.ap2.extra.expr <- toxo.ap2 %>% 
  dplyr::select(contains('6hes_FPKM')) %>% dplyr::select(-contains("B4")) %>%
  dplyr::select(-contains("intra")) %>% dplyr::select(-contains("_sd")) %>%
  dplyr::select(-contains("freshlyse"))  %>%  dplyr::select(-contains("_quant"))

colnames(toxo.ap2.extra.expr) <- gsub('RH', 'P1000', colnames(toxo.ap2.extra.expr))
  
  
p.times <- gsub("\\.", '', gsub("B\\.", "", gsub("B2\\.", "", gsub("6hes_FPKM", "", colnames(toxo.ap2.extra.expr)))))
colnames(toxo.ap2.extra.expr) <- p.times
col.order <- sort(as.numeric(gsub("P", "", (p.times))), index.return=TRUE)$ix
toxo.ap2.extra.expr <- toxo.ap2.extra.expr[,col.order]
#rownames(toxo.ap2.expr) <- toxo.ap2$GeneName
rownames(toxo.ap2.extra.expr) <- toxo.ap2$Product.Description
colnames(toxo.ap2.extra.expr) <- gsub('P1000', 'RH', colnames(toxo.ap2.extra.expr))

toxo.ap2.extra.expr <- toxo.ap2.extra.expr[rownames(toxo.ap2.extra.expr) %in% my.AP2s$AP2, ]
ap2.clusters <- getCluster(toxo.ap2.extra.expr[,1:4], nb_clusters = 3)

toxo.ap2.extra.expr.clusters <- data.frame(AP2 = rownames(toxo.ap2.extra.expr), 
                                           toxo.ap2.extra.expr,
                                           cluster = ap2.clusters$em_cluster_assignment)


toxo.ap2.extra.expr.clusters <- toxo.ap2.extra.expr.clusters %>% 
  gather(key = Passage, value = expression, P7:RH)

toxo.ap2.extra.expr.clusters$Passage <- factor(toxo.ap2.extra.expr.clusters$Passage,
                                               levels = unique(toxo.ap2.extra.expr.clusters$Passage))

toxo.ap2.extra.expr.clusters$cluster <- factor(toxo.ap2.extra.expr.clusters$cluster,
                                               levels = unique(sort(toxo.ap2.extra.expr.clusters$cluster)))

toxo.ap2.extra.expr.clusters$labs <- ""
toxo.ap2.extra.expr.clusters$labs[toxo.ap2.extra.expr.clusters$Passage == 'P7'] <- 
  as.character(toxo.ap2.extra.expr.clusters$AP2[toxo.ap2.extra.expr.clusters$Passage == 'P7'])
toxo.ap2.extra.expr.clusters$class <- ""
toxo.ap2.extra.expr.clusters$class[toxo.ap2.extra.expr.clusters$AP2 %in% activators] <- 'Brady Activator'
toxo.ap2.extra.expr.clusters$class[toxo.ap2.extra.expr.clusters$AP2 %in% repressors] <- 'Brady Represson'
toxo.ap2.extra.expr.clusters$class[toxo.ap2.extra.expr.clusters$AP2 %in% others] <- 'Activator'

toxo.ap2.extra.expr.clusters$class <- factor(toxo.ap2.extra.expr.clusters$class,
                                               levels = unique(toxo.ap2.extra.expr.clusters$class))


p <- ggplot(data = subset(toxo.ap2.extra.expr.clusters, Passage != 'RH'), 
                          aes(x = Passage, y = expression, group = AP2, color = class)) +
  geom_line() + theme_bw() + 
  geom_text(
    size    = 2,
    mapping = aes(x = Passage, y = expression, label = labs),
    hjust   = 1.05,
    vjust   = 1.5,
    fontface = "bold"
  )+
  
  facet_wrap(.~cluster) + 
  theme(strip.background = element_rect(colour="black", fill="white", 
                                        size=1, linetype="solid"))

plot(p)  
ggsave(filename="~/work/ToxoplasmaGondii/ToxoR/output/AP2_extra_cellular_curve_clusters.pdf", plot=p,
       width = 8, height = 8, 
       units = "in", # other options are "in", "cm", "mm" 
       dpi = 300
)
