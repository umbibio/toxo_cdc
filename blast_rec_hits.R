library(dplyr)
library(orthologr)
library(seqinr)
library(readxl)
library(openxlsx)

## This requires installing Blast from NCBI

rec_sn.vs.toxo <- blast_rec(query_file   = "../Input/rec_blast/ToxoDB-57_SneuronaSN3_AnnotatedProteins.fasta",
                            subject_file = "../Input/rec_blast/ToxoDB-57_TgondiiGT1_AnnotatedProteins.fasta",
                            delete_corrupt_cds = T, seq_type = "protein",
                            format = "fasta", blast_algorithm = "blastp",
                            comp_cores = 16,
                            eval = 0.0001)

rec_sn.vs.toxo <- rec_sn.vs.toxo %>% dplyr::select(query_id, subject_id)
rec_sn.vs.toxo$query_id <- gsub('-.*', '', rec_sn.vs.toxo$query_id)
rec_sn.vs.toxo$subject_id <- gsub('-.*', '', rec_sn.vs.toxo$subject_id)
colnames(rec_sn.vs.toxo) <- c('sn_gene_id', 'Tg_gene_id')
write.xlsx(rec_sn.vs.toxo, '../Input/toxo_genomics/Orthologs/SN_vs_TG.xlsx')

rec_sn.vs.pb <- blast_rec(query_file   = "../Input/rec_blast/ToxoDB-57_SneuronaSN3_AnnotatedProteins.fasta",
                            subject_file = "../Input/rec_blast/PlasmoDB-53_PbergheiANKA_AnnotatedProteins.fasta",
                            delete_corrupt_cds = T, seq_type = "protein",
                            format = "fasta", blast_algorithm = "blastp",
                            comp_cores = 16,
                            eval = 0.0001)

rec_sn.vs.pb <- rec_sn.vs.pb %>% dplyr::select(query_id, subject_id)
rec_sn.vs.pb$query_id <- gsub('-.*', '', rec_sn.vs.pb$query_id)
rec_sn.vs.pb$subject_id <- gsub('\\..*', '', rec_sn.vs.pb$subject_id)
colnames(rec_sn.vs.pb) <- c('sn_gene_id', 'pb_gene_id')
write.xlsx(rec_sn.vs.pb, '../Input/toxo_genomics/Orthologs/SN_vs_PB.xlsx')

rec_toxo.vs.pb <- blast_rec(query_file   = "../Input/rec_blast/ToxoDB-57_TgondiiGT1_AnnotatedProteins.fasta",
                            subject_file = "../Input/rec_blast/PlasmoDB-53_PbergheiANKA_AnnotatedProteins.fasta",
                            delete_corrupt_cds = T, seq_type = "protein",
                            format = "fasta", blast_algorithm = "blastp",
                            comp_cores = 16,
                            eval = 0.0001)

rec_toxo.vs.pb <- rec_toxo.vs.pb %>% dplyr::select(query_id, subject_id)
rec_toxo.vs.pb$query_id <- gsub('-.*', '', rec_toxo.vs.pb$query_id)
rec_toxo.vs.pb$subject_id <- gsub('\\..*', '', rec_toxo.vs.pb$subject_id)
colnames(rec_toxo.vs.pb) <- c('Tg_gene_id', 'pb_gene_id')
write.xlsx(rec_toxo.vs.pb, '../Input/toxo_genomics/Orthologs/TG_vs_PB.xlsx')


### Map to toxo orthologs and transfer labels
SN.prod.desc <- read.xlsx('../Input/sarcosistis_neurona/genes/SN_ProdDescription.xlsx')
colnames(SN.prod.desc)[1] <- "sn_gene_id"
colnames(SN.prod.desc)[2] <- "ProductDescription.SN"
TG.prod.desc <- read.xlsx('../Input/toxo_genomics/genes/ProductDescription_GT1.xlsx')
colnames(TG.prod.desc)[1] <- "Tg_gene_id"
colnames(TG.prod.desc)[2] <- "ProductDescription.TG"
PB.prod.desc <- read.csv('../Input/compScBdTgPb/genes/PBer_Prod_Desc.csv')
PB.prod.desc <- PB.prod.desc %>% transmute(GeneID = Gene.ID, ProductDescription.PB = Product.Description)
colnames(PB.prod.desc)[1] <- "pb_gene_id"
colnames(PB.prod.desc)[2] <- "ProductDescription.PB"


SN.prod.desc$Tg_gene_id <- rec_sn.vs.toxo$Tg_gene_id[match(SN.prod.desc$sn_gene_id, rec_sn.vs.toxo$sn_gene_id)]
SN.prod.desc$pb_gene_id <- rec_sn.vs.pb$pb_gene_id[match(SN.prod.desc$sn_gene_id, rec_sn.vs.pb$sn_gene_id)]
SN.prod.desc$ProductDescription.TG <- TG.prod.desc$ProductDescription.TG[match(SN.prod.desc$Tg_gene_id, TG.prod.desc$Tg_gene_id)]
SN.prod.desc$ProductDescription.PB <- PB.prod.desc$ProductDescription.PB[match(SN.prod.desc$pb_gene_id, PB.prod.desc$pb_gene_id)]
SN.prod.desc <- SN.prod.desc %>% transmute(sn_gene_id = sn_gene_id, Tg_gene_id = Tg_gene_id, pb_gene_id = pb_gene_id,
                                           ProductDescription.SN = ProductDescription.SN,
                                           ProductDescription.TG = ProductDescription.TG,
                                           ProductDescription.PB = ProductDescription.PB)

TG.prod.desc$sn_gene_id <- rec_sn.vs.toxo$sn_gene_id[match(TG.prod.desc$Tg_gene_id, rec_sn.vs.toxo$Tg_gene_id)]
TG.prod.desc$pb_gene_id <- rec_toxo.vs.pb$pb_gene_id[match(TG.prod.desc$Tg_gene_id, rec_toxo.vs.pb$Tg_gene_id)]
#TG.prod.desc <- left_join(left_join(TG.prod.desc, rec_sn.vs.toxo, by = 'Tg_gene_id'), rec_toxo.vs.pb, by = 'Tg_gene_id')
TG.prod.desc$ProductDescription.SN <- SN.prod.desc$ProductDescription.SN[match(TG.prod.desc$sn_gene_id, SN.prod.desc$sn_gene_id)]
TG.prod.desc$ProductDescription.PB <- PB.prod.desc$ProductDescription.PB[match(TG.prod.desc$pb_gene_id, PB.prod.desc$pb_gene_id)]
TG.prod.desc <- TG.prod.desc %>% transmute(sn_gene_id = sn_gene_id, Tg_gene_id = Tg_gene_id, pb_gene_id = pb_gene_id,
                                           ProductDescription.SN = ProductDescription.SN,
                                           ProductDescription.TG = ProductDescription.TG,
                                           ProductDescription.PB = ProductDescription.PB)

PB.prod.desc$sn_gene_id <- rec_sn.vs.pb$sn_gene_id[match(PB.prod.desc$pb_gene_id, rec_sn.vs.pb$pb_gene_id)]
PB.prod.desc$Tg_gene_id <- rec_toxo.vs.pb$Tg_gene_id[match(PB.prod.desc$pb_gene_id, rec_toxo.vs.pb$pb_gene_id)]
#PB.prod.desc <- left_join(left_join(PB.prod.desc, rec_sn.vs.pb, by = 'pb_gene_id'), rec_toxo.vs.pb, by = 'pb_gene_id')
PB.prod.desc$ProductDescription.SN <- SN.prod.desc$ProductDescription.SN[match(PB.prod.desc$sn_gene_id, SN.prod.desc$sn_gene_id)]
PB.prod.desc$ProductDescription.TG <- TG.prod.desc$ProductDescription.TG[match(PB.prod.desc$Tg_gene_id, TG.prod.desc$Tg_gene_id)]
PB.prod.desc <- PB.prod.desc %>% transmute(sn_gene_id = sn_gene_id, Tg_gene_id = Tg_gene_id, pb_gene_id = pb_gene_id,
                                           ProductDescription.SN = ProductDescription.SN,
                                           ProductDescription.TG = ProductDescription.TG,
                                           ProductDescription.PB = ProductDescription.PB)


ALL.prod.desc <- rbind(SN.prod.desc, TG.prod.desc, PB.prod.desc) %>% distinct()
write.xlsx(ALL.prod.desc, '../Input/toxo_genomics/genes/toxo_sn_pb_prod_desc.xlsx')
saveRDS(ALL.prod.desc, '../Input/toxo_cdc/rds/toxo_sn_pb_prod_desc.rds')






