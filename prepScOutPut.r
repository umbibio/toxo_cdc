# this scripts reads the filtered_feature_bc_matrix generated by Cellranger Count
# and generates the feature_count matrix and save it in RData folder.

library(Matrix)
                                                                                                
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
