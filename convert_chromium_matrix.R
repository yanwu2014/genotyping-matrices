### R script to convert mtx matrix to tab separated matrix ###

library(Seurat)
library(methods)

# Human genome
sample_path <- commandArgs(trailingOnly = T)[[1]]

matrix_path <- paste(sample_path,'/outs/filtered_gene_bc_matrices_mex/GRCh38/', sep = '')
out_file_name <- paste(sample_path, '.counts.tsv', sep = '')

# Load matrix
out_mat <- as.matrix(Read10X(matrix_path))

# Remove cells with duplicate barcodes
barcodes <- sapply(colnames(out_mat), function(x) strsplit(x, split = "-")[[1]][[1]])
dup <- duplicated(barcodes) | duplicated(barcodes, fromLast = T)
out_mat <- out_mat[, which(!dup)]

print(dim(out_mat))
write.table(out_mat, file = out_file_name, sep = '\t')
