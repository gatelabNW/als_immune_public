# Charles Zhang
# Gate Lab
# Northwestern University
# 11/28/2023
####################################################################
# Convert Seurat to AnnData, use R/4.1.1
####################################################################
dyn.load("/software/hdf5/1.8.19-serial/lib/libhdf5_hl.so.10")

library(Seurat)
library(SeuratDisk)
library(glue)
target_s <- "/projects/b1042/Gate_Lab/projects/als-project/crispr_clean_final/04.seurat/seurat_add_QC_04_04.rds"
out_root <- "/projects/b1042/Gate_Lab/projects/als-project/crispr_clean_final/04.anndata"
out_prefix <- "seurat_add_QC_04_04"
s <- readRDS(target_s)
DefaultAssay(s) <- "RNA"
SaveH5Seurat(s, filename = glue("{out_root}/{out_prefix}.h5Seurat"))
Convert(glue("{out_root}/{out_prefix}.h5Seurat"), dest = "h5ad")
