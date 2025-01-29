####################################################################
# Charles Zhang
# Gate Lab
# Northwestern University
# 11/21/2023
####################################################################
# Perform SCTransform and map reference labels onto the merged seurat dataset
# Expected output:
#       seurat_SCT_mapped_04_02, which has the additional SCTransform assay and cell type label
####################################################################
# INSTRUCTIONS FOR USE:
# Change input_seurat to the seurat object with added clonotype and contig information
# Change output_seurat_dir to the directory to output seurat object to
####################################################################
# See README.md for additional links
source("../00.ref/config/CRISPR_clean_config.R")
library(SeuratDisk)
seurat_output_dir <- glue("{out_dir_root}/04.seurat")
input_seurat <- glue("{seurat_output_dir}/seurat_harmony_04_01.rds")
reference_pbmc <- "/projects/b1169/projects/als-project/resources/reference/pbmc_multimodal.h5seurat"
reference <- LoadH5Seurat(reference_pbmc)
seurat <- readRDS(input_seurat)

# run SCTransform on all samples together since there was no batch
# and no need for integration
seurat <- SCTransform(seurat, vst.flavor = "v2")

# rerun a denovo umap in addition to the reference umap
seurat <- RunPCA(seurat, npcs = 30, verbose = FALSE)
seurat <- RunUMAP(seurat, dims = 1:30, reduction.key = "UMAP_SCT", reduction.name = "umap_sct")


anchors <- FindTransferAnchors(
  reference = reference,
  query = seurat,
  normalization.method = "SCT",
  reference.reduction = "spca",
  dims = 1:50
)
seurat <- TransferData(
  anchorset = anchors,
  reference = reference,
  query = seurat,
  refdata = list(
    celltype.l1 = "celltype.l1",
    celltype.l2 = "celltype.l2",
    predicted_ADT = "ADT")
)

saveRDS(seurat, glue("{seurat_output_dir}/seurat_SCT_mapped_04_02.rds"))
