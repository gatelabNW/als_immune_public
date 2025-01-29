####################################################################
# Charles Zhang
# Gate Lab
# Northwestern University
# 03/28/2023
####################################################################
# Perform standard normalization, scaling, PCA, integration
# Expected output:
#       seurat_MNN_04_02.rds
####################################################################
# See README.md for additional links
# reference: https://satijalab.org/seurat/articles/seurat5_integration

source("../00.ref/config/CRISPR_clean_config.R")
library(SeuratWrappers)
library(harmony)
seurat_output_dir <- glue("{out_dir_root}/04.seurat")
seurat_input_dir <- glue("{out_dir_root}/03.quality_control")
input_seurat <- glue("{seurat_input_dir}/seurat_merged_02.rds")
seurat <- readRDS(input_seurat)
dir.create(seurat_output_dir, showWarnings = FALSE, recursive = TRUE )


seurat <- subset(seurat, subset = nFeature_RNA > 250 & nCount_RNA > 500)

seurat <- JoinLayers(seurat)
seurat[["RNA"]]$data <- NULL
seurat[["RNA"]]$scale.data <- NULL
seurat[["RNA"]] <- split(seurat[["RNA"]], f = seurat$orig.ident)
seurat <- NormalizeData(seurat)
seurat <- FindVariableFeatures(seurat)
seurat <- ScaleData(seurat)
seurat <- RunPCA(seurat)

seurat <- IntegrateLayers(
  object = seurat, method = HarmonyIntegration,
  orig.reduction = "pca", new.reduction = "harmony",
  verbose = FALSE
)

seurat <- FindNeighbors(seurat, reduction = "harmony", dims = 1:30)
seurat <- FindClusters(seurat, resolution = 0.5, cluster.name = "harmony_clusters")
seurat <- RunUMAP(seurat, reduction = "harmony", dims = 1:30, reduction.name = "harmony.umap")

saveRDS(seurat, glue("{seurat_output_dir}/seurat_harmony_04_01.rds"))
