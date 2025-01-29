####################################################################
# Charles Zhang
# Gate Lab
# Northwestern University
# 03/28/2023
####################################################################
# Perform standard normalization, scaling, PCA, integration
# Expected output:
#       seurat_standard_04_01.rds
####################################################################
# See README.md for additional links

source("../00.ref/config/immune_panel_config.R")

seurat_input_dir <- glue("{out_dir_root}/03.quality_control")
input_seurat <- glue("{seurat_input_dir}/seurat_merged_02.rds")
seurat_output_dir <- glue("{out_dir_root}/04.seurat")

seurat <- readRDS(input_seurat)
dir.create(seurat_output_dir, showWarnings = FALSE, recursive = TRUE )

# smaller feature and umi pool?
seurat <- subset(seurat, subset = nFeature_RNA > 20 & nCount_RNA > 50)

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