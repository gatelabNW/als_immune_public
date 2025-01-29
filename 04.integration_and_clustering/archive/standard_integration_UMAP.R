####################################################################
# Charles Zhang
# Gate Lab
# Northwestern University
# 03/28/2023
####################################################################
# Perform standard normalization, scaling, PCA, integration
# Expected output:
#       seurat_standard_05_01.rds
####################################################################
# See README.md for additional links

source("../00.ref/config/immune_profiling_config.R")

seurat_output_dir <- glue("{out_dir_root}/04.seurat")
input_seurat <- glue("{seurat_output_dir}/seurat_04_01.rds")

seurat <- readRDS(input_seurat)

seurat_list <- SplitObject(seurat, split.by = "orig.ident")
seurat_list <- lapply(X = seurat_list, FUN = function(x) {
  x <- NormalizeData(x, verbose = FALSE)
  FindVariableFeatures(x, verbose = FALSE)
})

features <- SelectIntegrationFeatures(object.list = seurat_list)
seurat_list <- lapply(X = seurat_list, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
})

# Samples to use as reference, better way to write it?
for (i in c(1, 6, 11, 16, 21, 26, 31, 36)) {
  print(unique(seurat_list[[i]]$orig.ident))
}

anchors <- FindIntegrationAnchors(
  object.list = seurat_list,
  reference = c(1, 6, 11, 16, 21, 26, 31, 36), reduction = "rpca",
  dims = 1:50
)
print("INFO: Anchors Found!")
# Save this object just in case the script fails afer this point.
# Script takes a long time and a lot of memory.

seurat.integrated <- IntegrateData(anchorset = anchors, dims = 1:50)
seurat.integrated <- ScaleData(seurat.integrated, verbose = FALSE)
seurat.integrated <- RunPCA(seurat.integrated, verbose = FALSE)
seurat.integrated <- RunUMAP(seurat.integrated, dims = 1:50)

saveRDS(seurat.integrated, glue("{seurat_output_dir}/seurat_standard_04_02.rds"))
