# Charles Zhang
# Gate Lab
# Northwestern University
# 11/21/2023
####################################################################
####################################################################
dyn.load("/software/hdf5/1.8.19-serial/lib/libhdf5_hl.so.10")

library(Seurat)
library(SCP)
require(stats)
source("../00.ref/config/spatial_config.R")



# convert reference 
s <- readRDS("/projects/b1042/Gate_Lab/projects/als-project/spatial/00.ref/mg_ref.rds")
DefaultAssay(s) <- "RNA"
adata <- SCP::srt_to_adata(s, assay_X = "RNA", slot_X = "counts")
adata$write_h5ad(glue("/projects/b1042/Gate_Lab/projects/als-project/spatial/00.ref/mg_ref.h5ad"))









# convert data

# out_dir <- glue("{out_dir_root}/03.seurat_process_final")
# out_data_dir <- glue("{out_dir}/data")
# s <- readRDS(glue("{out_data_dir}/seurat_42.4__04.rds"))
# DefaultAssay(s) <- "Spatial"
# 
# s <- JoinLayers(s)
# 
# 
# cdr <- colMeans(GetAssayData(s, assay = "Spatial", layer = "counts") > 0)
# s@meta.data$cdr <- cdr
# s@meta.data$cdr_centered <- scale(s@meta.data$cdr)|>as.vector()
# 
# s <- JoinLayers(s)
# 
# assay_v3 <- CreateAssayObject(
#   counts = s[["Spatial"]]$counts
# )
# s@assays$Spatial <- assay_v3
# adata <- SCP::srt_to_adata(s, assay_X = "Spatial", slot_X = "counts")
# adata$write_h5ad(glue("{out_data_dir}/seurat_42.4__04.h5ad"))






#save a gray only version
# s_gray <- subset(s, subset = manual == "Grey_matter")
# s_gray <- JoinLayers(s_gray)
# s_gray@assays$Spatial <- as(s_gray@assays$Spatial, Class = "Assay")
# adata <- SCP::srt_to_adata(s_gray, assay_X = "Spatial", slot_X = "counts")
# adata$write_h5ad(glue("{out_data_dir}/all_samples_seurat_gray_only_0919_all_slides.h5ad"))

# sceasy::convertFormat(combined.sct.s, from="seurat", to="anndata",
#                       outFile="{out_data_dir}/all_samples_seurat_01.h5ad")
# combined.sct.s[["RNA"]] <- as(object = combined.sct.s[["Spatial"]], Class = "Assay")
# SaveH5Seurat(combined.sct.s, filename = glue("{out_data_dir}/all_samples_seurat_01.h5Seurat"))
# Convert(glue("{out_data_dir}/all_samples_seurat_01.h5Seurat"), dest = "h5ad", overwrite = TRUE)
