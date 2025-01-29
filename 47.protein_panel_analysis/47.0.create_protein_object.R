# ------------------------------------------------------------------------------
# -----                                                                    -----
# -----                           Gate Lab                                 -----
# -----                     Northwestern University                        -----
# -----                                                                    -----
# ------------------------------------------------------------------------------
#
# Date: 1-17-2024
# Written by: Charles Zhang
# Summary: Basic QC and SCTransform V2 integration
#
#-------------------------------------------------------------------------------

source("../00.ref/config/spatial_config.R")
out_dir <- glue("{out_dir_root}/03.seurat_process")
out_data_dir <- glue("{out_dir}/data")
out_plot_dir <- glue("{out_dir}/plot")
dir.create(out_data_dir, recursive = T, showWarnings = F)
# dir.create(out_plot_dir, recursive = T, showWarnings = F)

spaceranger_output_dir <- glue("{out_dir_root}/02.spaceranger_output_v4_0404")
`%notin%` <- Negate(`%in%`)

all_samples <- list.dirs(spaceranger_output_dir, recursive = F, full.names = T)

spatial_sample_df <- read.csv(spatial_sample_meta)|>
  dplyr::select(Donor.ID, Group.ID)
sample2group <- spatial_sample_df$Group.ID
names(sample2group) <- spatial_sample_df$Donor.ID

seurat_objects <- list()
for(cur_sample in all_samples){
  cur_slice <- basename(cur_sample)
  print(cur_slice)
  cur_dir <- glue("{cur_sample}/{cur_slice}/outs")
  # replace dash with underscore
  cur_slice <- str_replace_all(cur_slice, "-", ".")
  cur_s <- Seurat::Load10X_Spatial(cur_dir, slice = cur_slice)
  cur_s[["barcode"]] <- cur_s@meta.data|>rownames()
  cur_s@images[[cur_slice]]@coordinates[["tissue"]] <- as.integer(cur_s@images[[cur_slice]]@coordinates[["tissue"]])
  cur_s@images[[cur_slice]]@coordinates[["row"]] <- as.integer(cur_s@images[[cur_slice]]@coordinates[["row"]])
  cur_s@images[[cur_slice]]@coordinates[["col"]] <- as.integer(cur_s@images[[cur_slice]]@coordinates[["col"]])
  cur_s@images[[cur_slice]]@coordinates[["imagerow"]] <- as.integer(cur_s@images[[cur_slice]]@coordinates[["imagerow"]])
  cur_s@images[[cur_slice]]@coordinates[["imagecol"]] <- as.integer(cur_s@images[[cur_slice]]@coordinates[["imagecol"]])

  # add sample id
  cur_s[["sample"]] <- cur_slice

  # add sample id as prefix to cell id
  cur_s <- RenameCells(object = cur_s, add.cell.id = glue("{cur_slice}_"))

  # add condition
  cur_s_sample_short <- gsub("\\.", "-", strsplit(cur_slice, "___")[[1]][1])
  print(cur_s_sample_short)
  cur_group <- sample2group[[cur_s_sample_short]]
  cur_s[["sample_short"]] <- cur_s_sample_short
  cur_s[["condition"]] <- cur_group
  seurat_objects[[cur_slice]] <- cur_s
}

combined.sct.s <- merge(seurat_objects[[1]], seurat_objects[2:length(seurat_objects)])
saveRDS(combined.sct.s, glue("{out_data_dir}/all_samples_seurat_protein_01.rds"))