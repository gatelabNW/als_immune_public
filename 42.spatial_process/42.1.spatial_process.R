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
print(sessionInfo())
source("../00.ref/config/spatial_config.R")
out_dir <- glue("{out_dir_root}/03.seurat_process_final")
out_data_dir <- glue("{out_dir}/data")
out_plot_dir <- glue("{out_dir}/plot")
dir.create(out_data_dir, recursive = T, showWarnings = F)
# dir.create(out_plot_dir, recursive = T, showWarnings = F)
options(future.globals.maxSize = 8000 * 1024^2)

spaceranger_output_dir <- glue("{out_dir_root}/02.spaceranger_output_original_final")
`%notin%` <- Negate(`%in%`)


# load the meta data QC file, we are performing per sample QC
spatial_qc_meta_df <- read.csv(spatial_qc_meta)


# make a conversion to convert donor_ID to diagnosis
sample2group <- spatial_qc_meta_df$group_ID
names(sample2group) <- spatial_qc_meta_df$donor_ID

# keep only the samples that are pass qc
samples_to_use <- spatial_qc_meta_df[["spaceranger_id"]]
samples_to_use <- samples_to_use[samples_to_use != "discarded"]
all_samples <- glue("{spaceranger_output_dir}/{samples_to_use}")



# loop over each of the sample
seurat_objects <- list()
for(cur_sample in all_samples){

  # perform necessary processing
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

  # sample barcode
  cur_s@meta.data[["sample_barcode"]] <- rownames(cur_s@meta.data)

  # add condition
  cur_s_sample_short <- gsub("\\.", "-", strsplit(cur_slice, "___")[[1]][1])
  cur_group <- sample2group[[cur_s_sample_short]]
  cur_s[["condition"]] <- cur_group

  # lognormalize count data
  cur_s <- NormalizeData(cur_s, assay = "Spatial")

  # perform per sample QC
  num_spots_before_per_sample_QC <- dim(cur_s)[2]
  cur_slice_converted <- gsub("\\.", "-", cur_slice)
  if(!cur_slice_converted %in% spatial_qc_meta_df$spaceranger_id){
    print(glue("Skipped {cur_slice_converted} due to low quality!"))
    next
  }

  cur_qc_row <- dplyr::filter(spatial_qc_meta_df, spaceranger_id == cur_slice_converted)
  cur_max_umi <- cur_qc_row$max_umi
  cur_min_umi <- cur_qc_row$min_umi
  cur_max_feature <- cur_qc_row$max_feature
  cur_min_feature <- cur_qc_row$min_feature
  print(glue("{cur_slice}: {cur_max_umi}  {cur_min_umi}  {cur_max_feature}  {cur_min_feature} QC metric!"))


  if(cur_max_umi == 0){
    print(glue("Skipped {cur_slice_converted} due to low quality (0)!"))
    next
  }

  # filter for spots to be discarded
  spots_to_discard_df <-cur_s@meta.data[which((cur_s$nCount_Spatial<cur_min_umi) |
                                                (cur_s$nCount_Spatial>cur_max_umi) |
                                                (cur_s$nFeature_Spatial<cur_min_feature) |
                                                (cur_s$nFeature_Spatial>cur_max_feature)),]

  print(dim(spots_to_discard_df))

  # get the barcodes and remove them from the object
  spots_to_discard <- spots_to_discard_df|>rownames()
  cur_s <- subset(cur_s, subset = sample_barcode %notin% spots_to_discard)

  # remove MT
  cur_s[["percent.mt"]] <- PercentageFeatureSet(cur_s, pattern = "^MT-")
  cur_s <- subset(cur_s, subset = percent.mt < 20)

  # remove empty spot, if any
  keep.spots <- colnames(cur_s[, colSums(cur_s)!=0])
  cur_s <-subset(cur_s, cells = keep.spots)
  cur_s <- JoinLayers(cur_s)


  # SCTransform
  cur_s <- SCTransform(cur_s, assay = "Spatial",vars.to.regress = c("percent.mt"), variable.features.n = 8000)
  seurat_objects[[cur_slice]] <- cur_s
}

# seurat_objects <- readRDS(glue("{out_data_dir}/all_samples_seurat_list.rds"))
saveRDS(seurat_objects, glue("{out_data_dir}/all_samples_seurat_list.rds"))
s <- merge(x = seurat_objects[[1]],
           y = seurat_objects[2:length(seurat_objects)])
saveRDS(s, glue("{out_data_dir}/seurat_42.1__01.rds"))












# features <- SelectIntegrationFeatures(object.list = seurat_objects, nfeatures = 1000)
# seurat_objects <- PrepSCTIntegration(object.list = seurat_objects, anchor.features = features)
# anchors <- FindIntegrationAnchors(object.list = seurat_objects, normalization.method = "SCT",
#                                   anchor.features = features)
#
# combined.sct.s <- IntegrateData(anchorset = anchors, normalization.method = "SCT")
# combined.sct.s <- RunPCA(combined.sct.s, verbose = FALSE)
# combined.sct.s <- RunUMAP(combined.sct.s, reduction = "pca", dims = 1:30, verbose = FALSE)
# combined.sct.s <- FindNeighbors(combined.sct.s, reduction = "pca", dims = 1:30)
# combined.sct.s <- FindClusters(combined.sct.s, resolution = seq(0.1, 0.5, 0.05))
#
# saveRDS(combined.sct.s, glue("{out_data_dir}/seurat_42.1__01.rds"))