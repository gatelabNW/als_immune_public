####################################################################
# Charles Zhang
# Gate Lab
# Northwestern University
# 09/25/2023
####################################################################
# Run general QC on soupX corrected count matrix for ALS project
# Since immune panel has significantly less genes, the min feature threshold is set to 20 genes
# Expected output:
#     seurat_merged_01.rds: a seurat object before any QC is perfomed on
#     seurat_merged_02.rds: a seurat object with only singlets, MT<10
####################################################################
# INSTRUCTIONS FOR USE:
# Change soupX_count_mtx_root to the root directory holding all samples' soupX output
# Change out_root_dir to a directory for the new seurat objects
# Change meta to the path to the corrected sample meta data
# Change qc_plot_dir to output directory for basic qc plots
####################################################################
source("../00.ref/config/immune_panel_config.R")


# soupX corrected count matrices root directory
soupX_count_mtx_root <-  glue("{out_dir_root}/02.soupX/count")

# ouput root directory
out_root_dir <- glue("{out_dir_root}/03.quality_control")
dir.create(out_root_dir, showWarnings = FALSE, recursive = TRUE )

# plot directory for QC plots
qc_plot_dir <- glue("{out_dir_root}/03.quality_control/qc_plots/")
dir.create(qc_plot_dir, showWarnings = FALSE, recursive = TRUE )


run_doubletfinder <- function(s) {

  # TODO: Find out why this value was chosen
  # IIRC, those docs said to use 0.044
  # TODO: Link to this doc in a comment
  # TODO: Also do this for all other params in this file
  doublet_formation_rate <- 0.044
  print(paste0("Using doublet formation rate of ", doublet_formation_rate))

  s <- NormalizeData(s)
  s <- ScaleData(s)
  s <- FindVariableFeatures(s, selection.method = "vst", nfeatures = 2000)
  s <- RunPCA(s)
  s <- FindNeighbors(s, dims = 1:12)
  s <- FindClusters(s, resolution = 0.3)
  s <- RunTSNE(s, dims = 1:12)

  ## pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
  sweep.res.list <- paramSweep(s, PCs = 1:12, sct = FALSE)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn <- find.pK(sweep.stats) # mean-variance-normalized bimodality coefficient, find maxBCmvn per pK
  max_index <- which.max(bcmvn$BCmetric)
  optimal_pK <- as.numeric(as.character(bcmvn[max_index, "pK"]))

  ## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
  annotations <- s@meta.data$seurat_clusters
  homotypic.prop <- modelHomotypic(annotations)
  nExp_poi <- round(doublet_formation_rate*nrow(s@meta.data))


  ## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
  s <- doubletFinder(s, PCs = 1:12, pN = 0.25, pK = optimal_pK, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

  ## Rename column name for consistency
  colnames(s@meta.data)[ grep("DF.classifications*", colnames(s@meta.data)) ] <- "DF.classifications"
  return(s)
}


samples <- list.files(soupX_count_mtx_root)

seurat_object_list <- samples |> lapply(\(sample) {
  CreateSeuratObject(
    counts = "{soupX_count_mtx_root}/{sample}" |> glue() |> Read10X(),
    project =  sample,
    min.cells = 3,
    min.features = 20
  )
})
seurat_object_list <- seurat_object_list |> lapply(run_doubletfinder)

seurat <- merge(seurat_object_list[[1]], seurat_object_list[2:length(seurat_object_list)],
                add.cell.ids = samples,
                project = "ALS")

# add sample meta data
samples_metadata <- read.csv(meta)
rownames(samples_metadata) <- samples_metadata$gex_index
samples_metadata[["orig.ident"]] <- samples_metadata$gex_index
x_to_add <- left_join(x=seurat[["orig.ident"]], y=samples_metadata)
row.names(x_to_add) <- row.names(seurat[[]])
seurat <- AddMetaData(seurat, metadata = x_to_add)
print(table(seurat@meta.data[,c("sample_id","diagnosis")]))

# Save seurat before removing doublets
saveRDS(seurat, file = glue('{out_root_dir}/seurat_merged_01.rds'))

# Subset for singlets
seurat <- seurat |> subset(subset = DF.classifications == "Singlet")
print("After removing doublets...")
seurat

print("Before removing high mitochondrial content cells...")
seurat
seurat[["percent.mt"]] <- PercentageFeatureSet(seurat, pattern = "^MT-")
seurat <- subset(seurat, subset = percent.mt < 10)
print("After removing high mitochondrial content cells...")
seurat

saveRDS(seurat, file = glue('{out_root_dir}/seurat_merged_02.rds'))

# QC plots
# nFeatures per ID
Idents(seurat)<-"orig.ident"
myplot <- VlnPlot(seurat, features = c("nFeature_RNA"), pt.size = 0) +
  theme_Publication_blank() +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 8)) +
  labs(x = "Sample ID",
       y = "Number of Features",
       title = "Number of Features per Sample ID")
set_panel_size(myplot, file=paste0(qc_plot_dir, "nFeature_ID.pdf"),
               width=unit(8, "in"), height=unit(4, "in"))

# nCounts per ID
myplot <- VlnPlot(seurat, features = c("nCount_RNA"), pt.size = 0) +
  theme_Publication_blank() +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 8)) +
  labs(x = "Sample ID",
       y = "Number of Counts",
       title = "Number of Counts per Sample ID")
set_panel_size(myplot, file=paste0(qc_plot_dir, "nCount_ID.pdf"),
               width=unit(8, "in"), height=unit(4, "in"))

# percent.mt per ID
myplot <- VlnPlot(seurat, features = c("percent.mt"), pt.size = 0) +
  theme_Publication_blank() +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 8)) +
  labs(x = "Sample ID",
       y = "Percent Mitochondrial",
       title = "Percent Mitochondrial per Sample ID")
set_panel_size(myplot, file=paste0(qc_plot_dir, "percent.mt_ID.pdf"),
               width=unit(8, "in"), height=unit(4, "in"))

# nFeatures per Diagnosis
Idents(seurat) <- 'diagnosis_general'
myplot <- VlnPlot(seurat, features = c("nFeature_RNA"), pt.size = 0) +
  theme_Publication_blank() +
  theme(legend.position = "none") +
  labs(x = "Diagnosis",
       y = "Number of Features",
       title = "Number of Features per Diganosis")
set_panel_size(myplot, file=paste0(qc_plot_dir, "nFeature_Diagnosis_General.pdf"),
               width=unit(3, "in"), height=unit(4, "in"))

# nFeatures per sub diagnosis
Idents(seurat) <- 'diagnosis'
myplot <- VlnPlot(seurat, features = c("nFeature_RNA"), pt.size = 0) +
  theme_Publication_blank() +
  theme(legend.position = "none") +
  labs(x = "Diagnosis",
       y = "Number of Features",
       title = "Number of Features per Diganosis")
set_panel_size(myplot, file=paste0(qc_plot_dir, "nFeature_Diagnosis.pdf"),
               width=unit(3, "in"), height=unit(4, "in"))

Idents(seurat) <- 'diagnosis_general'
myplot <- VlnPlot(seurat, features = c("nCount_RNA"), pt.size = 0) +
  theme_Publication_blank() +
  theme(legend.position = "none") +
  labs(x = "Diagnosis",
       y = "Number of Counts",
       title = "Number of Counts per Diganosis")
set_panel_size(myplot, file=paste0(qc_plot_dir, "nCount_Diagnosis_General.pdf"),
               width=unit(3, "in"), height=unit(4, "in"))

# nFeatures per sub diagnosis
Idents(seurat) <- 'diagnosis'
myplot <- VlnPlot(seurat, features = c("nCount_RNA"), pt.size = 0) +
  theme_Publication_blank() +
  theme(legend.position = "none") +
  labs(x = "Diagnosis",
       y = "Number of Counts",
       title = "Number of Counts per Diganosis")
set_panel_size(myplot, file=paste0(qc_plot_dir, "nCount_Diagnosis.pdf"),
               width=unit(3, "in"), height=unit(4, "in"))

Idents(seurat) <- 'diagnosis_general'
myplot <- VlnPlot(seurat, features = c("percent.mt"), pt.size = 0) +
  theme_Publication_blank() +
  theme(legend.position = "none") +
  labs(x = "Diagnosis",
       y = "Percent Mitochondrial",
       title = "Percent Mitochondrial per Diganosis")
set_panel_size(myplot, file=paste0(qc_plot_dir, "percent.mt_Diagnosis_General.pdf"),
               width=unit(3, "in"), height=unit(4, "in"))

Idents(seurat) <- 'diagnosis'
myplot <- VlnPlot(seurat, features = c("percent.mt"), pt.size = 0) +
  theme_Publication_blank() +
  theme(legend.position = "none") +
  labs(x = "Diagnosis",
       y = "Percent Mitochondrial",
       title = "Percent Mitochondrial per Diganosis")
set_panel_size(myplot, file=paste0(qc_plot_dir, "percent.mt_Diagnosis.pdf"),
               width=unit(3, "in"), height=unit(4, "in"))