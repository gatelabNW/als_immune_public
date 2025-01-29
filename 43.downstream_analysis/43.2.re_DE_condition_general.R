source("../00.ref/config/spatial_config.R")
options(echo = TRUE)
source(glue("{repo_root}/00.lib/plot/ggplot_formatting.R"))
source(glue("{repo_root}/00.lib/plot/volcano_plot.R"))

DEG_out_dir <- glue("{out_dir_root}/04.downstream_analysis/DE_random_effect")
input_data_dir <- glue("{out_dir_root}/03.seurat_process")
out_data_dir <- glue("{DEG_out_dir}/condition_general/data")
out_plot_dir <- glue("{DEG_out_dir}/condition_general/plot")

dir.create(DEG_out_dir, recursive = T, showWarnings = F)
dir.create(out_data_dir, recursive = T, showWarnings = F)
dir.create(out_plot_dir, recursive = T, showWarnings = F)

args <- commandArgs(trailingOnly = TRUE)
# cluster <- 0
cluster <- as.numeric(args[1])

comparisons <- c("ALS_vs_Control")
cluster_col <- "integrated_snn_res.0.25"
p_thresh <- 0.05
fc_thresh <- log2(1.5)
condition_col <- "condition_general"

s <- readRDS(glue("{input_data_dir}/data/all_samples_seurat_manual_02.rds"))
DefaultAssay(s) <- "Spatial"
s <- JoinLayers(s)
cdr <- colMeans(GetAssayData(s, assay = "Spatial", layer = "counts") > 0)
s@meta.data$cdr <- cdr

# Make sure default assay is SCT
DefaultAssay(object = s) <- "SCT"

# generate condition_general column
clusters <- s@meta.data[[cluster_col]] |> unique() |> sort()
s[[condition_col]] <- sapply(s$condition, function(x){
  if(x == "Control"){
    out <- x
  }else if(x == "sALS"){
    out <- "ALS"
  }else if(x == "C9orf72"){
    out <- "ALS"
  }else{
    out <- "NA"
  }
  out
})

# MAST DE
print(cluster)

# Subset object to current cluster and recorrect SCT data - keep all NNC spots like we did for DESeq2 (to keep ROI the same)
clust_s <- subset(s, (integrated_snn_res.0.25 == cluster))

# TODO: figure out why this is necessary?
for (name in names(clust_s@assays$SCT@SCTModel.list)) {
  clust_s@assays$SCT@SCTModel.list[[name]]@umi.assay <- 'Spatial'
}
clust_s <- PrepSCTFindMarkers(object = clust_s)

clust_s@meta.data$cdr_centered <- as.numeric((clust_s@meta.data$cdr - mean(clust_s@meta.data$cdr)) / sd(clust_s@meta.data$cdr))
clust_s@meta.data$sample_id <- factor(clust_s@meta.data$sample)

# Get all expression data
expressionmat_full <- GetAssayData(clust_s, assay = "SCT", layer = 'data')
expressionmat_full <- as.matrix(expressionmat_full)

# set idents
ident.1 <- "ALS"
ident.2 <- "Control"

# Get avg_log2FC and percent expression for all genes
Idents(clust_s) <- condition_col
LFC <- FoldChange(object = clust_s, ident.1 = ident.1, ident.2 = ident.2, fc.name = "avg_log2FC", base = 2) %>% data.frame()

# Manual min.pct filtering
genes_keep_10pct <- row.names(LFC)[LFC$pct.1 >= 0.1 | LFC$pct.2 >= 0.1]
genes_keep_1pct <- row.names(LFC)[LFC$pct.1 >= 0.01 & LFC$pct.2 >= 0.01]
genes_keep <- intersect(genes_keep_10pct, genes_keep_1pct)

# Remove contamination
if (length(grep(pattern = "^RPS|^RPL|^MT-|^HB", x = genes_keep)) != 0) {
  genes_keep <- genes_keep[-grep(pattern = "^RPS|^RPL|^MT-|^HB", x = genes_keep)]
}

# Filter expression matrix
expressionmat <- expressionmat_full[row.names(expressionmat_full) %in% genes_keep,]

# Cell-level and feature-level meta data for MAST
cdat <- clust_s@meta.data
cdat$wellKey <- row.names(cdat)
fdat <- data.frame(primerid = row.names(expressionmat), row.names = row.names(expressionmat))

# Subset to conditions for comparison
cdat <- cdat[cdat[[condition_col]] %in% c(ident.1, ident.2),]
expressionmat <- expressionmat[,colnames(expressionmat) %in% row.names(cdat)]

# Create MAST object
sca <- FromMatrix(expressionmat, cdat, fdat, check_sanity = FALSE)

# Set reference level for condition
cond <- factor(colData(sca)[[condition_col]])
# TODO: figure out why setting reference as group1, which is not control, in the original script
cond <- relevel(cond, ident.2) # FindMarkers: reference level is "Group1"
colData(sca)[[condition_col]] <- cond

tryCatch(
{
  # MAST does not allow ebayes (empirical Bayes estimation to regularize variance) for glmer
  # zlm_condition <- zlm(~ condition + sex + age_centered + cdr_centered + gDNA_centered + (1 | sample_id), sca,
  #                      ebayes = FALSE, method = "glmer")
  print("INFO: fitting model!")
  zlm_condition <- zlm(~ condition_general + cdr_centered + (1 | sample_id), sca,
                       ebayes = FALSE, method = "glmer")
  lrt_name <- paste0("condition_general", ident.1)

  print("INFO: performing likelihood ratio test!")
  summary_condition <- summary(zlm_condition, doLRT = lrt_name) # FindMarkers: doLRT = "conditionGroup2"

  # From FindMarkers code
  summary_data <- summary_condition$datatable %>% data.frame()
  p_val <- summary_data[summary_data[, "component"] == "H", 4]
  genes.return <- summary_data[summary_data[, "component"] == "H", 1]

  # Compile results
  results <- data.frame(p_val = p_val, gene = genes.return, row.names = genes.return)
  results$BH <- p.adjust(results$p_val, method = "BH")

  LFC_use <- LFC[match(results$gene, row.names(LFC)),]
  results$avg_log2FC <- LFC_use$avg_log2FC

  results$DE <- "Not DE"
  results$DE[results$avg_log2FC > fc_thresh & results$BH < p_thresh] <- "Upregulated"
  results$DE[results$avg_log2FC < -fc_thresh & results$BH < p_thresh] <- "Downregulated"

  # For volcano plot labels
  results$DE_gene <- NA
  results$DE_gene[results$DE != "Not DE"] <- results$gene[results$DE != "Not DE"]

  out_file <- glue("{out_data_dir}/{cluster}__{ident.1}_vs._{ident.2}.csv")
  write.csv(results, out_file, quote = F)
  print(glue("INFO: saved output to {out_file}!"))
},
  error = function(e) {
    print(paste0("Error in ", cluster, ": ", comparison))
  }
)



