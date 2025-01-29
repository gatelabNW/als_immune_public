source("../00.ref/config/CRISPR_clean_config.R")
enrichment_output_dir <- glue("{out_dir_root}/08.enrichment")
dir.create(enrichment_output_dir, showWarnings = FALSE, recursive = TRUE )
decouple_output_dir <- glue("{out_dir_root}/08.enrichment/decoupleR")
dir.create(decouple_output_dir, showWarnings = FALSE, recursive = TRUE )

library(decoupleR)
library(dplyr)
library(pheatmap)

seurat_output_dir <- glue("{out_dir_root}/04.seurat")
input_seurat <- glue("{seurat_output_dir}/seurat_add_QC_04_04.rds")
s <- readRDS(input_seurat)

net <- get_progeny(organism = 'human', top = 500)
mat <- as.matrix(s[["RNA"]]$data)

# Run mlm
acts <- run_mlm(mat=mat, net=net, .source='source', .target='target',
                .mor='weight', minsize = 5)

write.csv(acts, file = glue("{decouple_output_dir}/mlm_progeny.csv"), row.names = F, quote = F)
