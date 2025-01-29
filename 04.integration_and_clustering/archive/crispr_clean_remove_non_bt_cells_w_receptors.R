####################################################################
# Charles Zhang
# Gate Lab
# Northwestern University
# 02/15/2023

# See README.md for additional links
source("../00.ref/config/CRISPR_clean_config.R")
seurat_output_dir <- glue("{out_dir_root}/04.seurat")
input_seurat <- glue("{seurat_output_dir}/seurat_SCT_mapped_04_02.rds")
seurat <- readRDS(input_seurat)

bcr_table <- dplyr::select(seurat@meta.data, predicted.celltype.l1, bcr_frequency)
non_b_cells_to_discard <- bcr_table[(bcr_table$predicted.celltype.l1!="B") & (!is.na(bcr_table$bcr_frequency)),]
non_b_cells_to_discard <- non_b_cells_to_discard|>rownames()

tcr_table <- dplyr::select(seurat@meta.data, predicted.celltype.l1, tcr_frequency)
non_t_cells_to_discard <- tcr_table[((tcr_table$predicted.celltype.l1!="CD4 T")&(tcr_table$predicted.celltype.l1!="CD8 T")&(tcr_table$predicted.celltype.l1!="other T")) & (!is.na(tcr_table$tcr_frequency)),]
non_t_cells_to_discard <- non_t_cells_to_discard|>rownames()
print(glue("INFO: Non B cells with BCR removed {length(non_b_cells_to_discard)}!"))
print(glue("INFO: Non T cells with TCR removed {length(non_t_cells_to_discard)}!"))

all_cells_to_discard <- c(non_b_cells_to_discard, non_t_cells_to_discard)

seurat_filtered <- seurat[,!colnames(seurat) %in% all_cells_to_discard]
saveRDS(seurat_filtered, glue("{seurat_output_dir}/seurat_add_QC_04_03.rds"))
