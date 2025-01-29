####################################################################
# Charles Zhang
# Gate Lab
# Northwestern University
# 11/14/2023
####################################################################

source("../00.ref/config/CRISPR_clean_config.R")
contig_clonotype_output <- glue("{out_dir_root}/04.contig_clonotype_output")
seurat_input_dir <- glue("{out_dir_root}/03.quality_control")
seurat_output_dir <- glue("{out_dir_root}/04.seurat")
seurat <- glue("{seurat_input_dir}/seurat_merged_02.rds")|>readRDS()

all_bcr_metadata <- read.csv(glue("{contig_clonotype_output}/all_bcr_meta.csv"))
all_tcr_metadata <- read.csv(glue("{contig_clonotype_output}/all_tcr_meta.csv"))

tcr_columns_to_add <- c('barcode','tcr_clonotype_id', 'tcr_frequency', 'tra_cdr3s', 'trb_cdr3s', 'inkt_evidence', 'mait_evidence')
tcr <- all_tcr_metadata|>dplyr::select(tcr_columns_to_add)
# since each cell has two entries, one for alpha chain and one for beta chain, duplicates are removed
tcr <- tcr[!duplicated(tcr),]
tcr <- as.data.frame(tcr)

# does not drop na as all cells are unique
bcr_columns_to_add <- c('barcode','bcr_clonotype_id', 'bcr_frequency', 'igl_cdr3s', 'igh_cdr3s', 'light_chain')
bcr <- all_bcr_metadata|>dplyr::select(bcr_columns_to_add)
bcr <- bcr[!duplicated(bcr),]
# some random BCR cell data have clonotypes that does not have any information in the clonotypes
# thus they are discarded
na_rows <- bcr[apply(is.na(bcr), 1, any), ]
bcr <- drop_na(bcr)
bcr <- as.data.frame(bcr)


# set rownames
rownames(bcr) <- bcr$barcode
bcr$barcode <- NULL
rownames(tcr) <- tcr$barcode
tcr$barcode <- NULL


seurat<-seurat |> AddMetaData(tcr)
seurat<-seurat |> AddMetaData(bcr)

saveRDS(seurat, glue("{seurat_output_dir}/seurat_SCT_mapped_adaptive_added_04_03.rds"))