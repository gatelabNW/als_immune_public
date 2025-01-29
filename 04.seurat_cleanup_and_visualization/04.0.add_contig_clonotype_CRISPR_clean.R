####################################################################
# Charles Zhang
# Gate Lab
# Northwestern University
# 11/14/2023
####################################################################

# TODO: add after mapping and correct clonality information

source("../00.ref/config/CRISPR_clean_config.R")
library(dplyr)
library(tidyr)
contig_clonotype_output <- glue("{out_dir_root}/04.contig_clonotype_output")
seurat_output_dir <- glue("{out_dir_root}/04.seurat")
seurat <- glue("{seurat_output_dir}/seurat_SCT_mapped_04_02.rds")|>readRDS()
all_bcr_metadata <- read.csv(glue("{contig_clonotype_output}/all_bcr_meta.csv"))
all_tcr_metadata <- read.csv(glue("{contig_clonotype_output}/all_tcr_meta.csv"))

# Take out tcr data to be added to the single cell object
tcr_columns_to_add <- c('barcode','tcr_clonotype_id', 'tcr_frequency', 'tra_cdr3s', 'trb_cdr3s', 'inkt_evidence',
                        'mait_evidence', 'v_gene', 'd_gene', 'j_gene', 'c_gene','length')
tcr <- all_tcr_metadata|>dplyr::select(all_of(tcr_columns_to_add))
colnames(tcr) <- c('barcode','tcr_clonotype_id', 'tcr_frequency', 'tra_cdr3s', 'trb_cdr3s', 'inkt_evidence',
                   'mait_evidence', 'tcr_v_gene', 'tcr_d_gene', 'tcr_j_gene', 'tcr_c_gene','tcr_length')
# split out alpha chains
alpha_tcr <- tcr |>
  dplyr::filter(str_detect(tcr_v_gene, "TRA"))
colnames(alpha_tcr) <- c('barcode','tcr_clonotype_id', 'tcr_frequency', 'tra_cdr3s', 'trb_cdr3s', 'inkt_evidence',
                         'mait_evidence', 'tcr_v_gene_alpha', 'tcr_d_gene_alpha', 'tcr_j_gene_alpha', 'tcr_c_gene_alpha','tcr_length_alpha')
# split out beta chains
beta_tcr <- tcr |>
  dplyr::filter(str_detect(tcr_v_gene, "TRB"))|>
  dplyr::select('barcode','tcr_v_gene', 'tcr_d_gene', 'tcr_j_gene', 'tcr_c_gene','tcr_length')
colnames(beta_tcr) <- c('barcode','tcr_v_gene_beta', 'tcr_d_gene_beta', 'tcr_j_gene_beta', 'tcr_c_gene_beta','tcr_length_beta')
# merge
tcr_wide <- merge(alpha_tcr, beta_tcr)|>
  arrange(barcode)


# DO the same thing for BCR
bcr_columns_to_add <- c('barcode','bcr_clonotype_id', 'bcr_frequency', 'igl_cdr3s', 'igh_cdr3s', 'light_chain',
                        'v_gene', 'd_gene', 'j_gene', 'c_gene','length')
bcr <- all_bcr_metadata|>dplyr::select(all_of(bcr_columns_to_add))
colnames(bcr) <- c('barcode','bcr_clonotype_id', 'bcr_frequency', 'igl_cdr3s', 'igh_cdr3s', 'light_chain',
                   'bcr_v_gene', 'bcr_d_gene', 'bcr_j_gene', 'bcr_c_gene','bcr_length')
heavy_bcr <- bcr |>
  dplyr::filter(str_detect(bcr_v_gene, "IGH"))
colnames(heavy_bcr) <- c('barcode','bcr_clonotype_id', 'bcr_frequency', 'igl_cdr3s', 'igh_cdr3s', 'light_chain',
                         'bcr_v_gene_heavy', 'bcr_d_gene_heavy', 'bcr_j_gene_heavy', 'bcr_c_gene_heavy','bcr_length_heavy')
light_bcr <- bcr |>
  dplyr::filter(str_detect(bcr_v_gene, "IGL|IGK"))|>
  dplyr::select('barcode','bcr_v_gene', 'bcr_d_gene', 'bcr_j_gene', 'bcr_c_gene','bcr_length')
colnames(light_bcr) <- c('barcode','bcr_v_gene_light', 'bcr_d_gene_light', 'bcr_j_gene_light', 'bcr_c_gene_light','bcr_length_light')
bcr_wide <- merge(heavy_bcr, light_bcr)|>
  arrange(barcode)

# format to be compatible for adding metadata
mt_base <- seurat@meta.data|>dplyr::select(orig.ident, sample_id)
mt_base[["barcode"]] <- rownames(mt_base)
bcr_wide_mt <- left_join(mt_base, bcr_wide)
bcr_wide_mt$orig.ident <- NULL
bcr_wide_mt$sample_id <- NULL
rownames(bcr_wide_mt) <- bcr_wide_mt$barcode
bcr_wide_mt$barcode <- NULL

tcr_wide_mt <- left_join(mt_base, tcr_wide)
tcr_wide_mt$orig.ident <- NULL
tcr_wide_mt$sample_id <- NULL
rownames(tcr_wide_mt) <- tcr_wide_mt$barcode
tcr_wide_mt$barcode <- NULL


seurat<-seurat |> AddMetaData(tcr_wide_mt)
seurat<-seurat |> AddMetaData(bcr_wide_mt)

saveRDS(seurat, glue("{seurat_output_dir}/seurat_SCT_mapped_adaptive_added_04_03.rds"))