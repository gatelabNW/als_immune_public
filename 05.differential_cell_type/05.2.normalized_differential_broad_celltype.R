# Charles Zhang
# Gate Lab
# Northwestern University
# 09/28/2023
####################################################################
# Analyze if there is enrichment of immune cell types
####################################################################

# source("../00.ref/config/immune_profiling_config.R")
# seurat_output_dir <- glue("{out_dir_root}/04.seurat")
# input_seurat <- glue("{seurat_output_dir}/seurat_SCT_mapped_04_03.rds")
# output_root <- glue("{out_dir_root}/05.differential_cell_type")
# dir.create(output_root, showWarnings = FALSE, recursive = TRUE)
#
#
# s <- readRDS(input_seurat)
#
# mt<-s@meta.data
# sample_celltype_tbl <- dplyr::select(mt, orig.ident, predicted.celltype.l1, diagnosis_general)
# sample_celltype_stats <- sample_celltype_tbl |>
#   group_by(orig.ident, predicted.celltype.l1) |>
#   summarise(count = n()) |>
#   left_join(sample_celltype_tbl|>
#   group_by(orig.ident) |>
#   summarise(total = n())) |>
#   mutate(percentage = (count / total) * 100) |>
#   select(orig.ident, predicted.celltype.l1, percentage)|>
#   as.data.frame()
#
# ref <- "../00.ref/meta/samples-metadata.csv"|>read.csv()
# gex_idx2diagnosis_general <- ref$diagnosis_general
# names(gex_idx2diagnosis_general) <- ref$gex_index
# sample_celltype_stats[["diagnosis_general"]] <- sapply(sample_celltype_stats$orig.ident,
# function(x) gex_idx2diagnosis_general[[x]])
#
# gex_idx2diagnosis <- ref$diagnosis
# names(gex_idx2diagnosis) <- ref$gex_index
# sample_celltype_stats[["diagnosis"]] <- sapply(sample_celltype_stats$orig.ident,
# function(x) gex_idx2diagnosis[[x]])
#
# p <- ggplot(sample_celltype_stats, aes(x = diagnosis_general, y = percentage, fill = diagnosis_general)) +
#   geom_boxplot() +
#   facet_wrap(~ predicted.celltype.l1, scales = "free") +
#   labs(title = "Diagnosis General Normalized Cell Type Proportion") +
#   theme_minimal()+
#   ggpubr::stat_compare_means(method = "wilcox.test")
# out_png <- glue("{output_root}/diagnosis_general_ct_lv1.png")
# png(out_png, width =1920, height = 1920)
# print(p)
# dev.off()
#
# p <- ggplot(sample_celltype_stats, aes(x = diagnosis, y = percentage, fill = diagnosis)) +
#   geom_boxplot() +
#   facet_wrap(~ predicted.celltype.l1, scales = "free") +
#   labs(title = "Diagnosis Normalized Cell Type Proportion") +
#   theme_minimal()+
#   ggpubr::stat_compare_means(method = "anova")
#
# out_png <- glue("{output_root}/diagnosis_ct_lv1.png")
# png(out_png, width =1920, height = 1920)
# print(p)
# dev.off()
#
# # save the calculated per sample normalized cell type data
# out_csv <- glue("{output_root}/normalized_celltype_count_lv1.csv")
# write.csv(sample_celltype_stats, file = out_csv, quote = FALSE, row.names = FALSE)


# CRISPR_clean
source("../00.ref/config/CRISPR_clean_config.R")
seurat_output_dir <- glue("{out_dir_root}/04.seurat")
input_seurat <- glue("{seurat_output_dir}/seurat_04_05.rds")
output_root <- glue("{out_dir_root}/05.differential_cell_type")
dir.create(output_root, showWarnings = FALSE, recursive = TRUE)


s <- readRDS(input_seurat)

mt<-s@meta.data
sample_celltype_tbl <- dplyr::select(mt, orig.ident, predicted.celltype.l1, diagnosis_general)
sample_celltype_stats <- sample_celltype_tbl |>
  group_by(orig.ident, predicted.celltype.l1) |>
  summarise(count = n()) |>
  left_join(sample_celltype_tbl|>
              group_by(orig.ident) |>
              summarise(total = n())) |>
  mutate(percentage = (count / total) * 100) |>
  select(orig.ident, predicted.celltype.l1, percentage)|>
  as.data.frame()

ref <- "../00.ref/meta/samples-metadata.csv"|>read.csv()
gex_idx2diagnosis_general <- ref$diagnosis_general
names(gex_idx2diagnosis_general) <- ref$gex_index
sample_celltype_stats[["diagnosis_general"]] <- sapply(sample_celltype_stats$orig.ident,
                                                       function(x) gex_idx2diagnosis_general[[x]])

gex_idx2diagnosis <- ref$diagnosis
names(gex_idx2diagnosis) <- ref$gex_index
sample_celltype_stats[["diagnosis"]] <- sapply(sample_celltype_stats$orig.ident,
                                               function(x) gex_idx2diagnosis[[x]])

p <- ggplot(sample_celltype_stats, aes(x = diagnosis_general, y = percentage, fill = diagnosis_general)) +
  geom_boxplot() +
  facet_wrap(~ predicted.celltype.l1, scales = "free") +
  labs(title = "Diagnosis General Normalized Cell Type Proportion") +
  theme_minimal()+
  ggpubr::stat_compare_means(method = "wilcox.test")
out_png <- glue("{output_root}/diagnosis_general_ct_lv1.png")
png(out_png, width =1920, height = 1920)
print(p)
dev.off()

p <- ggplot(sample_celltype_stats, aes(x = diagnosis, y = percentage, fill = diagnosis)) +
  geom_boxplot() +
  facet_wrap(~ predicted.celltype.l1, scales = "free") +
  labs(title = "Diagnosis Normalized Cell Type Proportion") +
  theme_minimal()+
  ggpubr::stat_compare_means(method = "anova")

out_png <- glue("{output_root}/diagnosis_ct_lv1.png")
png(out_png, width =1920, height = 1920)
print(p)
dev.off()

# save the calculated per sample normalized cell type data
out_csv <- glue("{output_root}/normalized_celltype_count_lv1.csv")
write.csv(sample_celltype_stats, file = out_csv, quote = FALSE, row.names = FALSE)

# Immune enriched
# source("../00.ref/config/immune_panel_config.R")
# seurat_output_dir <- glue("{out_dir_root}/04.seurat")
# input_seurat <- glue("{seurat_output_dir}/seurat_SCT_mapped_04_02.rds")
# output_root <- glue("{out_dir_root}/05.differential_cell_type")
# dir.create(output_root, showWarnings = FALSE, recursive = TRUE)
#
#
# s <- readRDS(input_seurat)
#
# mt<-s@meta.data
# sample_celltype_tbl <- dplyr::select(mt, orig.ident, predicted.celltype.l1, diagnosis_general)
# sample_celltype_stats <- sample_celltype_tbl |>
#   group_by(orig.ident, predicted.celltype.l1) |>
#   summarise(count = n()) |>
#   left_join(sample_celltype_tbl|>
#               group_by(orig.ident) |>
#               summarise(total = n())) |>
#   mutate(percentage = (count / total) * 100) |>
#   select(orig.ident, predicted.celltype.l1, percentage)|>
#   as.data.frame()
#
# ref <- "../00.ref/meta/samples-metadata.csv"|>read.csv()
# gex_idx2diagnosis_general <- ref$diagnosis_general
# names(gex_idx2diagnosis_general) <- ref$gex_index
# sample_celltype_stats[["diagnosis_general"]] <- sapply(sample_celltype_stats$orig.ident,
#                                                        function(x) gex_idx2diagnosis_general[[x]])
#
# gex_idx2diagnosis <- ref$diagnosis
# names(gex_idx2diagnosis) <- ref$gex_index
# sample_celltype_stats[["diagnosis"]] <- sapply(sample_celltype_stats$orig.ident,
#                                                function(x) gex_idx2diagnosis[[x]])
#
# p <- ggplot(sample_celltype_stats, aes(x = diagnosis_general, y = percentage, fill = diagnosis_general)) +
#   geom_boxplot() +
#   facet_wrap(~ predicted.celltype.l1, scales = "free") +
#   labs(title = "Diagnosis General Normalized Cell Type Proportion") +
#   theme_minimal()+
#   ggpubr::stat_compare_means(method = "wilcox.test")
# out_png <- glue("{output_root}/diagnosis_general_ct_lv1.png")
# png(out_png, width =1920, height = 1920)
# print(p)
# dev.off()
#
# p <- ggplot(sample_celltype_stats, aes(x = diagnosis, y = percentage, fill = diagnosis)) +
#   geom_boxplot() +
#   facet_wrap(~ predicted.celltype.l1, scales = "free") +
#   labs(title = "Diagnosis Normalized Cell Type Proportion") +
#   theme_minimal()+
#   ggpubr::stat_compare_means(method = "anova")
#
# out_png <- glue("{output_root}/diagnosis_ct_lv1.png")
# png(out_png, width =1920, height = 1920)
# print(p)
# dev.off()
#
# # save the calculated per sample normalized cell type data
# out_csv <- glue("{output_root}/normalized_celltype_count_lv1.csv")
# write.csv(sample_celltype_stats, file = out_csv, quote = FALSE, row.names = FALSE)

