####################################################################
# Charles Zhang
# Gate Lab
# Northwestern University
# 02/15/2023
####################################################################
# Format contig and clonotype information and add to meta data in seurat object
# Expected output:
#     filtered_{bcr/tcr}_contigs_merged: contig information for each cell
#     {bcr/tcr}_clonotypes.rds: clonotype information for each clonotype
####################################################################
# INSTRUCTIONS FOR USE:
# Change bcr_output_dir to the root directory holding all samples' 01.cellranger bcr output
# Change tcr_output_dir to the root directory holding all samples' 01.cellranger tcr output
# Change seurat_obj to input seurat object produced after general QC
# Change contig_clonotype_output to output directory for all clonotype and contig table output
# Change seurat_output_dir to where the seurat object with added meta data is stored
####################################################################
source("../00.ref/config/immune_profiling_config.R")
bcr_output_dir <- glue("{out_dir_root}/01.cellranger_bcr_vdj")
tcr_output_dir <- glue("{out_dir_root}/01.cellranger_tcr_vdj")
seurat_obj <- glue("{out_dir_root}/03.quality_control/seurat_merged_02.rds")
contig_clonotype_output <- glue("{out_dir_root}/04.contig_clonotype_output")
seurat_output_dir <- glue("{out_dir_root}/04.seurat")
dir.create(seurat_output_dir, showWarnings = FALSE, recursive = TRUE )


# separate out Alpha and Beta chain from the cdr3s_aa (amino acid sequence) into separate columns
# remove cells that does not have both alpha and beta chain
process_tcr_cdr3s <- function(clonotypes_merged) {
  clonotypes_merged$tra_cdr3s <- clonotypes_merged$cdr3s_aa |> lapply(get_chain, "TRA:") |> as.character()
  clonotypes_merged$trb_cdr3s <- clonotypes_merged$cdr3s_aa |> lapply(get_chain, "TRB:") |> as.character()
  clonotypes_merged[clonotypes_merged == ""] <- NA
  clonotypes_merged[!is.na(clonotypes_merged$tra_cdr3s) & !is.na(clonotypes_merged$trb_cdr3s),]
}

# separate out Light and Heavy Chain from the cdr3s_aa (amino acid sequence) into separate columns
# a separate column for different kind of light or heavy chain is created
# remove cells that does not have both a light and a heavy chain
# remove cells that has both L and K light chain
process_bcr_cdr3s <- function(clonotypes_merged) {
  igl <- clonotypes_merged$cdr3s_aa |> lapply(get_chain, "IGL:") |> as.character()
  igk <- clonotypes_merged$cdr3s_aa |> lapply(get_chain, "IGK:") |> as.character()
  light_chain <- seq_along(igl) |> lapply(\(i) {
    if (igk[[i]] == "" & igl[[i]] != "") { # if only L chain is present
      list(igl_crd3s = igl[[i]], light_chain = "L")
    } else if (igk[[i]] != "" & igl[[i]] == "" ) { # if only K chain is present
      list(igl_crd3s = igk[[i]], light_chain = "K")
    } else {
      list(igl_crd3s = NA, light_chain = NA)
    }
  }) |> bind_rows()
  clonotypes_merged$igl_cdr3s <- light_chain$igl_crd3s
  clonotypes_merged$light_chain <- light_chain$light_chain
  clonotypes_merged$igh_cdr3s <- clonotypes_merged$cdr3s_aa |> lapply(get_chain, "IGH:") |> as.character()
  clonotypes_merged[clonotypes_merged == ""] <- NA
  clonotypes_merged[!is.na(clonotypes_merged$igl_cdr3s) & !is.na(clonotypes_merged$igh_cdr3s),]
}

get_chain <- function(value, pattern) {
  sequences <- value |> str_split(';')
  grep(
    pattern,
    sequences[[1]],
    value = TRUE
  ) |> paste(collapse = ";") %>%
    gsub(pattern, "", .) |>
      as.character()
}

# Filtered contig information processing
samples <- tcr_output_dir |> list.files()
all_tcr_filtered_contigs <- "{tcr_output_dir}/{samples}/outs/filtered_contig_annotations.csv" |> glue()

tcr_df_list <- lapply(all_tcr_filtered_contigs, function(cur_file){
  cur_df <- read.csv(cur_file)
  print(cur_file)
  id <- unlist(strsplit(cur_file[[1]], "/")) |>
    tail(3) |>
    pluck(1)
  print(id)
  cur_df[["id"]] <- id
  cur_df
})

filtered_tcr_contigs_merged <- rbindlist(tcr_df_list)
print(nrow(filtered_tcr_contigs_merged))
filtered_tcr_contigs_merged <- filtered_tcr_contigs_merged[!is.na(filtered_tcr_contigs_merged$raw_clonotype_id)]
print(nrow(filtered_tcr_contigs_merged))
filtered_tcr_contigs_merged[["barcode"]]<- "{filtered_tcr_contigs_merged$id}_{filtered_tcr_contigs_merged$barcode}" |>
  glue()
filtered_tcr_contigs_merged[["raw_clonotype_id"]]<- "{filtered_tcr_contigs_merged$id}_{filtered_tcr_contigs_merged$raw_clonotype_id}" |>
  glue()


# merge bcr filtered contigs
samples <- bcr_output_dir |> list.files()
all_bcr_filtered_contigs <- "{bcr_output_dir}/{samples}/outs/filtered_contig_annotations.csv" |> glue()

bcr_df_list <- lapply(all_bcr_filtered_contigs, function(cur_file){
  cur_df <- read.csv(cur_file)
  id <- unlist(strsplit(cur_file[[1]], "/")) |>
    tail(3) |>
    pluck(1)
  cur_df[["id"]] <- id
  cur_df
})


filtered_bcr_contigs_merged <- rbindlist(bcr_df_list)
filtered_bcr_contigs_merged <- filtered_bcr_contigs_merged[!is.na(filtered_bcr_contigs_merged$raw_clonotype_id)]
filtered_bcr_contigs_merged[["barcode"]]<- "{filtered_bcr_contigs_merged$id}_{filtered_bcr_contigs_merged$barcode}" |>
  glue()
filtered_bcr_contigs_merged[["raw_clonotype_id"]]<- "{filtered_bcr_contigs_merged$id}_{filtered_bcr_contigs_merged$raw_clonotype_id}" |>
  glue()


# Clonotypes information Processing
# merge tcr clonotypes

samples <- tcr_output_dir |> list.files()
all_tcr_clonotypes <- "{tcr_output_dir}/{samples}/outs/clonotypes.csv" |> glue()
tcr_df_list <- lapply(all_tcr_clonotypes, function(cur_file){
  cur_df <- read.csv(cur_file)
  id <- unlist(strsplit(cur_file[[1]], "/")) |>
    tail(3) |>
    pluck(1)
  cur_df[["id"]] <-id
  cur_df
})
clonotypes_tcr_merged <- rbindlist(tcr_df_list)|>
  process_tcr_cdr3s()

# merge bcr clonotypes
samples <- bcr_output_dir |> list.files()
all_bcr_clonotypes <- "{bcr_output_dir}/{samples}/outs/clonotypes.csv" |> glue()

bcr_df_list <- lapply(all_bcr_clonotypes, function(cur_file){
  cur_df <- read.csv(cur_file)
  id <- unlist(strsplit(cur_file[[1]], "/")) |>
    tail(3) |>
    pluck(1)
  cur_df[["id"]] <-id
  cur_df
})
clonotypes_bcr_merged <- rbindlist(bcr_df_list)|>
  process_bcr_cdr3s()

clonotypes_tcr_merged[["raw_clonotype_id"]]<-glue("{clonotypes_tcr_merged$id}_{clonotypes_tcr_merged$clonotype_id}")
clonotypes_bcr_merged[["raw_clonotype_id"]] <- glue("{clonotypes_bcr_merged$id}_{clonotypes_bcr_merged$clonotype_id}")

# Save output
write.csv(filtered_bcr_contigs_merged, file = glue("{contig_clonotype_output}/filtered_bcr_contigs_merged.csv"), quote = F, row.names = F)
write.csv(filtered_tcr_contigs_merged, file = glue("{contig_clonotype_output}/filtered_tcr_contigs_merged.csv"), quote = F, row.names = F)
write.csv(clonotypes_bcr_merged, file = glue("{contig_clonotype_output}/bcr_clonotypes.csv"), quote = F, row.names = F)
write.csv(clonotypes_tcr_merged, file = glue("{contig_clonotype_output}/tcr_clonotypes.csv"), quote = F, row.names = F)



# add meta data
unique_tcr_clonotype_barcode_pair <- filtered_tcr_contigs_merged[,c("barcode", "raw_clonotype_id")]|>
  unique()

# filter for cells to have exactly one beta and alpha chain
# since we already filtered for both having alpha and beta chain, any cell ID that has more than two rows
# in this filtered contig dataframe has more than one alpha or beta chain
all_tcr_cell_id_freq <- table(filtered_tcr_contigs_merged$barcode)|>as.data.frame()
print(glue("INFO: TCR entries before filtering for dual chain barcodes {nrow(filtered_tcr_contigs_merged)}"))
all_tcr_cell_id_freq_dual_only <- all_tcr_cell_id_freq|>
  dplyr::filter(Freq==2)
print(glue("INFO: T cells before filtering for dual chain barcodes {length(filtered_tcr_contigs_merged$barcode|>unique())}"))

single_tra_trb_cells <- all_tcr_cell_id_freq_dual_only[,1]|>as.character()
filtered_tcr_contigs_merged_single_tra_trb <- subset(filtered_tcr_contigs_merged, barcode %in% single_tra_trb_cells)

print(glue("INFO: TCR entries after filtering for dual chain barcodes {nrow(filtered_tcr_contigs_merged_single_tra_trb)}"))
print(glue("INFO: T cells after filtering for dual chain barcodes {length(filtered_tcr_contigs_merged_single_tra_trb$barcode|>unique())}"))

names(clonotypes_tcr_merged)[names(clonotypes_tcr_merged) == "frequency"] <- "tcr_frequency"
tcr_columns <- c('raw_clonotype_id', 'tcr_frequency', 'tra_cdr3s', 'trb_cdr3s', 'inkt_evidence', 'mait_evidence')
clonotypes_tcr_merged <- dplyr::select(clonotypes_tcr_merged, tcr_columns)

# left join to filtered contigs
all_tcr_metadata <- merge(x = filtered_tcr_contigs_merged_single_tra_trb, y = clonotypes_tcr_merged, all.x=T)
names(all_tcr_metadata)[names(all_tcr_metadata) == "raw_clonotype_id"] <- "tcr_clonotype_id"
names(all_tcr_metadata)[names(all_tcr_metadata) == "raw_consensus_id"] <- "tcr_consensus_id"


##########################################################################################
# check if the dual chain is one beta and one alpha
# upon checking, all the dual ones have one alpha and one beta
# keep code, but commentted for faster runtime
##########################################################################################
# tcr_barcodes_with_only_single_TRB_A_chain<-pbsapply(all_tcr_cell_id_freq_dual_only[,1], function(cur_barcode){
#   cur_df <- dplyr::filter(filtered_tcr_contigs_merged, barcode == cur_barcode)
#   if(("TRB" %in% cur_df$chain) && ("TRA" %in% cur_df$chain)){
#     out <- cur_barcode
#   }else{
#     out <- NA
#   }
#   as.character(out)
# })
##########################################################################################


unique_bcr_clonotype_barcode_pair <- filtered_bcr_contigs_merged[,c("barcode", "raw_clonotype_id")]|>
  unique()

# filter for cells to have exactly one light and one heavy chain
all_bcr_cell_id_freq <- table(filtered_bcr_contigs_merged$barcode)|>as.data.frame()
print(glue("INFO: BCR entries before filtering for dual chain barcodes {nrow(filtered_bcr_contigs_merged)}"))


all_bcr_cell_id_freq_dual_only <- all_bcr_cell_id_freq|>
  dplyr::filter(Freq==2)
print(glue("INFO: B cells before filtering for dual chain barcodes {length(filtered_bcr_contigs_merged$barcode|>unique())}"))


single_heavy_light_cells <- all_bcr_cell_id_freq_dual_only[,1]|>as.character()
filtered_bcr_contigs_merged_single_heavy_light <- subset(filtered_bcr_contigs_merged, barcode %in% single_heavy_light_cells)

print(glue("INFO: BCR entries after filtering for dual chain barcodes {nrow(filtered_bcr_contigs_merged_single_heavy_light)}"))
print(glue("INFO: B cells after filtering for dual chain barcodes {length(filtered_bcr_contigs_merged_single_heavy_light$barcode|>unique())}"))

names(clonotypes_bcr_merged)[names(clonotypes_bcr_merged) == "frequency"] <- "bcr_frequency"
bcr_columns <- c('raw_clonotype_id', 'bcr_frequency', 'igl_cdr3s', 'igh_cdr3s', 'light_chain')
clonotypes_bcr_merged <- dplyr::select(clonotypes_bcr_merged, bcr_columns)
# left join to filtered contigs
all_bcr_metadata <- merge(x = filtered_bcr_contigs_merged_single_heavy_light, y = clonotypes_bcr_merged, all.x=T)
names(all_bcr_metadata)[names(all_bcr_metadata) == "raw_clonotype_id"] <- "bcr_clonotype_id"
names(all_bcr_metadata)[names(all_bcr_metadata) == "raw_consensus_id"] <- "bcr_consensus_id"

##########################################################################################
# check if the dual chain is one beta and one alpha
# upon checking, all the dual ones have one light and one heavy
# keep code, but commentted for faster runtime
##########################################################################################
# bcr_barcodes_with_only_single_heavy_light_chain<-pbsapply(all_bcr_cell_id_freq_dual_only[,1], function(cur_barcode){
#   cur_df <- dplyr::filter(filtered_bcr_contigs_merged, barcode == cur_barcode)
#   if(("IGL" %in% cur_df$chain) && ("IGH" %in% cur_df$chain)){
#     out <- cur_barcode
#   }else if(("IGK" %in% cur_df$chain) && ("IGH" %in% cur_df$chain)){
#     out <- cur_barcode
#   }else{
#     out <- NA
#   }
#   as.character(out)
# })
##########################################################################################

write.csv(all_bcr_metadata, file = glue("{contig_clonotype_output}/all_bcr_meta.csv"), quote = F, row.names = F)
write.csv(all_tcr_metadata, file = glue("{contig_clonotype_output}/all_tcr_meta.csv"), quote = F, row.names = F)

# add to seurat
seurat <- readRDS(seurat_obj)
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


saveRDS(seurat, glue("{seurat_output_dir}/seurat_04_01.rds"))