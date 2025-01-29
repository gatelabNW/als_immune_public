####################################################################
# Charles Zhang
# Gate Lab
# Northwestern University
# 02/13/2023
####################################################################
# Format clonotypes and contig information from 10X vdj output
# Expected output:
#     *contigs_merged.csv: a csv file that has the merged contig information, one for tcr and one for bcr
#     *clonotypes.csv: a csv file that has the clonoype information, one for tcr and one for bcr
####################################################################
# INSTRUCTIONS FOR USE:
# Change bcr_output_dir, tcr_output_dir to the root output directory for tcr/bcr cell ranger
# Change contig_clonotype_output to a directory for the merged clonotypes and contigs files
####################################################################
source("../00.ref/config/config.R")

# Input directory for cell ranger bcr output and cell ranger tcr output
bcr_output_dir <- glue("{out_dir_root}/01.cellranger_bcr_vdj")
tcr_output_dir <- glue("{out_dir_root}/01.cellranger_tcr_vdj")

# Output directory for merged contigs and clonotypes
contig_clonotype_output <- glue("{out_dir_root}/04.contig_clonotype_output")
dir.create(contig_clonotype_output, showWarnings = FALSE, recursive = TRUE )



process_tcr_cdr3s <- function(clonotypes_merged) {
  clonotypes_merged$tra_cdr3s <- clonotypes_merged$cdr3s_aa |> lapply(get_chain, "TRA:") |> as.character()
  clonotypes_merged$trb_cdr3s <- clonotypes_merged$cdr3s_aa |> lapply(get_chain, "TRB:") |> as.character()
  clonotypes_merged[clonotypes_merged == ""] <- NA
  clonotypes_merged[!is.na(clonotypes_merged$tra_cdr3s) & !is.na(clonotypes_merged$trb_cdr3s),]
}

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

# merge tcr filtered contigs
samples <- tcr_output_dir |> list.files()
all_tcr_filtered_contigs <- "{tcr_output_dir}/{samples}/outs/filtered_contig_annotations.csv" |> glue()

tcr_df_list <- lapply(all_tcr_filtered_contigs, function(cur_file){
  cur_df <- read.csv(cur_file)
  id <- unlist(strsplit(cur_file[[1]], "/")) |>
    tail(3) |>
    pluck(1)
  cur_df[["id"]] <- rep(id, dim(cur_df)[1])
  cur_df
})


filtered_tcr_contigs_merged <- rbindlist(tcr_df_list)
filtered_tcr_contigs_merged <- filtered_tcr_contigs_merged[!is.na(filtered_tcr_contigs_merged$raw_clonotype_id)]
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
  cur_df[["id"]] <- rep(id, dim(cur_df)[1])
  cur_df
})


filtered_bcr_contigs_merged <- rbindlist(bcr_df_list)
filtered_bcr_contigs_merged <- filtered_bcr_contigs_merged[!is.na(filtered_bcr_contigs_merged$raw_clonotype_id)]
filtered_bcr_contigs_merged[["barcode"]]<- "{filtered_bcr_contigs_merged$id}_{filtered_bcr_contigs_merged$barcode}" |>
  glue()
filtered_bcr_contigs_merged[["raw_clonotype_id"]]<- "{filtered_bcr_contigs_merged$id}_{filtered_bcr_contigs_merged$raw_clonotype_id}" |>
  glue()


# merge tcr clonotypes

samples <- tcr_output_dir |> list.files()
all_tcr_clonotypes <- "{tcr_output_dir}/{samples}/outs/clonotypes.csv" |> glue()
tcr_df_list <- lapply(all_tcr_clonotypes, function(cur_file){
  cur_df <- read.csv(cur_file)
  id <- unlist(strsplit(cur_file[[1]], "/")) |>
    tail(3) |>
    pluck(1)
  cur_df[["id"]] <- rep(id, dim(cur_df)[1])
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
  cur_df[["id"]] <- rep(id, dim(cur_df)[1])
  cur_df
})
clonotypes_bcr_merged <- rbindlist(bcr_df_list)|>
  process_bcr_cdr3s()


write.csv(filtered_bcr_contigs_merged, file = glue("{contig_clonotype_output}/filtered_bcr_contigs_merged.csv", quote = F, row.names = F))
write.csv(filtered_tcr_contigs_merged, file = glue("{contig_clonotype_output}/filtered_tcr_contigs_merged.csv", quote = F, row.names = F))
write.csv(clonotypes_bcr_merged, file = glue("{contig_clonotype_output}/bcr_clonotypes.csv", quote = F, row.names = F))
write.csv(clonotypes_tcr_merged, file = glue("{contig_clonotype_output}/tcr_clonotypes.csv", quote = F, row.names = F))