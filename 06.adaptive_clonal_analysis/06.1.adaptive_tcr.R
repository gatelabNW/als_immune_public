source("../00.ref/config/immune_profiling_config.R")
adaptive_in_dir <- "/projects/b1042/Gate_Lab/projects/als-project/TCRB_adaptive_seq/00.data"
adaptive_out_dir_data <- "/projects/b1042/Gate_Lab/projects/als-project/TCRB_adaptive_seq/01.vdj_comb"
adaptive_DE_genes_file <- "/projects/b1042/Gate_Lab/projects/als-project/scrna_original_final/06.receptors/data/adaptive_DE_table.tsv"
cur_beta_data_out_dir <- glue("{adaptive_out_dir_data}/{target_beta}/data")
cur_beta_plot_out_dir <- glue("{adaptive_out_dir_data}/{target_beta}/plots")
dir.create(cur_beta_data_out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(cur_beta_plot_out_dir, showWarnings = FALSE, recursive = TRUE)

# clonal cutoff
clonal_expansion_cutoff <- 2
clonal_expansion <- T

# format the data
filtered_table <- glue("{adaptive_in_dir}/filtered_rearrangements.tsv")
if (!file.exists(filtered_table)) {
  print("Reading raw and generating filtered table")
  raw_table <- read.table(glue("{adaptive_in_dir}/raw_rearrangements.tsv"), sep = "\t", header = T)

  # remove unproductive and incomplete chains
  filtered_table <- raw_table |>
    filter_all(all_vars(. != "unknown" & !is.na(.) & . != "na"))

  # format the gene nomanclature
  filtered_table[["V_10x"]] <- sapply(filtered_table$v_resolved, function(x){
    result <- sub("TCR", "TR", sub("\\*.*", "", gsub("0(\\d)", "\\1", x)))
    result
  })
  filtered_table[["D_10x"]] <- sapply(filtered_table$d_resolved, function(x){
    result <- sub("TCR", "TR", sub("\\*.*", "", gsub("0(\\d)", "\\1", x)))
    result
  })
  filtered_table[["J_10x"]] <- sapply(filtered_table$j_resolved, function(x){
    result <- sub("TCR", "TR", sub("\\*.*", "", gsub("0(\\d)", "\\1", x)))
    result
  })
  filtered_table[["vdj_recomb"]]<- paste(filtered_table$V_10x, filtered_table$D_10x, filtered_table$J_10x, sep="_")
  write.table(filtered_table, row.names = F, quote = F, sep="\t", file = filtered_table)
} else {
  filtered_table <- read.table(filtered_table, sep = "\t", header = T)
}

# subset only clonally expanded beta chains
if(clonal_expansion){
  filtered_table <- dplyr::filter(filtered_table, templates > clonal_expansion_cutoff)
}

# Get a list of input beta genes
adaptive_DE <- read.table(adaptive_DE_genes_file, sep = "\t", header = T)
all_DE_genes <- adaptive_DE$genes|>paste(collapse = ", ")|>strsplit(", ")
all_DE_genes <- all_DE_genes[[1]]
all_TRB_genes <- unique(all_DE_genes[grepl("TRB",all_DE_genes)])

# Normalize by total clonally expanded TCR count
# Since each expanded TCR is one row, we are how rows each sample has, which is how many unique TRB expanded each sample has
# there appears to be a small portion of rearrangement duplicates in the expanded TRBs, not sure why
# accounts for 0.0024, so shouldn't affect anything
total <- filtered_table$sample_name|>table()|>as.data.frame()
total[] <- lapply(total, function(col) {
  if (is.factor(col)) as.character(col) else col
})
colnames(total) <- c("sample_name", "total_beta_freq")

# for each DE TRB gene
for (target_beta in all_TRB_genes){
  # subset out all the entries for the current beta chain
  target_beta_df <- dplyr::filter(filtered_table, V_10x==target_beta)

  if(nrow(target_beta_df)==0){
    next
  }

  target_beta_freq <- table(dplyr::select(target_beta_df, sample_name, V_10x))|>
    as.data.frame()
  target_beta_freq[] <- lapply(target_beta_freq, function(col) {
    if (is.factor(col)) as.character(col) else col
  })
  target_beta_freq_merged <- merge(target_beta_freq, total)
  target_beta_freq_merged[["normed_vdj_freq"]] <- target_beta_freq_merged$Freq/target_beta_freq_merged$total_beta_freq
  target_beta_freq_merged[["diagnosis"]] <- sapply(target_beta_freq_merged$sample_name, function(x){
    x <- substr(x, start=1, stop=1)
    diagnosis <- switch(x,
                        "C" = "als_c9",
                        "F" = "als_fast",
                        "S" = "als_slow",
                        "H" = "healthy_control")
    diagnosis
  })
  target_beta_freq_merged[["diagnosis_general"]] <- sapply(target_beta_freq_merged$sample_name, function(x){
    x <- substr(x, start=1, stop=1)
    diagnosis <- switch(x,
                        "H" = "healthy_control",
                        "als")
    diagnosis
  })

  diagnosis_general_test <- wilcox.test(normed_vdj_freq ~ diagnosis_general, data = target_beta_freq_merged)
  diagnosis_general_pval <- diagnosis_general_test$p.value
  diagnosis_test <- aov(normed_vdj_freq ~ diagnosis, data = target_beta_freq_merged)
  diagnosis_pval <- summary(diagnosis_test)[[1]]$`Pr(>F)`[1]

  if(diagnosis_pval < 0.1 | diagnosis_general_pval < 0.1){
    print(target_beta)
  }

  cur_beta_out_dir <- glue("{adaptive_out_dir_data}/{target_beta}")
  cur_beta_data_out_dir <- glue("{cur_beta_out_dir}/data")
  cur_beta_plot_out_dir <- glue("{cur_beta_out_dir}/plots")
  dir.create(cur_beta_data_out_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(cur_beta_plot_out_dir, showWarnings = FALSE, recursive = TRUE)
  cur_gene_count_df_file <- glue("{cur_beta_data_out_dir}/{target_beta}.tsv")
  write.table(target_beta_freq_merged, file = cur_gene_count_df_file, quote = F, row.names = F, sep ="\t")

  # plot it
  # round pval
  aov_p <- diagnosis_pval|>round(3)
  rs_p <- diagnosis_general_pval|>round(3)
  gene <- target_beta


  # Create a box plot based on diagnosis
  p1 <- ggplot(target_beta_freq_merged, aes(x = diagnosis_general, y = normed_vdj_freq)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(aes(colour = diagnosis), size = 1.5, width = 0.2) +
    labs(title = glue("{gene}, rank sum {rs_p}"))+
    theme_minimal(base_size = 15)

  p2 <- ggplot(target_beta_freq_merged, aes(x = diagnosis, y = normed_vdj_freq)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(aes(colour = diagnosis), size = 1.5, width = 0.2) +
    labs(title = glue("{gene}, anova {aov_p}")) +
    theme_minimal(base_size = 15)

  out_fig_file <- glue("{cur_beta_plot_out_dir}/{gene}.diagnosis_general.pdf")
  pdf(out_fig_file, width = 9, height = 9)
  plot(p1)
  dev.off()

  out_fig_file <- glue("{cur_beta_plot_out_dir}/{gene}.diagnosis.pdf")
  pdf(out_fig_file, width = 9, height = 9)
  plot(p2)
  dev.off()
}