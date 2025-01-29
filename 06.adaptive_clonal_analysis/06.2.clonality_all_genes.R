source("../00.ref/config/immune_profiling_config.R")
adaptive_in_dir <- "/projects/b1042/Gate_Lab/projects/als-project/TCRB_adaptive_seq/00.data"
adaptive_out_dir_data <- "/projects/b1042/Gate_Lab/projects/als-project/TCRB_adaptive_seq/02.expansion"
DE_dir <- "/projects/b1042/Gate_Lab/projects/als-project/crispr_clean_final/05.DEG/degs/diagnosis_general/SCT/age_sex"
BH_thres <- 0.01
lfc_thres <- 0.585

library(ggnetwork)
library(network)
library(ggraph)
library(igraph)
library(RColorBrewer)
library(e1071)

dir.create(adaptive_out_dir_data, showWarnings = FALSE, recursive = TRUE)
filtered_table <- glue("{adaptive_in_dir}/filtered_rearrangements.tsv")
filtered_table <- read.table(filtered_table, sep = "\t", header = T)
clonal_expansion_cutoff <- 2
filtered_table <- dplyr::filter(filtered_table, templates > clonal_expansion_cutoff)



# Get all CRISPR clean upregulated TCR-B genes
all_DE_tables <- list.files(DE_dir, full.names = T)
TRB_DE_table <- list()

for(cur_file in all_DE_tables){
  cur_DE_table <- read.csv(cur_file)
  cur_DE <- dplyr::filter(cur_DE_table, BH < BH_thres & avg_log2FC > lfc_thres)
  cur_DE_TRB <- cur_DE[grepl("TRB", cur_DE$X), ]
  if(nrow(cur_DE_TRB)>0){
    cur_ct <- sub(".*___(.*)__.*", "\\1", cur_file)
    cur_DE_TRB[["cell_type"]] <- cur_ct
    TRB_DE_table[[cur_file]] <- cur_DE_TRB
  }
}

TRB_DE_table <- data.table::rbindlist(TRB_DE_table)|>
  as.data.frame()



# Examine those TRB for more clonal expansion in ALS patients
all_TRB <- TRB_DE_table$X
all_combs <- c("J_10x", "V_10x", "D_10x")
sk_df_ls <- list()
for(cur_comb in all_combs){
  cur_comb_TR <- filtered_table[[cur_comb]]|>unique()
  for(cur_TR in cur_comb_TR){
    cur_TRB_adaptive_table <- filtered_table|>dplyr::filter(.data[[cur_comb]] == cur_TR)
    print(glue("Working on {cur_TR}"))
    if(nrow(cur_TRB_adaptive_table) < 100){
      next
    }
    cur_TRB_adaptive_table_hc <- cur_TRB_adaptive_table[grepl("H", cur_TRB_adaptive_table$sample_name), ]
    cur_TRB_adaptive_table_als <- cur_TRB_adaptive_table[!grepl("H", cur_TRB_adaptive_table$sample_name), ]

    # calculate skweness
    als_sk_norm <- cur_TRB_adaptive_table_als$productive_frequency|>skewness(type = 2)
    hc_sk_norm <- cur_TRB_adaptive_table_hc$productive_frequency|>skewness(type = 2)
    als_sk_raw <- cur_TRB_adaptive_table_als$total_templates|>skewness(type = 2)
    hc_sk_raw <- cur_TRB_adaptive_table_hc$total_templates|>skewness(type = 2)
    if(cur_TR %in% all_TRB){
      DE <- "DE"
    }else{
      DE <- "Not DE"
    }
    als_hc_diff_norm <- als_sk_norm - hc_sk_norm
    als_hc_diff_raw <- als_sk_raw - hc_sk_raw
    # if(als_hc_diff > 20){
    #   color = "red"
    # }else if(als_hc_diff < -20){
    #   color = "blue"
    # }else{
    #   color = "gray50"
    # }

    output<-list(
      "hc_score" = hc_sk,
      "als_score" = als_sk,
      "DE" = DE,
      "gene" = cur_TR,
      "diff_norm" = als_hc_diff_norm,
      "diff_raw" = als_hc_diff_raw
    )
    sk_df_ls[[cur_TR]] <- output

    if(als_hc_diff_raw > 20 || als_hc_diff_norm > 20){
      result_hc <- cur_TRB_adaptive_table_hc |>
        group_by(amino_acid) |>
        dplyr::summarise(
          total_templates = sum(templates),
          combined_sample_name = paste(sample_name, collapse = ", ")
        )
      result_hc[["condition"]] <- "HC"
      result_als <- cur_TRB_adaptive_table_als |>
        group_by(amino_acid) |>
        dplyr::summarise(
          total_templates = sum(templates),
          combined_sample_name = paste(sample_name, collapse = ", ")
        )
      result_als[["condition"]] <- "als"

      # down sample smaller templates count for plotting efficiency
      if(nrow(result_als) >1500){
        result_als_low <- result_als |>
          filter(total_templates < 100) |>
          sample_n(min(1000, n()))  # Sample 1000 rows or less if there are fewer than 1000

        # Step 2: Select rows where templates >= 100
        result_als_high <- result_als |>
          filter(total_templates >= 100)
      }
      result_als <- bind_rows(result_als_low, result_als_high)

      if(nrow(result_hc) >1500){
        result_hc_low <- result_hc |>
          filter(total_templates < 100) |>
          sample_n(min(1000, n()))  # Sample 1000 rows or less if there are fewer than 1000

        # Step 2: Select rows where templates >= 100
        result_hc_high <- result_hc |>
          filter(total_templates >= 100)
      }
      result_hc <- bind_rows(result_hc_low, result_hc_high)

      result <- rbind(result_hc, result_als)


      result$num_strings <- sapply(strsplit(as.character(result$combined_sample_name), ","), length)

      # Create a network matrix
      network_matrix <- matrix(0, nrow = nrow(result), ncol = nrow(result))

      # Create a network
      net <- network(network_matrix, directed = FALSE)

      # Create a layout, adjusting the parameters to make the layout more circular
      net_layout <- ggnetwork(net, layout = 'fruchtermanreingold', arrow.gap = 0, niter = 5000, area = nrow(result)^2*1000)

      # Merge layout with your data
      df_net <- merge(net_layout, result, by.x = "vertex.names", by.y = "row.names")
      df_net[["condition"]] <- factor(df_net$condition, levels = c("HC", "als"))

      # Calculate quantiles and max value
      quantiles <- round(quantile(df_net$total_templates, probs = c(0, 0.9, 0.99, 0.999)))
      max_value <- round(max(df_net$total_templates))
      top_percentile <- round(quantile(df_net$total_templates, probs = 0.9999))
      breaks <- c(-Inf, quantiles, top_percentile, max_value + 1)
      labels <- c()
      for (i in 1:(length(breaks) - 1)) {
        if (is.infinite(breaks[i])) {
          # Format for -Inf
          labels <- c(labels, paste("<", as.integer(format(breaks[i + 1]))))
        } else if (is.infinite(breaks[i + 1])) {
          # Format for max_value + 1
          labels <- c(labels, paste(">", as.integer(format(breaks[i]))))
        } else {
          labels <- c(labels, paste(as.integer(format(breaks[i])), "-", as.integer(format(breaks[i + 1]))))
        }
      }
      df_net$size_category <- cut(df_net$total_templates, breaks = breaks, labels = labels)

      # Define sizes for each label (you may need to adjust these sizes)
      size_values <- setNames(c(1, 2, 3, 5, 15, 25), labels)

      # Plot
      p <- ggplot(data = df_net, aes(x = x, y = y, size = size_category, color = num_strings)) +
        geom_point(alpha=0.7) +
        scale_size_manual(values = size_values) +  # Use size_values here
        scale_color_gradient(low = "blue", high = "red") +
        facet_wrap(~condition, nrow = 1) +
        theme_void() +
        labs(color = "Number of Unique Recombinations", size = "Clonality Category") +
        theme(legend.position = "right")

      # Output the plot
      cur_plt_out <- glue("{adaptive_out_dir_data}/plots/{cur_TR}.pdf")
      pdf(file = cur_plt_out, height = 9, width = 15)
      print(p)  # Use print() to properly display ggplot plots in a graphics device
      dev.off()

      cur_plot_data_out <- glue("{adaptive_out_dir_data}/data/{cur_TR}_net.csv")
      write.csv(df_net, file = cur_plot_data_out, row.names = F, quote = F)
    }
  }
}

# sk table output
sk_data_out <- glue("{adaptive_out_dir_data}/data/!sk.csv")
sk_df <- rbindlist(sk_df_ls)|>as.data.frame()
write.csv(sk_df, file = sk_data_out, row.names = F, quote = F)