library(fgsea)
library(msigdbr)
library(tidyverse)
library(glue)
source("../00.ref/config/CRISPR_clean_config.R")

root_dir <- "/projects/b1042/Gate_Lab/projects/als-project/crispr_clean_final/05.DEG_717_filter/degs/female_c9_hc/SCT/age"

# root_dir <- "/projects/b1042/Gate_Lab/projects/als-project/crispr_clean_final/05.DEG_717_filter/degs/sALS_hc/SCT/age_sex"
# all_deg_files <- list.files(root_dir, full.names = T, pattern = "*control.csv$")

root_dir <- "/projects/b1042/Gate_Lab/projects/als-project/crispr_clean_final/05.DEG_717_filter/degs/diagnosis/SCT/age_sex"
all_deg_files <- list.files(root_dir, full.names = T, pattern = "als_fast_vs._healthy_control\\.csv$|als_slow_vs._healthy_control\\.csv$")


# root_dir <- "/projects/b1042/Gate_Lab/projects/als-project/crispr_clean_final/05.DEG_717_filter/degs/diagnosis_general/SCT/age_sex"
# all_deg_files <- list.files(root_dir, full.names = T, pattern = "*control.csv$")

category <- "H"
pathwaysDF <- msigdbr("human", category = category)
pathways <- split(pathwaysDF$gene_symbol, pathwaysDF$gs_name)


for(cur_deg_file in all_deg_files){
  cur_df <- read.csv(cur_deg_file)

  if(nrow(cur_df) == 0){
    next
  }

  na_rows <- cur_df[is.na(cur_df$BH), ]
  if(na_rows|>nrow()>0){
    print(na_rows)
  }
  cur_df <- cur_df[!is.na(cur_df$BH), ]



  cur_df[["sig"]] <- log10(cur_df[["BH"]]) * -1 * sign(cur_df[["avg_log2FC"]])
  max_finite_value <- max(cur_df$sig[cur_df$sig != Inf], na.rm = TRUE)
  cur_df$sig[cur_df$sig == Inf] <- max_finite_value + 1
  min_finite_value <- min(cur_df$sig[cur_df$sig != -Inf], na.rm = TRUE)
  cur_df$sig[cur_df$sig == -Inf] <- min_finite_value - 1
  cur_df <- cur_df[order(cur_df$sig, decreasing = F),]
  cur_cluster_rank <- cur_df$sig
  names(cur_cluster_rank) <- cur_df$gene_id


  fgseaRes <- fgsea(pathways = pathways,
                    stats    = cur_cluster_rank,
                    minSize  = 15,
                    maxSize  = 500)
  cur_fgsea_res_file <- str_replace(cur_deg_file, ".csv", "_fgsea.csv")
  fwrite(fgseaRes, file=cur_fgsea_res_file, sep="\t", sep2=c("", " ", ""))

  topPathwaysUp <- fgseaRes[ES > 0][head(order(pval), n=10), pathway]
  # topPathwaysDown <- fgseaRes[ES < 0][head(order(pval), n=10), pathway]
  # topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
  p <- plotGseaTable(pathways[topPathwaysUp], cur_cluster_rank, fgseaRes,
                     gseaParam=0.5)

  cur_out_file <- str_replace(cur_deg_file, ".csv", glue("__{category}_fgsea.pdf"))
  cur_out_file <- str_replace(cur_out_file, "degs", "plots")
  pdf(cur_out_file, width = 21)
  print(p)
  dev.off()
}