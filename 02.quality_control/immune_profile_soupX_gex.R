####################################################################
# Charles Zhang
# Gate Lab
# Northwestern University
# 02/07/2023
####################################################################
# Run soupX on 01.cellranger output for ALS project
# Details from: https://rawcdn.githack.com/constantAmateur/SoupX/204b602418df12e9fdb4b68775a8b486c6504fe4/inst/doc/pbmcTutorial.html
# Expected output:
#
####################################################################
# INSTRUCTIONS FOR USE:
# Change cell_ranger_count_dir to the root directory holding all samples' 01.cellranger output
# Change soupx_count_mtx_out_dir to a directory for the new soupX adjusted counts
# change soupx_plot_out_dir to a directory for the soupX plots
####################################################################
source("../00.ref/config/immune_profiling_config.R")


# Input directory of 01.cellranger count output
cell_ranger_count_dir <- glue("{out_dir_root}/01.cellranger_count")

# Output directory for all the clean expression matrix
soupx_count_mtx_out_dir <- glue("{out_dir_root}/02.soupX/count")
dir.create(soupx_count_mtx_out_dir, showWarnings = FALSE, recursive = TRUE )

# Output directory for contamination plots
soupx_plot_out_dir <- glue("{out_dir_root}/02.soupX/plot")
dir.create(soupx_plot_out_dir, showWarnings = FALSE, recursive = TRUE )

# List all the samples
samples <- list.files(cell_ranger_count_dir)
sample_dirs <- "{cell_ranger_count_dir}/{samples}/outs/" |> glue()

contamination_frac <- lapply(sample_dirs, function (cur_dir){
  cur_sample <- strsplit(cur_dir, "/")[[1]][9]
  sc <- load10X(cur_dir)
  sc <- autoEstCont(sc)
  out <- adjustCounts(sc)
  cur_estimated_rho <- sc$fit$rhoEst
  # write adjusted counts
  cur_out_dir <- glue("{soupx_count_mtx_out_dir}/{cur_sample}")
  dir.create(cur_out_dir, showWarnings = FALSE, recursive = TRUE )
  DropletUtils::write10xCounts(cur_out_dir, out, overwrite = TRUE)
  c(cur_sample, cur_estimated_rho)
})

contamination_frac <- contamination_frac|>as.data.frame()|>t()
colnames(contamination_frac) <- c("X", "contamination_frac")
rownames(contamination_frac) <- seq(1:nrow(contamination_frac))
contamination_frac <- as.data.frame(contamination_frac)
contamination_frac$contamination_frac<-contamination_frac$contamination_frac|>as.numeric()
contamination_frac$X <- "NA"

# contamination_frac <- contamination_frac |> bind_rows()
# contamination_frac$contamination_frac <- contamination_frac$contamination_frac |> as.character() |> as.numeric()

contamination_frac |> write.csv(file = "{soupx_plot_out_dir}/contamination_frac.csv" |> glue())
# Plot contamination fraction (box)
# p1 <- ggplot(contamination_frac, aes(x = X, y = contamination_frac)) +
#   geom_boxplot() +
#   geom_jitter(shape=16, position=position_jitter(0.2)) +
#   theme_Publication_blank() +
#   theme(axis.title.x=element_blank(),
#         axis.text.x=element_blank(),
#         axis.ticks.x=element_blank())

# Plot contamination fraction (violin)
p2 <- ggplot(contamination_frac, aes(x = X, y = contamination_frac)) +
  geom_violin() +
  geom_jitter(shape=16, position=position_jitter(0.2)) +
  theme_Publication_blank() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

# Export plots
# set_panel_size(p1, file="{soupx_plot_out_dir}/contaminationfrac_box.pdf" |> glue(),
#                width=unit(5, "in"), height=unit(4, "in"))

set_panel_size(p2, file="{soupx_plot_out_dir}/contaminationfrac_violin.pdf" |> glue(),
               width=unit(5, "in"), height=unit(4, "in"))

