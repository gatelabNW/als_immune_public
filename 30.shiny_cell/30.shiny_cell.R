# Charles Zhang
# Gate Lab
# Northwestern University
# 09/26/2023
####################################################################
# Generate ShinyCell App
####################################################################

# Load in libraries
suppressMessages({
  library("Seurat")
  library("plyr")
  library("tidyverse")
  library("ShinyCell")
  library("rsconnect")
  library("glue")
})

rsconnect::setAccountInfo(name='gatelabnu',
                          token='4BBF2E37652795198FBF498B31389548',
                          secret='3+KuwrOubk3mYRI79UyPzhseY7S2YUNmqeX/OxwQ')

# Load in data
# seurat_path <- "/projects/b1042/Gate_Lab/projects/als-project/scrna_original_final/04.seurat/seurat_SCT_mapped_04_03.rds"
# app_name <- "ALS_immune_profile"
# out_dir <- "/projects/b1042/Gate_Lab/projects/als-project/scrna_original_final/30.shiny_apps/immune_profile"

# seurat_path <- "/projects/b1042/Gate_Lab/projects/als-project/immune_enriched_final/04.seurat/seurat_SCT_mapped_04_02.rds"
# app_name <- "ALS_immune_enriched"
# out_dir <- "/projects/b1042/Gate_Lab/projects/als-project/immune_enriched_final/30.shiny_apps/immune_enriched"

# seurat_path <- "/projects/b1042/Gate_Lab/projects/als-project/crispr_clean_final/04.seurat/seurat_SCT_mapped_04_02.rds"
# app_name <- "ALS_crispr_clean"
# out_dir <- "/projects/b1042/Gate_Lab/projects/als-project/crispr_clean_final/30.shiny_apps/crispr_clean"

seurat_path <- "/projects/b1042/Gate_Lab/projects/als-project/crispr_clean_final/04.seurat/seurat_add_QC_04_04.rds"
app_name <- "ALS_crispr_clean"
out_dir <- "/projects/b1042/Gate_Lab/projects/als-project/crispr_clean_final/30.shiny_apps/crispr_clean"

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

s <- readRDS(seurat_path)

# modify celltype color
celltype_colors_path <- "/projects/b1169/nat/als/resources/metadata/cluster-metadata.csv"
celltype_colors <- read.csv(celltype_colors_path)
celltype_colors <- celltype_colors[which(celltype_colors$predicted.celltype.l2 %in% unique(s[["predicted.celltype.l2"]])[,1]),]
s$predicted.celltype.l2 <- factor(s$predicted.celltype.l2,
                                  levels = celltype_colors$predicted.celltype.l2)


# Create configuration
scConf <- createConfig(s)

# modify default display meta data
scConf <- modDefault(scConf, "predicted.celltype.l2", "age")

# Modify cell type colors
scConf <- modColours(scConf, meta.to.mod = "predicted.celltype.l2",
                     new.colours= celltype_colors$color)

# Show output
showLegend(scConf)
showOrder(scConf)

setwd(out_dir)
makeShinyApp(s, scConf, gene.mapping = FALSE, shiny.title = app_name)

# rename
shiny_folder <- glue("{out_dir}/shinyApp")
final_folder <- str_replace(shiny_folder, "shinyApp", app_name)
file.rename(shiny_folder, final_folder)

# deploy with notebook as it may ask for dependencies and stuff
# rsconnect::deployApp(final_folder)
# rsconnect::configureApp(app_name, size="xlarge")