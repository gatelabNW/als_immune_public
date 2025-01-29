# Charles Zhang
# Gate Lab
# Northwestern University
# 10/04/2023
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


seurat_path <- "/projects/b1169/zzhang/spatial_reference/als/snRNA-integrated-wMNs.rds"
app_name <- "ALS_Reference"
out_dir <- "/projects/b1169/zzhang/spatial_reference/als/shiny_app"

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

s <- readRDS(seurat_path)


# Create configuration
scConf <- createConfig(s)

# modify default display meta data
scConf <- modDefault(scConf, "AnnotationForDeconvolution", "top_level_annotation")


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