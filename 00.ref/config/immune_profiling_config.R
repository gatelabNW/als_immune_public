# ------------------------------------------------------------------------------
# -----                                                                    -----
# -----                           Gate Lab                                 -----
# -----                     Northwestern University                        -----
# -----                                                                    -----
# ------------------------------------------------------------------------------
#
# Date: 06-01-2023
# Written by: Ziyang Zhang
# Summary: Configuration file sourced by each analysis file for the ALS project
#
#-------------------------------------------------------------------------------

# ------------------------------ INPUT PARAMETERS ------------------------------
out_dir_root <- "/projects/b1042/Gate_Lab/projects/als-project/scrna_original_final"
repo_root <- "/projects/p31535/zzhang/als/als_repo"
meta <- "/projects/p31535/zzhang/als/als_repo/00.ref/meta/samples-metadata.csv"
color_path<-"/projects/p31535/zzhang/als/als_repo/00.ref/meta/cluster-metadata.csv"

dyn.load("/software/hdf5/1.8.19-serial/lib/libhdf5_hl.so.10")

# ---------------------------------- LIBRARIES ----------------------------------
suppressMessages({
  library(tidyverse)
  library(Seurat)
  library(DoubletFinder)
  library(glue)
  library(data.table)
  library(Matrix)
  library(SoupX)
  library(ggplot2)
  library(DropletUtils)
  library(purrr)
  library(Hmisc)
  library(parallel)
  library(dplyr)
  #library(SeuratDisk)
  library(MAST)
  library(gtools)
  library(argparse)
  library(stringr)
  library(ggrepel)
  library(ggthemes)
  library(grid)
  library(qpdf)
  library(UpSetR)
  library(stringr)
  library(tidyr)

})

options(echo = TRUE)
`%notin%` <- Negate(`%in%`)
source(glue("{repo_root}/00.lib/plot/ggplot_formatting.R"))
source(glue("{repo_root}/00.lib/plot/volcano_plot.R"))
source(glue("{repo_root}/00.lib/util/generate_degs.R"))
source(glue("{repo_root}/00.lib/util/run_doubletfinder.R"))

# --------------------------------- PARAMTERS ---------------------------------




