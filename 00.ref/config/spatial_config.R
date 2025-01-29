# ------------------------------------------------------------------------------
# -----                                                                    -----
# -----                           Gate Lab                                 -----
# -----                     Northwestern University                        -----
# -----                                                                    -----
# ------------------------------------------------------------------------------
#
# Date: 11-10-2023
# Written by: Ziyang Zhang
# Summary: Configuration file sourced by each analysis file for the ALS project
#
#-------------------------------------------------------------------------------

# ------------------------------ INPUT PARAMETERS ------------------------------
out_dir_root <- "/projects/b1042/Gate_Lab/projects/als-project/spatial"
repo_root <- "/projects/p31535/zzhang/als/als_repo"
meta <- "/projects/p31535/zzhang/als/als_repo/00.ref/meta/samples-metadata.csv"
spatial_qc_meta <- "/projects/p31535/zzhang/als/als_repo/00.ref/meta/final_spatial_meta.csv"
spatial_sample_meta <- "/projects/p31535/zzhang/als/als_repo/00.ref/meta/sequencing-spatial-samples-ALS-fall-2023.csv"
# color_path<-"/projects/p31535/zzhang/als/als_repo/00.ref/meta/cluster-metadata.csv"


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
  # library(DropletUtils)
  library(purrr)
  library(Hmisc)
  library(parallel)
  library(dplyr)
  # library(SeuratDisk)
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