#!/bin/bash
#SBATCH --account p31535
#SBATCH --partition short
#SBATCH --job-name shiny_cell
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 8
#SBATCH --mem 80G
#SBATCH --time 1:00:00
#SBATCH --output=/projects/b1042/Gate_Lab/zzhang/cluster_logs/%j.%N.txt
#SBATCH --error=/projects/b1042/Gate_Lab/zzhang/cluster_logs/%j.%N.err
#SBATCH --verbose

module load R/4.2.3
module load hdf5/1.8.19-serial

cd /projects/p31535/zzhang/als/als_repo/40.spatial_reference
Rscript 40.spatial_reference_shiny.R