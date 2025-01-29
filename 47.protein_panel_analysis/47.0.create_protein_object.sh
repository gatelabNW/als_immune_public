#!/bin/bash
#SBATCH --account p31535
#SBATCH --partition short
#SBATCH --job-name 47.0
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 16
#SBATCH --mem 180G
#SBATCH --time 4:00:00
#SBATCH --output /projects/b1042/Gate_Lab/zzhang/cluster_logs/%x_oe%j.log
#SBATCH --verbose

module purge all
module load R/4.2.3
module load hdf5/1.10.6-openmpi-3.1.3-gcc-8.4.0-R


cd /projects/p31535/zzhang/als/als_repo/47.protein_panel_analysis
Rscript 47.0.create_protein_object.R