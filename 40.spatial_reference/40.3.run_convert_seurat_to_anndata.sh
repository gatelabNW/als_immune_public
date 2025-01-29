#!/bin/bash
#SBATCH --account b1042
#SBATCH --partition genomicsguestA
#SBATCH --job-name cv_seu
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 8
#SBATCH --mem 80G
#SBATCH --time 2:00:00
#SBATCH --output=/projects/b1042/Gate_Lab/zzhang/cluster_logs/%j.%N.txt
#SBATCH --error=/projects/b1042/Gate_Lab/zzhang/cluster_logs/%j.%N.err
#SBATCH --verbose

date

# Load modules
module purge all
module load R/4.4.0
module load hdf5/1.10.6-openmpi-3.1.3-gcc-8.4.0-R


cd /projects/p31535/zzhang/als/als_repo/40.spatial_reference
Rscript 40.3.convert_seurat_to_anndata.R