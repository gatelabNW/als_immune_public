#!/bin/bash
#SBATCH --account p31535
#SBATCH --partition genhimem
#SBATCH --job-name 04.2
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 8
#SBATCH --mem 800G
#SBATCH --time 8:00:00
#SBATCH --output=/projects/b1042/Gate_Lab/zzhang/cluster_logs/%j.%N.txt
#SBATCH --error=/projects/b1042/Gate_Lab/zzhang/cluster_logs/%j.%N.err
#SBATCH --verbose


date

# Load modules
module purge all
module load R/4.2.3
module load hdf5/1.10.6-openmpi-3.1.3-gcc-8.4.0-R


####################################################################

# Navigate to working directory
cd /projects/p31535/zzhang/als/als_repo/04.integration_and_clustering

Rscript 04.2.crispr_clean_sct_and_assign_labels.R
