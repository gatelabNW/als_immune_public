#!/bin/bash
#SBATCH --account p31535
#SBATCH --partition genhimem
#SBATCH --job-name sct_map_lbl
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 8
#SBATCH --mem 800G
#SBATCH --time 8:00:00
#SBATCH --output=/projects/b1042/Gate_Lab/zzhang/cluster_logs/%j.%N.txt
#SBATCH --error=/projects/b1042/Gate_Lab/zzhang/cluster_logs/%j.%N.err
#SBATCH --verbose


date

# Load modules
module purge
module load geos/3.8.1
module load R/4.1.1


####################################################################

# Navigate to working directory
cd /projects/p31535/zzhang/als/als_repo/04.integration_and_clustering

Rscript crispr_clean_sct_and_assign_labels.R
