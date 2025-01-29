#!/bin/bash
#SBATCH --account p31535
#SBATCH --partition short
#SBATCH --job-name add_qc
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 8
#SBATCH --mem 180G
#SBATCH --time 3:00:00
#SBATCH --output=/projects/b1042/Gate_Lab/zzhang/cluster_logs/%j.%N.txt
#SBATCH --error=/projects/b1042/Gate_Lab/zzhang/cluster_logs/%j.%N.err
#SBATCH --verbose



date

# Load modules
module purge all
module load R/4.2.3

####################################################################

# Navigate to working directory
cd /projects/p31535/zzhang/als/als_repo/04.integration_and_clustering

Rscript crispr_clean_remove_non_bt_cells_w_receptors.R