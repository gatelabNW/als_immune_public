#!/bin/bash
#SBATCH --account b1042
#SBATCH --partition genomicsguestA
#SBATCH --job-name 06.2
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 8
#SBATCH --mem 64G
#SBATCH --time 8:00:00
#SBATCH --output=/projects/b1042/Gate_Lab/zzhang/cluster_logs/%j.%N.txt
#SBATCH --error=/projects/b1042/Gate_Lab/zzhang/cluster_logs/%j.%N.err
#SBATCH --verbose



date

# Module prep
module purge all

module load R/4.2.3

# Navigate to working directory
cd /projects/p31535/zzhang/als/als_repo/06.adaptive_clonal_analysis

Rscript 06.2.clonality_all_genes.R