#!/bin/bash
#SBATCH --account b1042
#SBATCH --partition genomicsguestA
#SBATCH --job-name 02
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 4
#SBATCH --mem 120G
#SBATCH --time 12:00:00
#SBATCH --output=/projects/b1042/Gate_Lab/zzhang/cluster_logs/02.%x_oe%j.txt
#SBATCH --error=/projects/b1042/Gate_Lab/zzhang/cluster_logs/02.%x_oe%j.err
#SBATCH --verbose



date

# Module prep
module purge all
module load R/4.2.3

# Navigate to working directory
cd /projects/p31535/zzhang/als/als_repo/05.DEG/crispr_clean_DEG_adaptive_only

Rscript 02.differential_expression_c9_female_healthy_only.R