#!/bin/bash
#SBATCH --account b1042
#SBATCH --partition genomicsguestA
#SBATCH --job-name 04
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 20
#SBATCH --mem 180G
#SBATCH --time 24:00:00
#SBATCH --output=/projects/b1042/Gate_Lab/zzhang/cluster_logs/04.%x_oe%j.txt
#SBATCH --error=/projects/b1042/Gate_Lab/zzhang/cluster_logs/04.%x_oe%j.err
#SBATCH --verbose



date

# Module prep
module purge all
module load R/4.4.0

# Navigate to working directory
cd /projects/p31535/zzhang/als/als_repo/05.DEG/crispr_clean_DEG

Rscript 04.differential_expression_sALS.R