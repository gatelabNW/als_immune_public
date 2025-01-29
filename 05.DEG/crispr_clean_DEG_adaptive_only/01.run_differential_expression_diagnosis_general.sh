#!/bin/bash
#SBATCH --account b1042
#SBATCH --partition genomicsguestA
#SBATCH --job-name 01
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 20
#SBATCH --mem 180G
#SBATCH --time 24:00:00
#SBATCH --output=/projects/b1042/Gate_Lab/zzhang/cluster_logs/01.%x_oe%j.log.txt
#SBATCH --error=/projects/b1042/Gate_Lab/zzhang/cluster_logs/01.%x_oe%j.log.err
#SBATCH --verbose



date

# Module prep
module purge all
module load R/4.2.3

# Navigate to working directory
cd /projects/p31535/zzhang/als/als_repo/05.DEG/crispr_clean_DEG_adaptive_only

Rscript 01.differential_expression_diagnosis_general.R