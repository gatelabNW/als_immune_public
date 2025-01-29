#!/bin/bash
#SBATCH --account p31535
#SBATCH --partition genhimem
#SBATCH --job-name DGGDE
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 20
#SBATCH --mem 500G
#SBATCH --time 24:00:00
#SBATCH --output=/projects/b1042/Gate_Lab/zzhang/cluster_logs/DGGDE.txt
#SBATCH --error=/projects/b1042/Gate_Lab/zzhang/cluster_logs/DGGDE.err
#SBATCH --verbose



date

# Module prep
module purge all
module load R/4.2.3

# Navigate to working directory
cd /projects/p31535/zzhang/als/als_repo/05.DEG/crispr_clean_DEG

Rscript differential_expression_diagnosis_general.R