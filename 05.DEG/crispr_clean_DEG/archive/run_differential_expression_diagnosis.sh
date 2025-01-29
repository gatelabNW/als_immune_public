#!/bin/bash
#SBATCH --account p31535
#SBATCH --partition genhimem
#SBATCH --job-name DGDE
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 20
#SBATCH --mem 500G
#SBATCH --time 36:00:00
#SBATCH --output=/projects/b1042/Gate_Lab/zzhang/cluster_logs/%j.%N.txt
#SBATCH --error=/projects/b1042/Gate_Lab/zzhang/cluster_logs/%j.%N.err
#SBATCH --verbose



date

# Module prep
module purge all
module load R/4.2.3

# Navigate to working directory
cd /projects/p31535/zzhang/als/als_repo/05.DEG/crispr_clean_DEG

Rscript differential_expression_diagnosis.R