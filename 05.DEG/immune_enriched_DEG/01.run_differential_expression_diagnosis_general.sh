#!/bin/bash
#SBATCH --account b1042
#SBATCH --partition genomicsguestA
#SBATCH --job-name differential_expression
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 20
#SBATCH --mem 150G
#SBATCH --time 12:00:00
#SBATCH --output=/projects/b1042/Gate_Lab/zzhang/cluster_logs/%j.%N.txt
#SBATCH --error=/projects/b1042/Gate_Lab/zzhang/cluster_logs/%j.%N.err
#SBATCH --verbose



date

# Module prep
module purge all
module load R/4.4.0

# Navigate to working directory
cd /projects/p31535/zzhang/als/als_repo/05.DEG/immune_enriched_DEG

Rscript 01.differential_expression_diagnosis_general.R