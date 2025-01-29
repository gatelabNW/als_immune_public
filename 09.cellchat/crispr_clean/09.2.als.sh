#!/bin/bash
#SBATCH --account b1042
#SBATCH --partition genomicsguestA
#SBATCH --job-name 09.2
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 20
#SBATCH --mem 180G
#SBATCH --time 24:00:00
#SBATCH --output=/projects/b1042/Gate_Lab/zzhang/cluster_logs/%j.%N_09.2.txt
#SBATCH --error=/projects/b1042/Gate_Lab/zzhang/cluster_logs/%j.%N_09.2.err
#SBATCH --verbose



date

# Module prep
module load R/4.4.0
module load glpk/4.65-gcc-12.3.0

# Navigate to working directory
cd /projects/p31535/zzhang/als/als_repo/09.cellchat

Rscript 09.2.als.R