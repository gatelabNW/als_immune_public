#!/bin/bash
#SBATCH --account p31535
#SBATCH --partition short
#SBATCH --job-name CT_DIFF
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 8
#SBATCH --mem 80G
#SBATCH --time 2:00:00
#SBATCH --output=/projects/b1042/Gate_Lab/zzhang/cluster_logs/%j.%N.txt
#SBATCH --error=/projects/b1042/Gate_Lab/zzhang/cluster_logs/%j.%N.err
#SBATCH --verbose

date
module purge all
module load R/4.2.3

cd /projects/p31535/zzhang/als/als_repo/05.differential_cell_type
Rscript 05.1.normalized_differential_celltype.R
