#!/bin/bash
#SBATCH --account p31535
#SBATCH --partition short
#SBATCH --job-name COR_CT
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 8
#SBATCH --mem 80G
#SBATCH --time 1:00:00
#SBATCH --output=/projects/b1042/Gate_Lab/zzhang/cluster_logs/%j.%N.txt
#SBATCH --error=/projects/b1042/Gate_Lab/zzhang/cluster_logs/%j.%N.err
#SBATCH --verbose

module load R/4.2.3

cd /projects/p31535/zzhang/als/als_repo/06.sc_clonal_analysis
Rscript 06.1.correct_post_QC_receptor_count.R