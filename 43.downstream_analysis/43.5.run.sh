#!/bin/bash
#SBATCH --account b1042
#SBATCH --partition genomicsguestA
#SBATCH --job-name 43.5
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 16
#SBATCH --mem 80G
#SBATCH --time 12:00:00
#SBATCH --output /projects/b1042/Gate_Lab/zzhang/cluster_logs/%x_oe%j.log
#SBATCH --verbose

date

# Module prep
module purge all
module load R/4.4.0
module load nlopt/2.7.1-gcc-12.3.0


cluster_id=$1
echo "$cluster_id"

cd /projects/p31535/zzhang/als/als_repo/43.downstream_analysis
Rscript 43.5.re_DE_sALS_hc_manual_annotation.R "$cluster_id"
#Rscript 43.5.re_DE_c9_hc_manual_annotation.R "$cluster_id"
