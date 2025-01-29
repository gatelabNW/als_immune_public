#!/bin/bash
#SBATCH --account b1042
#SBATCH --partition genomicsguestA
#SBATCH --job-name 43.6
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 16
#SBATCH --mem 80G
#SBATCH --time 4:00:00
#SBATCH --output /projects/b1042/Gate_Lab/zzhang/cluster_logs/%x_oe%j.log
#SBATCH --verbose

date

# Module prep
module purge all
module load R/4.2.3

cluster_id=$1
echo "$cluster_id"

cd /projects/p31535/zzhang/als/als_repo/43.downstream_analysis
Rscript 43.6.c2l_enriched_RE_DE.R "$cluster_id"

