#!/bin/bash
#SBATCH --account b1042
#SBATCH --partition genomicsguestA
#SBATCH --job-name 43.2
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 16
#SBATCH --mem 150G
#SBATCH --time 24:00:00
#SBATCH --output /projects/b1042/Gate_Lab/zzhang/cluster_logs/%x_oe%j.log
#SBATCH --verbose

date

# Module prep
module purge all
module load R/4.2.3

cluster_id=$1
echo "$cluster_id"

cd /projects/p31535/zzhang/als/als_repo/43.downstream_analysis
Rscript 43.2.re_DE_condition_general.R "$cluster_id"

