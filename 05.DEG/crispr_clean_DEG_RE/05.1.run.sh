#!/bin/bash
#SBATCH --account b1042
#SBATCH --partition genomicsguestA
#SBATCH --job-name 05.1
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 16
#SBATCH --mem 180G
#SBATCH --time 18:00:00
#SBATCH --output /projects/b1042/Gate_Lab/zzhang/cluster_logs/%x_oe%j.log
#SBATCH --verbose

date

# Module prep
module purge all
module load R/4.2.3

cluster_id=$1
echo "$cluster_id"

cd /projects/p31535/zzhang/als/als_repo/05.DEG/crispr_clean_DEG_RE
Rscript 05.1.differential_expression_diagnosis_general.R "$cluster_id"

