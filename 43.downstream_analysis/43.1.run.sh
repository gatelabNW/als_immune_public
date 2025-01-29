#!/bin/bash
#SBATCH --account b1042
#SBATCH --partition genomics
#SBATCH --job-name 43.1
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 16
#SBATCH --mem 80G
#SBATCH --time 24:00:00
#SBATCH --output /projects/b1042/Gate_Lab/zzhang/cluster_logs/43.1.%x_oe%j.log
#SBATCH --verbose

date

# Module prep
module purge all
module load R/4.2.3

cd /projects/p31535/zzhang/als/als_repo/43.downstream_analysis
Rscript 43.1.spatial_DE.R