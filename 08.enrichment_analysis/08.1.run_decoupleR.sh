#!/bin/bash
#SBATCH --account b1042
#SBATCH --partition genomicsguestA
#SBATCH --job-name DCP
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 20
#SBATCH --mem 180G
#SBATCH --time 48:00:00
#SBATCH --output=/projects/b1042/Gate_Lab/zzhang/cluster_logs/%j.%N.txt
#SBATCH --error=/projects/b1042/Gate_Lab/zzhang/cluster_logs/%j.%N.err
#SBATCH --verbose

date

# Module prep
module purge
module load geos/3.8.1
module load R/4.2.3

cd /projects/p31535/zzhang/als/als_repo/08.enrichment_analysis
Rscript 08.1.decoupleR.R