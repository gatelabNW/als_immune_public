#!/bin/bash
#SBATCH --account p31535
#SBATCH --partition genhimem
#SBATCH --job-name general_QC
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 8
#SBATCH --mem 600G
#SBATCH --time 18:00:00
#SBATCH --output=/projects/b1042/Gate_Lab/zzhang/cluster_logs/%j.%N.txt
#SBATCH --error=/projects/b1042/Gate_Lab/zzhang/cluster_logs/%j.%N.err
#SBATCH --verbose



date

# Load modules
module purge
module load R/4.2.3

####################################################################

# Navigate to working directory
cd /projects/p31535/zzhang/als/als_repo/02.quality_control

Rscript 02.6.immune_profile_general_qc.R