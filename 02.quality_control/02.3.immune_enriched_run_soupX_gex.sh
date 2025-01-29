#!/bin/bash
#SBATCH --account b1042
#SBATCH --partition genomicsguestA
#SBATCH --job-name imm_enr_soupX
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 8
#SBATCH --mem 64G
#SBATCH --time 8:00:00
#SBATCH --output=/projects/b1042/Gate_Lab/zzhang/cluster_logs/%j.%N.txt
#SBATCH --error=/projects/b1042/Gate_Lab/zzhang/cluster_logs/%j.%N.err
#SBATCH --verbose



date

# Load modules
module purge all
module load R/4.4.0

####################################################################

# Navigate to working directory
cd /projects/p31535/zzhang/als/als_repo/02.quality_control

# Run script
Rscript 02.3.immune_enriched_soupX_gex.R