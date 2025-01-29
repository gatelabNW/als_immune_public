#!/bin/bash
#SBATCH --account p31535
#SBATCH --partition short
#SBATCH --job-name CLU_VIS
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 8
#SBATCH --mem 80G
#SBATCH --time 1:00:00
#SBATCH --output=/projects/b1042/Gate_Lab/zzhang/cluster_logs/%j.%N.txt
#SBATCH --error=/projects/b1042/Gate_Lab/zzhang/cluster_logs/%j.%N.err
#SBATCH --verbose


date

# Load modules
module purge
module load geos/3.8.1
module load R/4.1.1


####################################################################

# Navigate to working directory
cd /projects/p31535/zzhang/als/als_repo/04.integration_and_clustering

Rscript cluster_visualizations.R
