#!/bin/bash
#SBATCH --account b1042
#SBATCH --partition genomicsguestA
#SBATCH --job-name fmt_ctg_cltp
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 8
#SBATCH --mem 16G
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
cd /projects/p31535/zzhang/als/als_repo/03.process_vdj_output

Rscript format_contigs_and_clonotypes.R