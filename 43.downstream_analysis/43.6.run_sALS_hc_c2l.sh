#!/bin/bash
#SBATCH --account p31535
#SBATCH --partition short
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
module load R/4.4.0
module load hdf5/1.14.1-2-gcc-12.3.0
module load gsl/2.7.1-gcc-12.3.0
module load fftw/3.3.10-gcc-12.3.0
module load gdal/3.7.0-gcc-12.3.0
module load nlopt/2.7.1-gcc-12.3.0
cluster_id=$1
echo "$cluster_id"

cd /projects/p31535/zzhang/als/als_repo/43.downstream_analysis
Rscript 43.6.sALS_hc__c2l_enriched_RE_DE.R "$cluster_id"

