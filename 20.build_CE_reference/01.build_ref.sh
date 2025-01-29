#!/bin/bash
#SBATCH --account b1042
#SBATCH --partition genomicsguestA
#SBATCH --job-name hg38_ref
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 8
#SBATCH --mem 150G
#SBATCH --time 24:00:00
#SBATCH --output=/projects/b1042/Gate_Lab/zzhang/cluster_logs/%j.%N.txt
#SBATCH --error=/projects/b1042/Gate_Lab/zzhang/cluster_logs/%j.%N.err
#SBATCH --verbose

date
source activate scNanoGPS

cd /projects/b1169/zzhang/ont_reference
minimap2 -d ref_hg38.mmi hg38.fa.gz

