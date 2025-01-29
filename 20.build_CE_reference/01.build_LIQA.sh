#!/bin/bash
#SBATCH --account b1042
#SBATCH --partition genomicsguestA
#SBATCH --job-name LIQA_ref
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
liqa -task refgene -ref gencode.v38.chr_patch_hapl_scaff.annotation.gtf -format gtf -out gencode.v38.chr_patch_hapl_scaff.annotation.refgene

