#!/bin/bash
#SBATCH --account b1042
#SBATCH --partition genomicsguestA
#SBATCH --job-name 25.5
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 52
#SBATCH --mem 100G
#SBATCH --time 2:00:00
#SBATCH --output=/projects/b1042/Gate_Lab/zzhang/cluster_logs/%j.%N.txt
#SBATCH --error=/projects/b1042/Gate_Lab/zzhang/cluster_logs/%j.%N.err
#SBATCH --verbose

conda activate isoquant
cd /projects/b1042/Gate_Lab/projects/als-project/ont_seq/04.isoquant/04.isoquant_vis/CD16_Mono_batch
python /home/zzj4347/softwares/anaconda3/envs/isoquant/bin/visualize.py /projects/b1042/Gate_Lab/projects/als-project/ont_seq/04.isoquant/03.isoquant_out/CD16_Mono_batch --gene_list /projects/b1042/Gate_Lab/projects/als-project/ont_seq/04.isoquant/04.isoquant_vis/CD16_Mono_batch/genes_to_explore.txt