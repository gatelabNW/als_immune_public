#!/bin/bash
#SBATCH --account b1042
#SBATCH --partition genomicsguestA
#SBATCH --job-name FASTQC
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 52
#SBATCH --mem 150G
#SBATCH --time 12:00:00
#SBATCH --output=/projects/b1042/Gate_Lab/zzhang/cluster_logs/%j.%N.txt
#SBATCH --error=/projects/b1042/Gate_Lab/zzhang/cluster_logs/%j.%N.err
#SBATCH --verbose


module load fastqc
module load multiqc

input_dir="/projects/b1042/Gate_Lab/projects/als-project/ont_seq/00.raw/E4/20230907_1711_P2S-00717-A_PAO96600_0d01e666/fastq_pass"
output_dir="/projects/b1042/Gate_Lab/projects/als-project/ont_seq/00.raw/E4/20230907_1711_P2S-00717-A_PAO96600_0d01e666"

find "${input_dir}" -name '*.fastq.gz' -exec fastqc -t 52 {} \;
