#!/bin/bash
#SBATCH --account b1042
#SBATCH --partition genomicsguestA
#SBATCH --job-name FASTQC
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 52
#SBATCH --mem 64G
#SBATCH --time 12:00:00
#SBATCH --output=/projects/b1042/Gate_Lab/zzhang/cluster_logs/%j.%N.txt
#SBATCH --error=/projects/b1042/Gate_Lab/zzhang/cluster_logs/%j.%N.err
#SBATCH --verbose


module load fastqc
module load multiqc

input_dir="/projects/b1042/Gate_Lab/projects/als-project/ont_seq/01.isolated_reads"
output_dir="/projects/b1042/Gate_Lab/projects/als-project/ont_seq/01.isolated_reads"

find "${input_dir}" -name '*merged_GGGGCC_reads.fastq.gz' -exec fastqc -t 52 {} \;
multiqc "${input_dir}" -o "${output_dir}"

