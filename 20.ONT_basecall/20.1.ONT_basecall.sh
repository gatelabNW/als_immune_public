#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics-gpu
#SBATCH --gres=gpu:a100:1
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 48:00:00
#SBATCH --mem=32G
#SBATCH --output=/projects/b1042/Gate_Lab/zzhang/cluster_logs/%j.%N.txt
#SBATCH --error=/projects/b1042/Gate_Lab/zzhang/cluster_logs/%j.%N.err
#SBATCH --verbose

model="/projects/b1169/zzhang/ont_reference/dorado_models/dna_r10.4.1_e8.2_400bps_sup@v4.2.0"
input_dir="/projects/b1042/Gate_Lab/projects/als-project/ont_seq/00.raw/F3/20230926_1425_P2S-00717-B_PAS47105_441fb3c4"
pod5_dir="${input_dir}/pod5_skip/"
out_dir="${input_dir}/pod5_skip_fastq"

if [ ! -d "$out_dir" ]; then
    mkdir "$out_dir"
fi

dorado basecaller "${model}" "${pod5_dir}" --emit-fastq > "${out_dir}/calls.fastq"