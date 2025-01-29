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

# Charles Zhang
# Gate Lab
# Northwestern University
# 10/10/2023
####################################################################
# Batch script to submit to Quest for basecalling skiped pod5 files

# Check that exactly 3 arguments are provided
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <model> <pod5_dir> <out_dir>"
    exit 1
fi

# Assign arguments to named variables for clarity
model="$1"
pod5_dir="$2"
out_dir="$3"

# Print out provided arguments

echo "model = ${model}"
echo "Input pod5 Directory = ${pod5_dir}"
echo "Output dir = ${out_dir}"

# Activate environment
source activate scNanoGPS

if [ ! -d "$out_dir" ]; then
    mkdir "$out_dir"
fi

dorado basecaller "${model}" "${pod5_dir}" --emit-fastq > "${out_dir}/calls.fastq"
gzip "${out_dir}/calls.fastq"
