#!/bin/bash
# Charles Zhang
# Gate Lab
# Northwestern University
# 10/10/2023
####################################################################
# Batch run the basecaller script for each sample

model="/projects/b1169/zzhang/ont_reference/dorado_models/dna_r10.4.1_e8.2_400bps_sup@v4.2.0"
out_root="/projects/b1042/Gate_Lab/projects/als-project/ont_seq"
raw_root="/projects/b1042/Gate_Lab/projects/als-project/ont_seq/00.raw"

while IFS= read -r dir; do
    # Your code to handle each directory goes here
    cur_sample=$(echo "$dir" | rev | awk -F'/' '{print $3}' | rev)
    cur_out_dir=$(echo "$dir" | sed 's/pod5_skip/fastq_pass/')
    echo "Current sample: $cur_sample"
    echo "Processing pod5 directory: $dir"
    echo "Output directory: $cur_out_dir"
    # submit jobs to quest
#    sbatch 20.2.batch_basecall.sh "${model}" "${dir}" "${cur_out_dir}"
done < <(find "$raw_root" -type d -name "*pod5_skip" -exec realpath {} \;)