#!/bin/bash
# Charles Zhang
# Gate Lab
# Northwestern University
# 10/10/2023
####################################################################
# Batch run the scanner script for each sample in the raw_root directory

SC_NANO="/home/zzj4347/softwares/scNanoGPS/scNanoGPS"
out_root="/projects/b1042/Gate_Lab/projects/als-project/ont_seq"
#raw_root="/projects/b1042/Gate_Lab/projects/als-project/ont_seq/00.raw"

raw_root="/projects/b1042/Gate_Lab/projects/als-project/ont_seq/00.raw"
#raw_root="/projects/b1042/Gate_Lab/projects/als-project/ont_seq/01.isolated_reads"

while IFS= read -r dir; do
    # Your code to handle each directory goes here
    cur_sample=$(echo "$dir" | awk -F'/' '{print $9}')
    # check if empty
    if [ -z "$cur_sample" ]; then
        # Skip to the next iteration of the loop
        continue
    fi
    echo "Current sample: $cur_sample"
    echo "Processing fastq directory: $dir"
    # submit jobs to quest
    sbatch 22.1.batch_scanner.sh ${SC_NANO} ${dir} ${out_root} ${cur_sample}
done < <(find "$raw_root" -type d -name "*fastq_pass*" -exec realpath {} \;)