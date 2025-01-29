#!/bin/bash
# Charles Zhang
# Gate Lab
# Northwestern University
# 10/11/2023
####################################################################
# Batch run the scanner script for each sample in the raw_root directory

SC_NANO="/home/zzj4347/softwares/scNanoGPS/scNanoGPS"
out_root="/projects/b1042/Gate_Lab/projects/als-project/ont_seq"
ONT_preprocess_dir="${out_root}/01.preprocess"

for dir in $(find "${ONT_preprocess_dir}" -maxdepth 1 -type d); do
    # exclude the root directory if you don't want to include it
    if [ "${dir}" != "${ONT_preprocess_dir}" ]; then
        # Your operations on the directory
        echo "Current dir is ${dir}"
        cur_sample=$(basename "${dir}")
        sbatch 23.2.batch_assigner.sh "${SC_NANO}" "${out_root}" "${cur_sample}"
    fi
done