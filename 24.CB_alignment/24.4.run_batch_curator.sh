#!/bin/bash
# Charles Zhang
# Gate Lab
# Northwestern University
# 10/10/2023
####################################################################
# Batch run the scanner script for each sample in the raw_root directory

SC_NANO="/home/zzj4347/softwares/scNanoGPS/scNanoGPS"
output_root="/projects/b1042/Gate_Lab/projects/als-project/ont_seq"
input_CB_root="/projects/b1042/Gate_Lab/projects/als-project/ont_seq/02.CB_collapse"
input_FQ_root="/projects/b1042/Gate_Lab/projects/als-project/ont_seq/01.preprocess"
ref_g="/projects/b1169/zzhang/ont_reference/hg38.fa.gz"
ref_g_idx="/projects/b1169/zzhang/ont_reference/ref_hg38.mmi"

while IFS= read -r dir; do
    cur_sample=$(echo "$dir" | rev | awk -F'/' '{print $1}' | rev)
    echo "Current sample: $cur_sample"

    # Skip the sample if it is A4 or B2
    if [ "$cur_sample" = "A4" ] || [ "$cur_sample" = "B2" ]; then
        echo "Skipping sample $cur_sample"
        continue
    fi

    # Skip the sample if it is A4 or B2
    if [ "$cur_sample" = "E4" ]; then
        echo "Skipping sample $cur_sample"
        continue
    fi

    # Skip the sample if it is A4 or B2
    if [ "$cur_sample" = "A4" ] || [ "$cur_sample" = "B2" ]; then
        echo "Skipping sample $cur_sample"
        continue
    fi

    # Skip the sample if it is A4 or B2
    if [ "$cur_sample" = "F3" ]; then
        echo "Skipping sample $cur_sample"
        continue
    fi

    # Skip the sample if it is A4 or B2
    if [ "$cur_sample" = "H5" ]; then
        echo "Skipping sample $cur_sample"
        continue
    fi

    cur_FQ="${input_FQ_root}/${cur_sample}"
    cur_CB="${input_CB_root}/${cur_sample}"
    # submit jobs to quest
    sbatch 24.1.run_curator_highmem.sh "${SC_NANO}" "${cur_FQ}" "${cur_CB}" "${output_root}" "${cur_sample}" "${ref_g}" "${ref_g_idx}"
done < <(find "${input_FQ_root}" -mindepth 1 -type d -exec realpath {} \;)