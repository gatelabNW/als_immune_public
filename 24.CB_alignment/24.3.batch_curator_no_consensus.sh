#!/bin/bash
#SBATCH --account b1042
#SBATCH --partition genomicsguestA
#SBATCH --job-name CURATOR
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 52
#SBATCH --mem 180G
#SBATCH --time 48:00:00
#SBATCH --output=/projects/b1042/Gate_Lab/zzhang/cluster_logs/%j.%N.txt
#SBATCH --error=/projects/b1042/Gate_Lab/zzhang/cluster_logs/%j.%N.err
#SBATCH --verbose

# Charles Zhang
# Gate Lab
# Northwestern University
# 10/16/2023
####################################################################
# Batch script to submit to Quest for aligning each sample

source activate scNanoGPS

if [ "$#" -ne 7 ]; then
    echo "Usage: $0 <SC_NANO> <in_FQ_dir> <in_CB_dir> <out_root> <sample_id> <ref_g> <ref_g_idx>"
    exit 1
fi

SC_NANO="$1"
in_FQ_dir="$2"
in_CB_dir="$3"
out_root="$4"
cur_sample="$5"
ref_g="$6"
ref_g_idx="$7"

FQ_file="${in_FQ_dir}/processed.fastq.gz"
BC_file="${in_FQ_dir}/barcode_list.tsv.gz"
CB_count_file="${in_CB_dir}/CB_counting.tsv.gz"
CB_merged_file="${in_CB_dir}/CB_merged_list.tsv.gz"
out_dir="${out_root}/03.CB_alignment_no_consensus_v2/${cur_sample}"
temp_out_dir="${out_root}/03.CB_alignment_no_consensus_v2/${cur_sample}/tmp"

echo "Output dir is ${out_dir}"

if [ ! -d "${out_dir}" ]; then
    mkdir -p "${out_dir}"
    echo "Created ${out_dir}"

fi

if [ ! -d "${temp_out_dir}" ]; then
    mkdir -p "${temp_out_dir}"
    echo "Created ${temp_out_dir}"
fi

python3 ${SC_NANO}/curator.py --fq_name ${FQ_file} -d ${out_dir} -t 52 -b ${BC_file} --CB_count ${CB_count_file} --CB_list ${CB_merged_file} --ref_genome ${ref_g} --idx_genome ${ref_g_idx} --tmp_dir ${temp_out_dir} --skip_curation 1