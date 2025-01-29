#!/bin/bash
#SBATCH --account b1042
#SBATCH --partition genomicsburst
#SBATCH --job-name CURATOR
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 52
#SBATCH --mem 600G
#SBATCH --time 168:00:00
#SBATCH --output=/projects/b1042/Gate_Lab/zzhang/cluster_logs/%j.%N.txt
#SBATCH --error=/projects/b1042/Gate_Lab/zzhang/cluster_logs/%j.%N.err
#SBATCH --verbose


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
out_dir="${out_root}/03.CB_alignment/${cur_sample}"
temp_out_dir="${out_root}/03.CB_alignment/${cur_sample}/tmp"

echo "Output dir is ${out_dir}"

source activate scNanoGPS

if [ ! -d "$out_dir" ]; then
    mkdir "$out_dir"
fi

if [ ! -d "$temp_out_dir" ]; then
    mkdir "$temp_out_dir"
fi

python3 ${SC_NANO}/curator.py --fq_name ${FQ_file} -d ${out_dir} -t 52 -b ${BC_file} --CB_count ${CB_count_file} --CB_list ${CB_merged_file} --ref_genome ${ref_g} --idx_genome ${ref_g_idx} --tmp_dir=${temp_out_dir}
