#!/bin/bash
#SBATCH --account p31535
#SBATCH --partition short
#SBATCH --job-name B_ASNER
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 16
#SBATCH --mem 180G
#SBATCH --time 4:00:00
#SBATCH --output=/projects/b1042/Gate_Lab/zzhang/cluster_logs/%j.%N.txt
#SBATCH --error=/projects/b1042/Gate_Lab/zzhang/cluster_logs/%j.%N.err
#SBATCH --verbose

# Charles Zhang
# Gate Lab
# Northwestern University
# 10/11/2023
####################################################################
# Batch script to submit to Quest for running scNanoGPS assigner
# Default minimum number of cells is 3000

source activate scNanoGPS

# Assign arguments to named variables for clarity
SC_NANO="$1"
out_root="$2"
sample_id="$3"

in_dir="${out_root}/01.preprocess"
input_cell_barcodes="${out_root}/01.preprocess/${sample_id}/barcode_list.tsv.gz"
out_dir="${out_root}/02.CB_collapse/${sample_id}"

# get 10X estimate of cell number from crispr clean data
CRISPR_clean_dir="/projects/b1042/Gate_Lab/projects/als-project/crispr_clean_final/01.cellranger_count"
metrics_10X="${CRISPR_clean_dir}/${sample_id}/outs/metrics_summary.csv"
cur_cell_num=$(awk -F'^"|","|"$' 'NR==2 {gsub(/^"|"$/,"",$2); gsub(/,/,"",$2); print $2; exit}' "${metrics_10X}")

echo "Output dir is ${out_dir}"
echo "Estimated cell number is ${cur_cell_num}"

if [ ! -d "$out_dir" ]; then
    mkdir "$out_dir"
fi
python3 "${SC_NANO}"/assigner.py -i "${input_cell_barcodes}" -d "${out_dir}" -t 16 --forced_no "${cur_cell_num}" --tmp_dir "${out_dir}/tmp"
