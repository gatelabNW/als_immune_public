#!/bin/bash
#SBATCH --account p31535
#SBATCH --partition normal
#SBATCH --job-name CT_SUB
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 8
#SBATCH --mem 80G
#SBATCH --time 48:00:00
#SBATCH --output=/projects/b1042/Gate_Lab/zzhang/cluster_logs/%j.%N.txt
#SBATCH --error=/projects/b1042/Gate_Lab/zzhang/cluster_logs/%j.%N.err
#SBATCH --verbose

# Define Constant Vairables
RAW_DIR="/projects/b1042/Gate_Lab/projects/als-project"

# check for the four arguments passed in
if [ "$#" -ne 4 ]; then
    echo "Usage: $0 <ct_path> <sample_id> <assay_id> <out_dir>"
    exit 1
fi

# Assign arguments to named variables for clarity
ct_path="$1"
sample_id="$2"
assay_id="$3"
out_dir="$4"

# specify input dir
input_root_dir="${RAW_DIR}/${assay_id}/01.cellranger_count"

# l1 first
output_dir="${out_dir}/${sample_id}"
if [ ! -d "$output_dir" ]; then
  mkdir -p "$output_dir"
fi

# loop through the cell types barcodes to create a dictionary
declare -A barcode_to_celltype
for barcode_file in "${ct_path}"/*; do
  cell_type=$(basename "${barcode_file}" .txt) # Assuming file extensions are .txt. Adjust if needed.
  while read -r barcode; do
    barcode_to_celltype["${barcode}"]="${cell_type}"
  done <"$barcode_file"
done

echo "Begin ${sample_id}!"
input_bam="${input_root_dir}/${sample_id}/outs/possorted_genome_bam.bam"

# write header first for all cell type
echo "Writing headers to sam files!"
for barcode_file in "${ct_path}"/*; do
  cur_cell_type=$(basename "${barcode_file}" .txt) # Assuming file extensions are .txt. Adjust if needed.
  samtools view -H "${input_bam}" >"${output_dir}/${cur_cell_type}_${sample_id}.sam"
done

# Process the BAM file and split based on barcodes
echo "Processing input bam for ${sample_id}!"
samtools view "${input_bam}" | while read -r line; do
  if [[ ${line} =~ CB:Z:([A-Za-z0-9-]+) ]]; then
    # keep the regex match and store it
    barcode="${BASH_REMATCH[1]}"
    if [[ ${barcode_to_celltype["${barcode}"]+exists} ]]; then
      cur_cell_type="${barcode_to_celltype[${barcode}]}"
      echo "${line}" >>"${output_dir}/${cur_cell_type}_${sample_id}.sam"
    fi
  fi
done
echo "Processing input bam for ${sample_id} finished!"




