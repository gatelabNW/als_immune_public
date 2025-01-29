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

# barcodes output from 21.2
path_l1="/projects/b1042/Gate_Lab/projects/als-project/crispr_clean_final/21.celltype_specific_reads/barcodes_l1"
path_l2="/projects/b1042/Gate_Lab/projects/als-project/crispr_clean_final/21.celltype_specific_reads/barcodes_l2"

# define samples identity
HC_samples=("B5" "B4" "H4" "E1" "G5" "C4" "A4" "E4" "F2" "A5" "A3" "E3" "D2" "D1" "A2" "B2" "F1" "H1")
ALS_samples=("H2" "G3" "E5" "B3" "H3" "A1" "G1" "C5" "E2" "F5" "D5" "D4" "F4" "G4" "C3" "G2" "B1" "C1" "F3" "C2" "H5" "D3")

assays=("crispr_clean_final" "immune_enriched_final" "scrna_original_final")
for cur_assay in "${assays[@]}"; do
  # specify input dir
  input_root_dir="/projects/b1042/Gate_Lab/projects/als-project/${cur_assay}/01.cellranger_count"

  # l1 first
  output_dir_root="/projects/b1042/Gate_Lab/projects/als-project/${cur_assay}/21.celltype_specific_reads/bams_l1"
  if [ ! -d "$output_dir" ]; then
    mkdir -p "$output_dir"
  fi

  ALS_out_dir="${output_dir_root}/ALS"
  HC_out_dir="${output_dir_root}/HC"

  if [ ! -d "$ALS_out_dir" ]; then
    mkdir "$ALS_out_dir"
  fi
  if [ ! -d "$HC_out_dir" ]; then
    mkdir "$HC_out_dir"
  fi

  # loop through the cell types barcodes to create a dictionary
  declare -A barcode_to_celltype
  for barcode_file in "${path_l1}"/*; do
    cell_type=$(basename "${barcode_file}" .txt) # Assuming file extensions are .txt. Adjust if needed.
    while read -r barcode; do
      barcode_to_celltype["${barcode}"]="${cell_type}"
    done <"$barcode_file"
  done
  # Iterate over barcode files
  for cur_sample in "${ALS_samples[@]}"; do
    echo "Begin ${cur_sample}!"
    output_dir="${output_dir_root}/ALS"
    input_bam="${input_root_dir}/${cur_sample}/outs/possorted_genome_bam.bam"

    # write header first for all cell type
    echo "Writing headers to sam files!"
    for barcode_file in "${path_l1}"/*; do
      cur_cell_type=$(basename "${barcode_file}" .txt) # Assuming file extensions are .txt. Adjust if needed.
      samtools view -H "${input_bam}" >"${output_dir}/${cur_cell_type}_${cur_sample}.sam"
    done

    # Process the BAM file and split based on barcodes
    echo "Processing input bam for ${cur_sample}!"
    samtools view "${input_bam}" | while read -r line; do
      if [[ ${line} =~ CB:Z:([A-Za-z0-9-]+) ]]; then
        # keep the regex match and store it
        barcode="${BASH_REMATCH[1]}"
        if [[ ${barcode_to_celltype["${barcode}"]+exists} ]]; then
          cur_cell_type="${barcode_to_celltype[${barcode}]}"
          echo "${line}" >>"${output_dir}/${cur_cell_type}_${cur_sample}.sam"
        fi
      fi
    done
    echo "Processing input bam for ${cur_sample} finished!"
  done

  # convert to bams first
  echo "Converting all sams to bams!"
  for sam_file in "${output_dir}"/*.sam; do
    bam_file="${output_dir}/$(basename "$sam_file" .sam).bam"
    samtools view -Sb "$sam_file" >"$bam_file"
  done

  # merge
  echo "Merging and then sorting!"
  for barcode_file in "${path_l1}"/*; do
    cur_cell_type=$(basename "${barcode_file}" .txt) # Assuming file extensions are .txt. Adjust if needed.
    samtools merge -f "ALS_${cur_cell_type}.bam" "*${cur_cell_type}*.bam"
    samtools sort "ALS_${cur_cell_type}.bam" -o "ALS_${cur_cell_type}_sorted.bam"
    samtools index "ALS_${cur_cell_type}_sorted.bam"
  done

  # Work on HC controls
  for cur_sample in "${HC_samples[@]}"; do
    echo "${cur_sample}"
    output_dir="${output_dir_root}/HC"
    input_bam="${input_root_dir}/${cur_sample}/outs/possorted_genome_bam.bam"

    # write header first for all cell type
    for barcode_file in "${path_l1}"/*; do
      cur_cell_type=$(basename "${barcode_file}" .txt) # Assuming file extensions are .txt. Adjust if needed.
      samtools view -H "${input_bam}" >"${output_dir}/${cur_cell_type}_${cur_sample}.sam"
    done

    # Process the BAM file and split based on barcodes
    samtools view "${input_bam}" | while read -r line; do
      if [[ ${line} =~ CB:Z:([A-Za-z0-9-]+) ]]; then
        # keep the regex match and store it
        barcode="${BASH_REMATCH[1]}"
        if [[ ${barcode_to_celltype["${barcode}"]+exists} ]]; then
          cur_cell_type="${barcode_to_celltype[${barcode}]}"
          echo "${line}" >>"${output_dir}/${cur_cell_type}_${cur_sample}.sam"
        fi
      fi
    done
  done

  # convert to bams first
  for sam_file in "${output_dir}"/*.sam; do
    bam_file="${output_dir}/$(basename "$sam_file" .sam).bam"
    samtools view -Sb "$sam_file" >"$bam_file"
  done

  # merge
  echo "Merging and then sorting!"
  for barcode_file in "${path_l1}"/*; do
    cur_cell_type=$(basename "${barcode_file}" .txt) # Assuming file extensions are .txt. Adjust if needed.
    samtools merge -f "HC_${cur_cell_type}.bam" "*${cur_cell_type}*.bam"
    samtools sort "HC_${cur_cell_type}.bam" -o "HC_${cur_cell_type}_sorted.bam"
    samtools index "HC_${cur_cell_type}_sorted.bam"
  done

done
