#!/bin/bash

# barcodes output from 21.2
path_l1="/projects/b1042/Gate_Lab/projects/als-project/crispr_clean_final/21.celltype_specific_reads/barcodes_l1"
path_l2="/projects/b1042/Gate_Lab/projects/als-project/crispr_clean_final/21.celltype_specific_reads/barcodes_l2"

# define samples identity
HC_samples=("B5" "B4" "H4" "E1" "G5" "C4" "A4" "E4" "F2" "A5" "A3" "E3" "D2" "D1" "A2" "B2" "F1" "H1")
ALS_samples=("H2" "G3" "E5" "B3" "H3" "A1" "G1" "C5" "E2" "F5" "D5" "D4" "F4" "G4" "C3" "G2" "B1" "C1" "F3" "C2" "H5" "D3")

#assays=("crispr_clean_final" "immune_enriched_final" "scrna_original_final")
assays=("crispr_clean_final")
for cur_assay in "${assays[@]}"; do

  # l1 first
  output_dir_root="/projects/b1042/Gate_Lab/projects/als-project/${cur_assay}/21.celltype_specific_reads/bams_l1"
  cur_ct_path="${path_l1}"
  if [ ! -d "$output_dir_root" ]; then
    mkdir -p "$output_dir_root"
  fi

  ALS_out_dir="${output_dir_root}/ALS"
  HC_out_dir="${output_dir_root}/HC"

  if [ ! -d "$ALS_out_dir" ]; then
    mkdir "$ALS_out_dir"
  fi
  if [ ! -d "$HC_out_dir" ]; then
    mkdir "$HC_out_dir"
  fi


  for cur_sample in "${ALS_samples[@]}"; do
    sbatch 21.3.batch_subset_celltype_specific_reads.sh "${cur_ct_path}" "${cur_sample}" "${cur_assay}" "${ALS_out_dir}"
  done

#   Work on HC controls
  for cur_sample in "${HC_samples[@]}"; do
    sbatch 21.3.batch_subset_celltype_specific_reads.sh "${cur_ct_path}" "${cur_sample}" "${cur_assay}" "${HC_out_dir}"
  done


  # l2 second
  output_dir_root="/projects/b1042/Gate_Lab/projects/als-project/${cur_assay}/21.celltype_specific_reads/bams_l2"
  cur_ct_path="${path_l2}"

  if [ ! -d "$output_dir_root" ]; then
      mkdir -p "$output_dir_root"
  fi

  ALS_out_dir="${output_dir_root}/ALS"
  HC_out_dir="${output_dir_root}/HC"

  if [ ! -d "$ALS_out_dir" ]; then
    mkdir "$ALS_out_dir"
  fi
  if [ ! -d "$HC_out_dir" ]; then
    mkdir "$HC_out_dir"
  fi


  for cur_sample in "${ALS_samples[@]}"; do
    sbatch 21.3.batch_subset_celltype_specific_reads.sh "${cur_ct_path}" "${cur_sample}" "${cur_assay}" "${ALS_out_dir}"
  done

  # Work on HC controls
  for cur_sample in "${HC_samples[@]}"; do
    sbatch 21.3.batch_subset_celltype_specific_reads.sh "${cur_ct_path}" "${cur_sample}" "${cur_assay}" "${HC_out_dir}"
  done
done
