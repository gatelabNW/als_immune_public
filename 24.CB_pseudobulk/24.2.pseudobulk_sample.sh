#!/bin/bash
#SBATCH --account b1042
#SBATCH --partition genomicsguestA
#SBATCH --job-name BM_MG
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 8
#SBATCH --mem 32G
#SBATCH --time 12:00:00
#SBATCH --output=/projects/b1042/Gate_Lab/zzhang/cluster_logs/%j.%N.txt
#SBATCH --error=/projects/b1042/Gate_Lab/zzhang/cluster_logs/%j.%N.err
#SBATCH --verbose


search_root="/projects/b1042/Gate_Lab/projects/als-project/ont_seq/03.CB_alignment_no_consensus/03.CB_alignment_no_consensus"

# list all sample dir
all_sample_dir=$(find "${search_root}" -maxdepth 1 -type d | grep -v "^$search_root$")
echo "${all_sample_dir}"
# Iterate through each sample
for cur_sample_dir in ${all_sample_dir}; do
  cur_sample_name="${cur_sample_dir##*/}"
  echo "Info: working on sample ${cur_sample_name}"
  temp_dir="${cur_sample_dir}/temp_merge_dir"
  echo "temp folder ${temp_dir}"
  mkdir -p "$temp_dir"
  # create output dir
  cur_sample_bam_out="${cur_sample_dir}/all_sample"
  if [ ! -d "$cur_sample_bam_out" ]; then
    mkdir -p "$cur_sample_bam_out"
  else
    echo "Folder already exists: $cur_sample_bam_out"
  fi


  # initialize a bam array
  cur_sample_dir_bam_files=($(find "${cur_sample_dir}" -type f -name '*curated.minimap2.bam'))
  num_cell=${#cur_sample_dir_bam_files[@]}
  echo "INFO: ${cur_sample_name} has ${num_cell} matching cells."

  # samtools merge
  cur_sample_bam="${cur_sample_bam_out}/${cur_sample_name}_all.bam"

  batch_size=100
  batch_number=1
  batch_files=()

  for bam_file in "${cur_sample_dir_bam_files[@]}"; do
    batch_files+=("$bam_file")

    # When batch is full, merge and clear batch
    if [ ${#batch_files[@]} -eq $batch_size ]; then
      samtools merge "$temp_dir/batch_$batch_number.bam" "${batch_files[@]}"
      samtools index "$temp_dir/batch_$batch_number.bam"
      batch_files=()
      ((batch_number++))
    fi
  done
  # Merge any remaining files in the last batch
  if [ ${#batch_files[@]} -gt 0 ]; then
    samtools merge "$temp_dir/batch_$batch_number.bam" "${batch_files[@]}"
    samtools index "$temp_dir/batch_$batch_number.bam"
  fi

  # Merge all the batches
  samtools merge "$cur_sample_bam" "$temp_dir/*.bam"
  samtools index "$cur_sample_bam"

  # Clean up temporary directory
  rm -r "$temp_dir"
  echo "Merging complete. Merged file: $cur_sample_bam"

done
