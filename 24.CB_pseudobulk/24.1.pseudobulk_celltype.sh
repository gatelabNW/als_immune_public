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

# barcodes output from 21.2
#barcode_file_dir="/projects/b1042/Gate_Lab/projects/als-project/crispr_clean_final/21.celltype_specific_reads/barcodes_l1"
barcode_file_dir="/projects/b1042/Gate_Lab/projects/als-project/crispr_clean_final/21.celltype_specific_reads/barcodes_l2"
#path_l2="/projects/b1042/Gate_Lab/projects/als-project/crispr_clean_final/21.celltype_specific_reads/barcodes_l2"

search_root="/projects/b1042/Gate_Lab/projects/als-project/ont_seq/03.CB_alignment_no_consensus_v2"
# Define the directory path
dir_path="${barcode_file_dir}"

# Initialize an empty array
files=()

# list all barcode files
while IFS= read -r -d $'\0' file; do
  files+=("$file")
done < <(find "$dir_path" -maxdepth 1 -type f -print0)

# list all sample dir
all_sample_dir=$(find "${search_root}" -maxdepth 1 -type d | grep -v "^$search_root$")
echo "${all_sample_dir}"
# Iterate through each sample
for cur_sample_dir in ${all_sample_dir}; do
  cur_sample_name="${cur_sample_dir##*/}"
  echo "Info: working on sample ${cur_sample_name}"

  # only rerunning on sample
  if [ "$cur_sample_name" != "F1" ]; then
    # The strings are equal, so replace original_string with replacement_string
    echo "Skipping ${cur_sample_name}"
    continue
  fi

  # replace F1 with E1 barcodes
  if [ "$cur_sample_name" = "F1" ]; then
    # The strings are equal, so replace original_string with replacement_string
    cur_sample_name="E1"
    echo "F1 to E1 replacement"
  fi

  # create log file
  log_file="${cur_sample_dir}/log.txt"
  > "${log_file}"

  # create output dir
  cur_level="${barcode_file_dir##*_}"
  cur_sample_bam_out="${cur_sample_dir}/${cur_level}"
  if [ ! -d "$cur_sample_bam_out" ]; then
    mkdir -p "$cur_sample_bam_out"
  else
    echo "Folder already exists: $cur_sample_bam_out"
  fi

  # Iterate through each celltype level file
  for file in "${files[@]}"; do
    echo "Info: working on ${file}"

    # initialize a bam array
    cur_sample_dir_bam_files=()

    # parse the cell type level name from the cell type specific barcodes file name
    cur_barcodes_ct_name="${file##*/}"
    cur_barcodes_ct_name="${cur_barcodes_ct_name%.*}"

    # keep track of the current cell type total count
    counter=0
    while IFS= read -r line
    do
      # extract barcode and find target bam file
      cur_barcode="${line//-1/}"
      # parse the sample underscore barcode format into two separate entities
      cur_barcode_no_sample_id="${cur_barcode#*_}"
      cur_sample_in_barcode="${cur_barcode%%_*}"

      # skip unmathcing sample to save time
      if [[ "$cur_sample_in_barcode" != "$cur_sample_name" ]]; then
        continue
      fi

      # update counter
      ((counter++))

      result=$(find "${cur_sample_dir}" -type f -name "${cur_barcode_no_sample_id}.curated.minimap2.bam")
      # Check if the result is not empty
      if [ -n "$result" ]; then
        cur_sample_dir_bam_files+=("$result")
      fi
    done < "${file}"
    num_cell=${#cur_sample_dir_bam_files[@]}
    echo "INFO: ${cur_barcodes_ct_name} has ${num_cell} matching cells out of ${counter} cells." >> "${log_file}"

    # samtools merge
    cur_sample_celltype_bam="${cur_sample_bam_out}/${cur_barcodes_ct_name}.bam"
    # Merge the BAM files
    bam_array_length=${#cur_sample_dir_bam_files[@]}
    if [ "$bam_array_length" -gt 2 ]; then
        samtools merge "${cur_sample_celltype_bam}" "${cur_sample_dir_bam_files[@]}"
    fi


    # Check if samtools merge was successful
    if [ $? -eq 0 ]; then
      echo "Merge successful, indexing the BAM file..."

      # Index the merged BAM file
      samtools index "${cur_sample_celltype_bam}"

      # Check if samtools index was successful
      if [ $? -eq 0 ]; then
        echo "Indexing successful"
      else
        echo "Indexing failed"
      fi
    else
      echo "Merge failed"
    fi
  done
done
