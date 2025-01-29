#!/bin/bash
#SBATCH --account b1042
#SBATCH --partition genomicsburst
#SBATCH --job-name 25.4_array
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 52
#SBATCH --mem 1000G
#SBATCH --time 72:00:00
#SBATCH --output=/projects/b1042/Gate_Lab/zzhang/cluster_logs/%j.%N.txt
#SBATCH --error=/projects/b1042/Gate_Lab/zzhang/cluster_logs/%j.%N.err
#SBATCH --verbose

source activate isoquant
conda list IsoQuant

ref_gtf="/projects/b1169/zzhang/ont_reference/gencode.v46.chr_patch_hapl_scaff.annotation.gtf.gz"
ref_g="/projects/b1169/zzhang/ont_reference/hg38.fa.gz"

#yaml="/projects/b1042/Gate_Lab/projects/als-project/ont_seq/04.isoquant/02.ymal/A4_CD8_TEM_batch.yaml"
#output="/projects/b1042/Gate_Lab/projects/als-project/ont_seq/04.isoquant/03.isoquant_out/A4_CD8_TEM_batch"

# Base directory for inputs
input_dir="/projects/b1042/Gate_Lab/projects/als-project/ont_seq/04.isoquant/02.ymal"
# Base directory for outputs
output_base_dir="/projects/b1042/Gate_Lab/projects/als-project/ont_seq/04.isoquant/03.isoquant_out_3.6.1"
input_files=($input_dir/*)            # Create an array of all files in the directory


# Get the file corresponding to the current task ID
input_file=${input_files[$SLURM_ARRAY_TASK_ID-1]}  # Access the file for this task ID

# Extract the base file name (without directory path and suffix)
file_name=$(basename "$input_file")                 # Get the file name with extension
base_name="${file_name%.*}"                         # Remove the suffix (extension)

# Create an output directory for this file (based on the base file name)
output_dir="${output_base_dir}/${base_name}"        # Output directory based on file name
mkdir -p "$output_dir"                              # Create the directory if it doesn't exist



# Save information to the log file
log_file="$output_dir/process_log.txt"
echo "Job started at: $(date)" > "$log_file"
echo "Processing input file: $input_file" >> "$log_file"
echo "Output directory: $output_dir" >> "$log_file"

isoquant.py -d nanopore --yaml ${input_file} --complete_genedb --genedb ${ref_gtf} --reference ${ref_g} --output ${output_dir} -t 52