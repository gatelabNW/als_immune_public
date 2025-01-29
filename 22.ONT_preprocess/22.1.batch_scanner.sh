#!/bin/bash
#SBATCH --account p31535
#SBATCH --partition long
#SBATCH --job-name B_SCANNER
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 16
#SBATCH --mem 32G
#SBATCH --time 120:00:00
#SBATCH --output=/projects/b1042/Gate_Lab/zzhang/cluster_logs/%j.%N.txt
#SBATCH --error=/projects/b1042/Gate_Lab/zzhang/cluster_logs/%j.%N.err
#SBATCH --verbose

# Charles Zhang
# Gate Lab
# Northwestern University
# 10/10/2023
####################################################################
# Batch script to submit to Quest for scanning each sample

# Check that exactly 5 arguments are provided
if [ "$#" -ne 4 ]; then
    echo "Usage: $0 <SC_NANO> <input_fastq_dir> <out_root> <sample_id>"
    exit 1
fi

# Assign arguments to named variables for clarity
SC_NANO="$1"
input_fastq_dir="$2"
out_root="$3"
sample_id="$4"

# Print out provided arguments

echo "SC_NANO = ${SC_NANO}"
echo "Input FASTQ Directory = ${input_fastq_dir}"
echo "Output Root = ${out_root}"
echo "Sample ID = ${sample_id}"

# Activate environment
source activate scNanoGPS

# create current sample output dir
out_dir="${out_root}/01.preprocess/${sample_id}"
if [ ! -d "$out_dir" ]; then
    mkdir -p "$out_dir"
fi
echo "Output Directory = ${out_dir}"

# Update 7/12/2024 with benny's command
python3 ${SC_NANO}/scanner.py -i ${input_fastq_dir} -d ${out_dir} --a5 AAGCAGTGGTATCAACGCAGAG --a3 CTACACGACGCTCTTCCGATCT --pT TTTCTTATATGGG --lUMI 10 -t 16 --scanning_region 200 --debug_mode True