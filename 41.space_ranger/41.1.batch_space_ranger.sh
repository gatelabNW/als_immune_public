#!/bin/bash
#SBATCH --account b1042
#SBATCH --partition genomicsguestA
#SBATCH --job-name spaceranger
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 16
#SBATCH --mem 128G
#SBATCH --time 4:00:00
#SBATCH --output /projects/b1042/Gate_Lab/zzhang/cluster_logs/%x_oe%j.log
#SBATCH --verbose

if [[ "$#" -lt 10 ]]; then
    echo "Error: Too few arguments."
fi

# Define directories
#SAMPLE_ID="B_NMA22-205-B1"
#main_dir="/projects/b1042/Gate_Lab/lynn/AN1792/data/B1"
#ref_dir="/projects/p31535/spatial-transcriptomics/spacerange-ref-files-before-transfer"
#cytassist_image_dir="${main_dir}"
#hi_res_image_dir="${main_dir}"
#alignment_dir="${main_dir}"
#fastq_dir="${main_dir}/fastq/B1/"
#output_dir="/projects/b1042/Gate_Lab/lynn/AN1792/results/spaceranger/cohort5-LCMB"

id=$1
transcriptome_ref=$2
probe_set=$3
feature_ref=$4
lib=$5
cyt_img=$6
dk_img=$7
json=$8
slide=$9
area=${10}
output_dir=${11}



# Navigate to output directory
cd $output_dir

# Run spaceranger
/projects/b1169/zzhang/software/spaceranger-2.1.1/bin/spaceranger count \
--id="${id}" \
--transcriptome="${transcriptome_ref}" \
--probe-set="${probe_set}" \
--feature-ref="${feature_ref}" \
--libraries="${lib}" \
--cytaimage="${cyt_img}" \
--darkimage="${dk_img}" \
--loupe-alignment="${json}" \
--slide="${slide}" \
--area="${area}" \
--localcores=16 \
--localmem=128