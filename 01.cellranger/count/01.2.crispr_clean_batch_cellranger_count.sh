####################################################################
# Natalie Piehl, Charles Zhang
# Gate Lab
# Northwestern University
# Date: 09-14-2023
####################################################################
# Run 01.cellranger count on GEX fastq files from ALS project
# Details from: https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/output/overview
# Expected output:
#    Run summary HTML: web_summary.html
#    Run summary CSV: metrics_summary.csv
#    BAM: possorted_genome_bam.bam
#    BAM index: possorted_genome_bam.bam.bai
#    Filtered feature-barcode matrices MEX: filtered_feature_bc_matrix
#    Filtered feature-barcode matrices HDF5: filtered_feature_bc_matrix.h5
#    Unfiltered feature-barcode matrices MEX: raw_feature_bc_matrix
#    Unfiltered feature-barcode matrices HDF5: raw_feature_bc_matrix_h5.h5
#    Secondary analysis output CSV: analysis
#    Per-molecule read information: molecule_info.h5
#    Loupe Browser file: cloupe.cloupe
####################################################################
# INSTRUCTIONS FOR USE:
# Change fastq_dir to the directory holding fastq files to process
# Change output_dir to directory you want to direct output
# Change csv to path to a file containing cell number per sample
#   Note: For the ALS project, we did not provide a cell count file
# Execute "bash batch_cellranger_count.sh"
####################################################################

date

# Define and check directory containing GEX fastq files
fastq_dir="/projects/b1042/Gate_Lab/projects/als-project/crispr_clean_final/00.raw"
[ -d "$fastq_dir" ] && echo "Directory to fastq files exists."

# Define, make (if necessary), and navigate to output directory
output_dir="/projects/b1042/Gate_Lab/projects/als-project/crispr_clean_final/01.cellranger_count"
[ -d "$output_dir" ] || mkdir "$output_dir"

# Create array of fastq_dir subdirectories
sample_arr=()
while IFS=  read -r -d $'\0'; do
    sample_arr+=("$REPLY")
done < <(find "$fastq_dir" -maxdepth 1 -mindepth 1 -type d -print0)


# Perform 01.cellranger count on each sample
for sample_dir in "${sample_arr[@]}"
do
    # Isolate sample ID
    id=$(basename "$sample_dir")

    # Get names of files in sample_dir
    files=()
    while IFS=  read -r -d $'\0'; do
        files+=("$REPLY")
    done < <(find "$sample_dir" -type f -name "*.fastq.gz" -print0)

    # Identify if need to cutoff last 24 or 25 characters
    sample_num=$(echo "$files[0]" | grep -o -P '(?<=_S).*(?=_L)')
    if [ ${#sample_num} -eq 2 ]; then
      cutoff=25; else
      cutoff=24
    fi

    # Isolate sample name
    sample=${files[0]::-$cutoff}
    sample=$(basename "$sample")
    echo --------------------------------------------------------------------------------------------
    echo "Submitting $id: $sample now"

    cell_num=10000

    # Run 01.cellranger count script
    sbatch 01.1.cellranger_count.sh "$id" "$sample_dir" "$sample" "$cell_num" "$output_dir"
    echo $id,$sample,$sample_dir,$output_dir
    sleep 1
done