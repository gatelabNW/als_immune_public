####################################################################
# Natalie Piehl, Charles Zhang
# Gate Lab
# Northwestern University
# 02/08/2023
####################################################################
# Run 01.cellranger vdj on TCR fastq files from ALS project
# Details from: https://support.10xgenomics.com/single-cell-vdj/software/pipelines/latest/output/annotation
# Expected output:
#   clonotypes.csv
#       High-level descriptions of each clonotype.
#   consensus_annotations.csv
#       High-level and detailed annotations of each clonotype consensus sequence.
#   filtered_contig_annotations.csv
#       High-level annotations of each high-confidence, cellular contig. This is a subset of all_contig_annotations.csv.
#   all_contig_annotations.{csv,bed,json}
#       High-level and detailed annotations of each contig.
#   airr_rearrangement.tsv
#       Annotated contigs and consensus sequences of VDJ rearrangements in the AIRR format.
####################################################################
# INSTRUCTIONS FOR USE:
# Change fastq_dir to the directory holding fastq files to process
# Change output_dir to directory you want to direct output
# Execute "bash batch_cellranger_tcr_vdj.sh"
####################################################################

date

# Define and check directory containing GEX fastq files
fastq_dir="/projects/b1042/Gate_Lab/projects/als-project/scrna_original_final/00.raw"
[ -d "$fastq_dir" ] && echo "Directory to fastq files exists."

# Define, make (if necessary), and navigate to output directory
output_dir="/projects/b1042/Gate_Lab/projects/als-project/scrna_original_final/01.cellranger_tcr_vdj"
[ -d "$output_dir" ] || mkdir "$output_dir"


find "$fastq_dir" -maxdepth 1 -mindepth 1 -type d -name "*TCR*" -print0 | while read -d $'\0' sample_dir
do
  id=$(basename "$sample_dir")
  id=${id%"TCRALS"}
  # Check if output directory already exists, and if so skip to next sample
#  if [ -d "$output_dir/$id" ]; then
#    continue
#  fi

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
  # Run 01.cellranger vdj script
  sbatch 01.1.cellranger_vdj.sh "$id" "$sample_dir" "$sample" "$output_dir"
  echo "cellranger_vdj.sh" "$id" "$sample_dir" "$sample" "$output_dir"
done
