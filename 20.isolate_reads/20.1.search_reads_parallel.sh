#!/bin/bash
#SBATCH --account p31535
#SBATCH --partition short
#SBATCH --job-name ST_RD
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 52
#SBATCH --mem 16G
#SBATCH --time 1:00:00
#SBATCH --output=/projects/b1042/Gate_Lab/zzhang/cluster_logs/%j.%N.txt
#SBATCH --error=/projects/b1042/Gate_Lab/zzhang/cluster_logs/%j.%N.err
#SBATCH --verbose
# Directory containing fastq.gz files
FASTQ_DIR="/projects/b1042/Gate_Lab/projects/als-project/ont_seq/00.raw/D3/20230922_1539_P2S-00717-A_PAS47180_d87f833c/fastq_pass"

# Output directory where the merged fastq.gz will be saved
OUTPUT_DIR="/projects/b1042/Gate_Lab/projects/als-project/ont_seq/01.isolated_reads/D3/reads"

# Temporary directory for intermediate files
TEMP_DIR="/projects/b1042/Gate_Lab/projects/als-project/ont_seq/01.isolated_reads/D3/tmp"

# Ensure the provided directories exist
if [ ! -d "$FASTQ_DIR" ]; then
    echo "Error: Input directory '$FASTQ_DIR' does not exist."
    exit 1
fi

if [ ! -d "$OUTPUT_DIR" ]; then
    echo "Error: Output directory '$OUTPUT_DIR' does not exist."
    mkdir -p "$OUTPUT_DIR"
fi

if [ ! -d "$TEMP_DIR" ]; then
    echo "Creating temporary directory '$TEMP_DIR'."
    mkdir -p "$TEMP_DIR"
fi

# Find all fastq.gz files in the provided directory
FASTQ_FILES=("$FASTQ_DIR"/*.fastq.gz)

# Array to hold the PIDs of the background processes
PIDS=()

# Total number of files
TOTAL_FILES=${#FASTQ_FILES[@]}
# Counter for progress
COUNT=0

# Function to search for the pattern in a single file
search_reads() {
  local file=$1
  local temp=$2
  local outfile="$temp"/$(basename "$file" .fastq.gz)_GGGGCC.fastq
  # Use awk to extract reads containing the sequence and count the occurrences
  zcat "$file" | awk -v outfile="$outfile" '
  BEGIN {OFS = "\n"}
  {
    header = $0;
    getline seq;
    getline sep;
    getline qual;
    count = gsub(/GGGGCC/, "&", seq); # Count occurrences of GGGGCC in the sequence
    if(count > 10) {
      print header, seq, sep, qual > outfile;
      print count > outfile "_counts.txt"; # Output the count to a separate file
    }
  }'
  # Update progress
  echo -ne "Progress: $COUNT/$TOTAL_FILES files processed.\r"
}

# Export the function so it can be used by background processes
export -f search_reads

# Loop over files and search for reads containing "GGGGCC" in the background
for file in "${FASTQ_FILES[@]}"; do
  search_reads "$file" "$TEMP_DIR" &
  # Save the PID of the background process
  PIDS+=($!)
  # We need to limit the number of background processes to avoid overloading the system
  # This will wait until we have less than the number of CPU cores processes running
  while [ $(jobs -r | wc -l) -ge 52 ]; do
    sleep 1
  done
done

# Wait for all background processes to finish
for pid in "${PIDS[@]}"; do
  wait $pid
done

echo -ne "\nAll files processed. Merging results...\n"

# Concatenate all the search results into a single fastq file
cat "$TEMP_DIR"/*_GGGGCC.fastq > "$TEMP_DIR"/merged_GGGGCC_reads.fastq

# Compress the concatenated file into fastq.gz
gzip -c "$TEMP_DIR"/merged_GGGGCC_reads.fastq > "$OUTPUT_DIR"/merged_GGGGCC_reads.fastq.gz
cat "$TEMP_DIR"/*_counts.txt > "$OUTPUT_DIR"/merged_counts.txt


# Clean up the temporary directory and files, if it was created by this script
#if [ ! -d "$3" ]; then
#  rm -r "$TEMP_DIR"
#fi

echo "Merging completed. Compressed file is ready for alignment."

