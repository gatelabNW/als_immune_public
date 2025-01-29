# Charles Zhang
# Gate Lab
# Northwestern University
# 10/03/2023
####################################################################
# Extract known genes with CE in ALS patients and HC controls
# Then merge them into a bams per condition
# Genes of interest are:
#    - STMN2: chr8:79611117-79666158
#    - KALRN: chr3:124033369-124726325
#    - UNC13A: chr19:17601336-17688354
#    - C9orf72: chr9:27546546-27573481

# declare dictionary
declare -A gene_dictionary
gene_dictionary["STMN2"]="chr8:79611117-79666158"
gene_dictionary["KALRN"]="chr3:124033369-124726325"
gene_dictionary["UNC13A"]="chr19:17601336-17688354"
gene_dictionary["C9orf72"]="chr9:27546546-27573481"

# assays

assays=("crispr_clean_final" "immune_enriched_final" "scrna_original_final")
for cur_assay in "${assays[@]}"
do
  # specify input dir
  input_root_dir="/projects/b1042/Gate_Lab/projects/als-project/${cur_assay}/01.cellranger_count"
  output_dir="/projects/b1042/Gate_Lab/projects/als-project/${cur_assay}/20.cryptic_exon/21.bams_of_interest"
  if [ ! -d "$output_dir" ]; then
      mkdir -p "$output_dir"
  fi

  # define patients and HC grouping
  HC_samples=("B5" "B4" "H4" "E1" "G5" "C4" "A4" "E4" "F2" "A5" "A3" "E3" "D2" "D1" "A2" "B2" "F1" "H1")
  ALS_samples=("H2" "G3" "E5" "B3" "H3" "A1" "G1" "C5" "E2" "F5" "D5" "D4" "F4" "G4" "C3" "G2" "B1" "C1" "F3" "C2" "H5" "D3")
  ALS_out_dir="/projects/b1042/Gate_Lab/projects/als-project/${cur_assay}/20.cryptic_exon/21.bams_of_interest/ALS"
  HC_out_dir="/projects/b1042/Gate_Lab/projects/als-project/${cur_assay}/20.cryptic_exon/21.bams_of_interest/HC"

  if [ ! -d "$ALS_out_dir" ]; then
      mkdir "$ALS_out_dir"
  fi
  if [ ! -d "$HC_out_dir" ]; then
      mkdir "$HC_out_dir"
  fi

  # Do ALS first
  echo "INFO: Working on ALS samples!"
  cd ${ALS_out_dir}
  for cur_sample in "${ALS_samples[@]}"
  do
    echo "$cur_sample"
    for cur_gene in "${!gene_dictionary[@]}"; do
        echo "${cur_gene}: ${gene_dictionary[$cur_gene]}"
        cur_bam="${input_root_dir}/${cur_sample}/outs/possorted_genome_bam.bam"
        out_sam="${ALS_out_dir}/${cur_gene}_${cur_sample}.sam"
        out_bam="${ALS_out_dir}/${cur_gene}_$cur_sample.bam"
        out_sorted_bam="${ALS_out_dir}/${cur_gene}_${cur_sample}_sorted.bam"
        samtools view -h "${cur_bam}" "${gene_dictionary[$cur_gene]}" > "${out_sam}"
        samtools view -bhS "${out_sam}" > "${out_bam}"
        samtools sort "${out_bam}" -o "${out_sorted_bam}"
    done
  done

  # merge
  for cur_gene in "${!gene_dictionary[@]}"; do
    samtools merge -f ALS_${cur_gene}.bam *${cur_gene}*sorted.bam
    samtools sort "ALS_${cur_gene}.bam" -o "ALS_${cur_gene}_sorted.bam"
    samtools index "ALS_${cur_gene}_sorted.bam"
  done


  # HC Second
  echo "INFO: Working on HC samples!"
  cd ${HC_out_dir}
  for cur_sample in "${HC_samples[@]}"
  do
    echo "$cur_sample"
    for cur_gene in "${!gene_dictionary[@]}"; do
        echo "${cur_gene}: ${gene_dictionary[$cur_gene]}"
        cur_bam="${input_root_dir}/${cur_sample}/outs/possorted_genome_bam.bam"
        out_sam="${HC_out_dir}/${cur_gene}_${cur_sample}.sam"
        out_bam="${HC_out_dir}/${cur_gene}_$cur_sample.bam"
        out_sorted_bam="${HC_out_dir}/${cur_gene}_${cur_sample}_sorted.bam"
        samtools view -h "${cur_bam}" "${gene_dictionary[$cur_gene]}" > "${out_sam}"
        samtools view -bhS "${out_sam}" > "${out_bam}"
        samtools sort "${out_bam}" -o "${out_sorted_bam}"
    done
  done

  # merge
  for cur_gene in "${!gene_dictionary[@]}"; do
    samtools merge -f HC_${cur_gene}.bam *${cur_gene}*sorted.bam
    samtools sort "HC_${cur_gene}.bam" -o "HC_${cur_gene}_sorted.bam"
    samtools index "HC_${cur_gene}_sorted.bam"
  done

done

