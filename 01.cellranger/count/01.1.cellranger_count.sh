#!/bin/bash
#SBATCH --account b1042
#SBATCH --partition genomicsguestA
#SBATCH --job-name cellranger_count
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 8
#SBATCH --mem 32G
#SBATCH --time 24:00:00
#SBATCH --output=/projects/b1042/Gate_Lab/zzhang/cluster_logs/%j.%N.txt
#SBATCH --error=/projects/b1042/Gate_Lab/zzhang/cluster_logs/%j.%N.err
#SBATCH --verbose

# ${1} = id
# ${2} = sample_dir
# ${3} = sample
# ${4} = cell_num
# ${5} = output_dir

date

# Module prep
module purge

# Navigate to output directory
cd ${5}

# Run 01.cellranger count
/projects/b1169/zzhang/software/cellranger-8.0.0/cellranger count \
--id ${1} \
--fastqs ${2} \
--transcriptome "/projects/b1169/projects/als-project/resources/reference/refdata-gex-GRCh38-2020-A" \
--sample ${3} \
--expect-cells ${4} \
--localcores 8 \
--create-bam true