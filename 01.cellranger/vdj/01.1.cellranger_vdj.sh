#!/bin/bash
#SBATCH --account p31535
#SBATCH --partition normal
#SBATCH --job-name cellranger_vdj
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 8
#SBATCH --mem 32G
#SBATCH --time 16:00:00
#SBATCH --output /projects/b1042/Gate_Lab/zzhang/cluster_logs/%j.%N.txt
#SBATCH --verbose

# ${1} = id
# ${2} = sample_dir
# ${3} = sample
# ${4} = output_dir

date

# Module prep
module purge


# Navigate to output directory
cd ${4}

# Run cellranger vdj
/projects/b1169/zzhang/software/cellranger-8.0.0/cellranger vdj \
--id ${1} \
--fastqs ${2} \
--reference "/projects/b1169/zzhang/cellranger_ref/refdata-cellranger-vdj-GRCh38-alts-ensembl-7.1.0" \
--sample ${3} \
--localcores 8