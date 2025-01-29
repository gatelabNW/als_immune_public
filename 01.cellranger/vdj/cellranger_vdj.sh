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
module load cellranger/6.0.0


# Navigate to output directory
cd ${4}

# Run cellranger vdj
cellranger vdj \
--id ${1} \
--fastqs ${2} \
--reference "/projects/p31535/als-project/resources/reference/refdata-cellranger-vdj-GRCh38-alts-ensembl-5.0.0" \
--sample ${3} \
--localcores 8