#!/bin/bash
#SBATCH --account b1042
#SBATCH --partition genomicsburst
#SBATCH --job-name 25.3
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 52
#SBATCH --mem 2000G
#SBATCH --time 240:00:00
#SBATCH --output=/projects/b1042/Gate_Lab/zzhang/cluster_logs/%j.%N.txt
#SBATCH --error=/projects/b1042/Gate_Lab/zzhang/cluster_logs/%j.%N.err
#SBATCH --verbose

source activate isoquant

ref_gtf="/projects/b1169/zzhang/ont_reference/gencode.v46.chr_patch_hapl_scaff.annotation.gtf.gz"
ref_g="/projects/b1169/zzhang/ont_reference/hg38.fa.gz"

yaml="/projects/b1042/Gate_Lab/projects/als-project/ont_seq/04.isoquant/02.ymal/A4.yaml"
output="/projects/b1042/Gate_Lab/projects/als-project/ont_seq/04.isoquant/03.isoquant_out/A4"

if [ ! -d "$output" ]; then
    mkdir "$output"
fi

isoquant.py -d nanopore --yaml ${yaml} --complete_genedb --genedb ${ref_gtf} --reference ${ref_g} --output ${output} -t 52