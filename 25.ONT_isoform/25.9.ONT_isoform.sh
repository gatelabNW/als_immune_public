#!/bin/bash
#SBATCH --account b1042
#SBATCH --partition genomicsguestA
#SBATCH --job-name ISO_DET
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 52
#SBATCH --mem 200G
#SBATCH --time 48:00:00
#SBATCH --output=/projects/b1042/Gate_Lab/zzhang/cluster_logs/%j.%N.txt
#SBATCH --error=/projects/b1042/Gate_Lab/zzhang/cluster_logs/%j.%N.err
#SBATCH --verbose


source activate scNanoGPS

SC_NANO="/home/zzj4347/softwares/scNanoGPS/scNanoGPS"
out_root="/projects/b1042/Gate_Lab/projects/als-project/ont_seq"
cur_sample="B2"
out_dir="${out_root}/04.novel_isoforms/${cur_sample}"
REF_GENE="/projects/b1169/zzhang/ont_reference/gencode.v38.chr_patch_hapl_scaff.annotation.refgene"
temp_out_dir="${out_root}/03.CB_alignment_no_consensus/03.CB_alignment_no_consensus/${cur_sample}/tmp"

echo "Output dir is ${out_dir}"

if [ ! -d "${out_dir}" ]; then
    mkdir -p "${out_dir}"
    echo "Created ${out_dir}"
fi

if [ ! -d "${temp_out_dir}" ]; then
    mkdir -p "${temp_out_dir}"
    echo "Created ${temp_out_dir}"
fi

python3 ${SC_NANO}/reporter_isoform.py -d ${out_dir} -t 52 --liqa_ref ${REF_GENE} --tmp_dir ${temp_out_dir}