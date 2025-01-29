#!/bin/bash
#SBATCH --account p31535
#SBATCH --partition long
#SBATCH --job-name CURATOR
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 52
#SBATCH --mem 150G
#SBATCH --time 168:00:00
#SBATCH --output=/projects/b1042/Gate_Lab/zzhang/cluster_logs/%j.%N.txt
#SBATCH --error=/projects/b1042/Gate_Lab/zzhang/cluster_logs/%j.%N.err
#SBATCH --verbose


source activate scNanoGPS

SC_NANO="/home/zzj4347/softwares/scNanoGPS/scNanoGPS"
out_root="/projects/b1042/Gate_Lab/projects/als-project/ont_seq"

in_FQ_dir="${out_root}/01.preprocess"
in_CB_dir="${out_root}/02.CB_collapse"

FQ_file="${in_FQ_dir}/processed.fastq.gz"
BC_file="${in_FQ_dir}/barcode_list.tsv.gz"
CB_count_file="${in_CB_dir}/CB_counting.tsv.gz"
CB_merged_file="${in_CB_dir}/CB_merged_list.tsv.gz"
ref_g="/projects/b1169/zzhang/ont_reference/hg38.fa.gz"
ref_g_idx="/projects/b1169/zzhang/ont_reference/ref_hg38.mmi"
out_dir="${out_root}/03.CB_alignment"
temp_out_dir="${out_root}/03.CB_alignment/tmp"

echo "Output dir is ${out_dir}"

if [ ! -d "$out_dir" ]; then
    mkdir "$out_dir"
fi

if [ ! -d "$temp_out_dir" ]; then
    mkdir "$temp_out_dir"
fi

python3 ${SC_NANO}/curator.py --fq_name ${FQ_file} -d ${out_dir} -t 64 -b ${BC_file} --CB_count ${CB_count_file} --CB_list ${CB_merged_file} --ref_genome ${ref_g} --idx_genome ${ref_g_idx} --tmp_dir=${temp_out_dir}
