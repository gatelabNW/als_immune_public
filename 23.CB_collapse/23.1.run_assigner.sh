#!/bin/bash
#SBATCH --account b1042
#SBATCH --partition genomicsguestA
#SBATCH --job-name ASSIGNER
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 16
#SBATCH --mem 100G
#SBATCH --time 24:00:00
#SBATCH --output=/projects/b1042/Gate_Lab/zzhang/cluster_logs/%j.%N.txt
#SBATCH --error=/projects/b1042/Gate_Lab/zzhang/cluster_logs/%j.%N.err
#SBATCH --verbose


source activate scNanoGPS

SC_NANO="/home/zzj4347/softwares/scNanoGPS/scNanoGPS"
out_root="/projects/b1042/Gate_Lab/projects/als-project/ont_seq"

in_dir="${out_root}/01.preprocess"
input_cell_barcodes="${out_root}/01.preprocess/barcode_list.tsv.gz"
out_dir="${out_root}/02.CB_collapse"

echo "Output dir is ${out_dir}"

if [ ! -d "$out_dir" ]; then
    mkdir "$out_dir"
fi
python3 ${SC_NANO}/assigner.py -i ${input_cell_barcodes} -d ${out_dir} -t 16