#!/bin/bash
#SBATCH --account p31535
#SBATCH --partition gengpu
#SBATCH --job-name step2_all
#SBATCH --gres=gpu:a100:1
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 1
#SBATCH --mem 60GB
#SBATCH --time 24:00:00
#SBATCH --output=/projects/b1042/Gate_Lab/zzhang/cluster_logs/%j.%N.txt
#SBATCH --error=/projects/b1042/Gate_Lab/zzhang/cluster_logs/%j.%N.err
#SBATCH --verbose


# State date and name
date
source activate cell2loc_env

cd /projects/p31535/zzhang/als/als_repo/44.spatial_model
python3 44.1.step2.py