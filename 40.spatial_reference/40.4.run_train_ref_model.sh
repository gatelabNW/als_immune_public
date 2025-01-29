#!/bin/bash
#SBATCH --account p31535
#SBATCH --partition gengpu
#SBATCH --job-name C2LT
#SBATCH --gres=gpu:a100:1
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 1
#SBATCH --mem 60GB
#SBATCH --time 4:00:00
#SBATCH --output=/projects/b1042/Gate_Lab/zzhang/cluster_logs/%j.%N.txt
#SBATCH --error=/projects/b1042/Gate_Lab/zzhang/cluster_logs/%j.%N.err
#SBATCH --verbose


# State date and name
date
source activate cell2loc_env

cd /projects/p31535/zzhang/als/als_repo/40.spatial_reference
python3 40.4.train_ref_model.py