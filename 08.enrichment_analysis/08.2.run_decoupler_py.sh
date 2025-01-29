#!/bin/bash
#SBATCH --account b1042
#SBATCH --partition genomicsguestA
#SBATCH --job-name DCPPY
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 20
#SBATCH --mem 60G
#SBATCH --time 12:00:00
#SBATCH --output=/projects/b1042/Gate_Lab/zzhang/cluster_logs/DCPPY_%j.%N.txt
#SBATCH --error=/projects/b1042/Gate_Lab/zzhang/cluster_logs/DCPPY_%j.%N.err
#SBATCH --verbose

date

# Module prep
module purge
source activate decoupler

cd /projects/p31535/zzhang/als/als_repo/08.enrichment_analysis
python3 08.2.decoupler_py.py