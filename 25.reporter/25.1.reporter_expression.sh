#!/bin/bash
#SBATCH --account b1042
#SBATCH --partition genomicsguestA
#SBATCH --job-name 25.1
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 16
#SBATCH --mem 180G
#SBATCH --time 48:00:00
#SBATCH --output=/projects/b1042/Gate_Lab/zzhang/cluster_logs/25.1%j.%N.txt
#SBATCH --error=/projects/b1042/Gate_Lab/zzhang/cluster_logs/25.1%j.%N.err
#SBATCH --verbose


source activate scNanoGPS



cd /projects/b1042/Gate_Lab/projects/als-project/ont_seq/03.CB_alignment/B2

# TODO write this programatically later
# mkdir scNanoGPS_res


python3 /home/zzj4347/softwares/scNanoGPS/scNanoGPS/reporter_expression.py -t 16 \
                                       --gtf /projects/b1169/zzhang/ont_reference/gencode.v38.annotation.gtf \
                                       --featurecounts ~/softwares/anaconda3/envs/scNanoGPS/bin/featureCounts

