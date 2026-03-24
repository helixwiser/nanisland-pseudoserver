#!/bin/bash

#SBATCH -p intel-sc3,amd-ep2,intel-sc3-32c
#SBATCH -q normal
#SBATCH --cpus-per-task=20
#SBATCH --output=/storage/caishangLab/fangxiunan/log/exec_%A_%a.log
#SBATCH --error=/storage/caishangLab/fangxiunan/log/exec_%A_%a.err
#SBATCH --mem=100G

source /home/caishangLab/fangxiunan/anaconda3/etc/profile.d/conda.sh
script_path=/storage/caishangLab/fangxiunan/ASAP/test1202-cluster/script
conda activate metagenomics
python $script_path/step1_getstringdbsubnet.py
python $script_path/step1_networksummary.py