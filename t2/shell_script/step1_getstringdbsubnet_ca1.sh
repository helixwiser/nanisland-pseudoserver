#!/bin/bash

#SBATCH -p intel-sc3,amd-ep2,intel-sc3-32c
#SBATCH -q normal
#SBATCH --cpus-per-task=20
#SBATCH --output=/storage/caishangLab/fangxiunan/log/exec_%A_%a.log
#SBATCH --error=/storage/caishangLab/fangxiunan/log/exec_%A_%a.err
#SBATCH --mem=100G

source /home/caishangLab/fangxiunan/anaconda3/etc/profile.d/conda.sh
conda activate metagenomics


process_dir=/storage/caishangLab/fangxiunan/ASAP/test1223_proteincluster_foldseek
script_path=$process_dir/script


cd $process_dir/ca1_foldseek
cp $process_dir/ca1_allids_prokaryotes.txt UniProt_ids.txt
python $script_path/step1_getstringdbsubnet.py
python $script_path/step1_networksummary.py