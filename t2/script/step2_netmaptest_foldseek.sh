#!/bin/bash

#SBATCH -p intel-sc3,amd-ep2,intel-sc3-32c
#SBATCH -q normal
#SBATCH --cpus-per-task=32
#SBATCH --output=/storage/caishangLab/fangxiunan/log/exec_%A_%a.log
#SBATCH --error=/storage/caishangLab/fangxiunan/log/exec_%A_%a.err
#SBATCH --mem=100G

source /home/caishangLab/fangxiunan/anaconda3/etc/profile.d/conda.sh
conda activate metagenomics

mkdir cluster_data
cp string_uniprot_subnet_renamed.csv cluster_data/
cp foldseek-res_cluster.tsv cluster_data/
script_path=/storage/caishangLab/fangxiunan/ASAP/test1205_cluster_network/script

python $script_path/step2_mapppis_foldseek.py
python $script_path/step2_countppis.py
python $script_path/step2_poissontest.py
python $script_path/step2_buildweightednetwork.py
python $script_path/step2_networkpyvis.py