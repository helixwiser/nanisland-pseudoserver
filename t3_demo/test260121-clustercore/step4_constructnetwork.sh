#!/bin/bash

#SBATCH -p intel-sc3,amd-ep2,intel-sc3-32c
#SBATCH -q normal
#SBATCH --cpus-per-task=20
#SBATCH --output=/storage/caishangLab/fangxiunan/log/exec_%A_%a.log
#SBATCH --error=/storage/caishangLab/fangxiunan/log/exec_%A_%a.err
#SBATCH --mem=100G

source /home/caishangLab/fangxiunan/anaconda3/etc/profile.d/conda.sh
conda activate metagenomics

process_dir=/storage/caishangLab/fangxiunan/ASAP/test260121-clustercore
script_path=/storage/caishangLab/fangxiunan/ASAP/test1223_proteincluster_foldseek/script
###################

cd $process_dir/ca1_foldseek

cp $process_dir/ca1_allids_prokaryotes.txt UniProt_ids.txt
python $script_path/step1_getstringdbsubnet.py
#python $script_path/step1_networksummary.py

mkdir cluster_data
cp string_uniprot_subnet_renamed.csv cluster_data/
cp foldseek-res_c0.5_cluster.tsv cluster_data/foldseek-res_cluster.tsv


python $script_path/step2_mapppis_foldseek.py
python $script_path/step2_countppis.py
python $script_path/step2_poissontest.py
python $script_path/step2_buildweightednetwork.py
python $script_path/step2_networkpyvis.py


###################
cd $process_dir/no1_foldseek

cp $process_dir/no1_allids_prokaryotes.txt UniProt_ids.txt
python $script_path/step1_getstringdbsubnet.py
#python $script_path/step1_networksummary.py

mkdir cluster_data
cp string_uniprot_subnet_renamed.csv cluster_data/
cp foldseek-res_c0.5_cluster.tsv cluster_data/foldseek-res_cluster.tsv


python $script_path/step2_mapppis_foldseek.py
python $script_path/step2_countppis.py
python $script_path/step2_poissontest.py
python $script_path/step2_buildweightednetwork.py
python $script_path/step2_networkpyvis.py