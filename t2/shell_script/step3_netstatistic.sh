#!/bin/bash

#SBATCH -p intel-sc3,amd-ep2,intel-sc3-32c
#SBATCH -q normal
#SBATCH --cpus-per-task=16
#SBATCH --output=/storage/caishangLab/fangxiunan/log/exec_%A_%a.log
#SBATCH --error=/storage/caishangLab/fangxiunan/log/exec_%A_%a.err
#SBATCH --mem=50G

source /home/caishangLab/fangxiunan/anaconda3/etc/profile.d/conda.sh
conda activate metagenomics


process_dir=/storage/caishangLab/fangxiunan/ASAP/test1223_proteincluster_foldseek
script_path=$process_dir/script


cd $process_dir/ca1_foldseek
cp $process_dir/merged_details_unique.tsv .

python $script_path/step3_test_fa.py
#python $script_path/step3_test_kk.py
#
##
#cd $process_dir/no1_foldseek
#cp $process_dir/merged_details_unique.tsv .
#python $script_path/step3_test_fa.py
#python $script_path/step3_test_kk.py
#
#
#cd $process_dir/ca2_foldseek
#cp $process_dir/merged_details_unique.tsv .
#python $script_path/step3_test_fa.py
#python $script_path/step3_test_kk.py
#
#cd $process_dir/no2_foldseek
#cp $process_dir/merged_details_unique.tsv .
#python $script_path/step3_test_fa.py
#python $script_path/step3_test_kk.py



#python step3_networkbone_pyvis.py
#python step3_networkbone_2d.py
#cd "/storage/caishangLab/fangxiunan/ASAP/test1205_cluster_network/ca1_mmseq/"
#python ../step3_test_fa.py
#python ../step3_test_kk.py
#cd "/storage/caishangLab/fangxiunan/ASAP/test1205_cluster_network/no1_mmseq/"
#python ../step3_test_fa.py
#python ../step3_test_kk.py
#cd "/storage/caishangLab/fangxiunan/ASAP/test1205_cluster_network/ca1_foldseek/"
#python ../step3_test_fa.py
#python ../step3_test_kk.py
#cd "/storage/caishangLab/fangxiunan/ASAP/test1205_cluster_network/no1_foldseek/"
#python ../step3_test_fa.py
#python ../step3_test_kk.py
#cd "/storage/caishangLab/fangxiunan/ASAP/test1205_cluster_network/ca2_mmseq/"
#python ../step3_test_fa.py
#python ../step3_test_kk.py
#cd "/storage/caishangLab/fangxiunan/ASAP/test1205_cluster_network/no2_mmseq/"
#python ../step3_test_fa.py
#python ../step3_test_kk.py
#cd "/storage/caishangLab/fangxiunan/ASAP/test1205_cluster_network/ca2_foldseek/"
#python ../step3_test_fa.py
#python ../step3_test_kk.py
#cd "/storage/caishangLab/fangxiunan/ASAP/test1205_cluster_network/no2_foldseek/"
#python ../step3_test_fa.py
#python ../step3_test_kk.py