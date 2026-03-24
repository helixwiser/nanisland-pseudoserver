#!/bin/bash

#SBATCH -p intel-sc3,amd-ep2,intel-sc3-32c
#SBATCH -q normal
#SBATCH --cpus-per-task=32
#SBATCH --output=/storage/caishangLab/fangxiunan/log/exec_%A_%a.log
#SBATCH --error=/storage/caishangLab/fangxiunan/log/exec_%A_%a.err
#SBATCH --mem=100G

source /home/caishangLab/fangxiunan/anaconda3/etc/profile.d/conda.sh
conda activate metagenomics
#python step3_networkbone_pyvis.py
#python step3_networkbone_2d.py
cd "/storage/caishangLab/fangxiunan/ASAP/test1205_cluster_network/ca1_mmseq/"
python ../step3_test_fa.py
python ../step3_test_kk.py
cd "/storage/caishangLab/fangxiunan/ASAP/test1205_cluster_network/no1_mmseq/"
python ../step3_test_fa.py
python ../step3_test_kk.py
cd "/storage/caishangLab/fangxiunan/ASAP/test1205_cluster_network/ca1_foldseek/"
python ../step3_test_fa.py
python ../step3_test_kk.py
cd "/storage/caishangLab/fangxiunan/ASAP/test1205_cluster_network/no1_foldseek/"
python ../step3_test_fa.py
python ../step3_test_kk.py
cd "/storage/caishangLab/fangxiunan/ASAP/test1205_cluster_network/ca2_mmseq/"
python ../step3_test_fa.py
python ../step3_test_kk.py
cd "/storage/caishangLab/fangxiunan/ASAP/test1205_cluster_network/no2_mmseq/"
python ../step3_test_fa.py
python ../step3_test_kk.py
cd "/storage/caishangLab/fangxiunan/ASAP/test1205_cluster_network/ca2_foldseek/"
python ../step3_test_fa.py
python ../step3_test_kk.py
cd "/storage/caishangLab/fangxiunan/ASAP/test1205_cluster_network/no2_foldseek/"
python ../step3_test_fa.py
python ../step3_test_kk.py