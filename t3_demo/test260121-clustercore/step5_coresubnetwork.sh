source /home/caishangLab/fangxiunan/anaconda3/etc/profile.d/conda.sh
conda activate metagenomics

cd /storage/caishangLab/fangxiunan/ASAP/test260121-clustercore/ca1_foldseek/
python "/storage/caishangLab/fangxiunan/ASAP/test260108-cluster2round/step5_leiden_fa.py"

cd /storage/caishangLab/fangxiunan/ASAP/test260121-clustercore/no1_foldseek/
python "/storage/caishangLab/fangxiunan/ASAP/test260108-cluster2round/step5_leiden_fa.py"
