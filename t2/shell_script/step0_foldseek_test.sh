#!/bin/bash

#SBATCH -p intel-sc3,amd-ep2,intel-sc3-32c
#SBATCH -q normal
#SBATCH --cpus-per-task=20
#SBATCH --output=/storage/caishangLab/fangxiunan/log/exec_%A_%a.log
#SBATCH --error=/storage/caishangLab/fangxiunan/log/exec_%A_%a.err
#SBATCH --mem=100G


source /home/caishangLab/fangxiunan/anaconda3/etc/profile.d/conda.sh
###########################################################################
cd /storage/caishangLab/fangxiunan/ASAP/test1223_proteincluster_foldseek/
conda activate foldseek-env
# 回到工作目录（根据你的路径调整）
cd ca1_foldseek

# 定义要测试的 -c 值列表
coverage_thresholds=(0.3 0.4 0.5 0.6 0.7 0.8 0.9)


# 为每个 -c 值创建独立文件夹并运行 foldseek
for c in "${coverage_thresholds[@]}"; do
    # 创建以 c 值命名的输出目录
    out_dir="foldseek-res_c${c}"
    tmp_dir="tmp_c${c}"
    mkdir -p "$out_dir" "$tmp_dir"


    echo "Running foldseek with -c $c ..."
    foldseek easy-cluster sw_sub/ "$out_dir" "$tmp_dir" -c "$c"


    echo "Result saved in: $out_dir"
done


echo "All foldseek clustering tasks completed."
#######################################################################
cd /storage/caishangLab/fangxiunan/ASAP/test1223_proteincluster_foldseek/
cd no1_foldseek

# 定义要测试的 -c 值列表
coverage_thresholds=(0.3 0.4 0.5 0.6 0.7 0.8 0.9)


# 为每个 -c 值创建独立文件夹并运行 foldseek
for c in "${coverage_thresholds[@]}"; do
    # 创建以 c 值命名的输出目录
    out_dir="foldseek-res_c${c}"
    tmp_dir="tmp_c${c}"
    mkdir -p "$out_dir" "$tmp_dir"


    echo "Running foldseek with -c $c ..."
    foldseek easy-cluster sw_sub/ "$out_dir" "$tmp_dir" -c "$c"


    echo "Result saved in: $out_dir"
done


echo "All foldseek clustering tasks completed."
#mkdir ca1_foldseek
#sed 's/.*/AF-&-F1-model_v6.pdb.gz/' ca1_allids_prokaryotes.txt > ca1_foldseek/ca1_pdb_ids.txt
#cd ca1_foldseek
#tar -xf "/storage/caishangLab/fangxiunan/database/ALPHAFOLD_DB/swissprot_pdb_v6.test.tar" -T ca1_pdb_ids.txt
#mkdir -p sw_sub
#find . -name "AF*" -maxdepth 1 -print0 | xargs -0 mv -t sw_sub/
#foldseek easy-cluster sw_sub/ foldseek-res tmp -c 0.7 
#
#
##
##cd /storage/caishangLab/fangxiunan/ASAP/test1202-cluster/
##conda activate mmseq2-env 
### 1. 创建序列数据库
##mkdir ca1_mmseq
##cd ca1_mmseq
##mmseqs createdb "/storage/caishangLab/fangxiunan/ASAP/test1202-cluster/ca1_protein_all.fasta" inDB
##mmseqs cluster inDB outDB tmp
### 3. 转换聚类结果为 TSV
##mmseqs createtsv inDB inDB outDB cluster_results.tsv
#
############################################################################
#cd /storage/caishangLab/fangxiunan/ASAP/test1223_proteincluster_foldseek/
#conda activate  foldseek-env
#
#mkdir ca2_foldseek
#sed 's/.*/AF-&-F1-model_v6.pdb.gz/' ca2_allids_prokaryotes.txt > ca2_foldseek/ca2_pdb_ids.txt
#cd ca2_foldseek
#tar -xf "/storage/caishangLab/fangxiunan/database/ALPHAFOLD_DB/swissprot_pdb_v6.test.tar" -T ca2_pdb_ids.txt
#mkdir -p sw_sub
#find . -name "AF*" -maxdepth 1 -print0 | xargs -0 mv -t sw_sub/
#foldseek easy-cluster sw_sub/ foldseek-res tmp -c 0.7 
#
#
##
##cd /storage/caishangLab/fangxiunan/ASAP/test1202-cluster/
##conda activate mmseq2-env 
### 1. 创建序列数据库
##mkdir ca2_mmseq
##cd ca2_mmseq
##mmseqs createdb "/storage/caishangLab/fangxiunan/ASAP/test1202-cluster/ca2_protein_all.fasta" inDB
##mmseqs cluster inDB outDB tmp
### 3. 转换聚类结果为 TSV
##mmseqs createtsv inDB inDB outDB cluster_results.tsv
#
############################################################################
#cd /storage/caishangLab/fangxiunan/ASAP/test1223_proteincluster_foldseek/
#conda activate  foldseek-env
#
#mkdir no2_foldseek
#sed 's/.*/AF-&-F1-model_v6.pdb.gz/' no2_allids_prokaryotes.txt > no2_foldseek/no2_pdb_ids.txt
#cd no2_foldseek
#tar -xf "/storage/caishangLab/fangxiunan/database/ALPHAFOLD_DB/swissprot_pdb_v6.test.tar" -T no2_pdb_ids.txt
#mkdir -p sw_sub
#find . -name "AF*" -maxdepth 1 -print0 | xargs -0 mv -t sw_sub/
#foldseek easy-cluster sw_sub/ foldseek-res tmp -c 0.7 
#
#
##
##cd /storage/caishangLab/fangxiunan/ASAP/test1202-cluster/
##conda activate mmseq2-env 
### 1. 创建序列数据库
##mkdir no2_mmseq
##cd no2_mmseq
##mmseqs createdb "/storage/caishangLab/fangxiunan/ASAP/test1202-cluster/no2_protein_all.fasta" inDB
##mmseqs cluster inDB outDB tmp
### 3. 转换聚类结果为 TSV
##mmseqs createtsv inDB inDB outDB cluster_results.tsv
#
#
############################################################################
#cd /storage/caishangLab/fangxiunan/ASAP/test1223_proteincluster_foldseek/
#conda activate  foldseek-env
#
#mkdir no1_foldseek
#sed 's/.*/AF-&-F1-model_v6.pdb.gz/' no1_allids_prokaryotes.txt > no1_foldseek/no1_pdb_ids.txt
#cd no1_foldseek
#tar -xf "/storage/caishangLab/fangxiunan/database/ALPHAFOLD_DB/swissprot_pdb_v6.test.tar" -T no1_pdb_ids.txt
#mkdir -p sw_sub
#find . -name "AF*" -maxdepth 1 -print0 | xargs -0 mv -t sw_sub/
#foldseek easy-cluster sw_sub/ foldseek-res tmp -c 0.7 
#
#
#
##cd /storage/caishangLab/fangxiunan/ASAP/test1202-cluster/
##conda activate mmseq2-env 
### 1. 创建序列数据库
##mkdir no1_mmseq
##cd no1_mmseq
##mmseqs createdb "/storage/caishangLab/fangxiunan/ASAP/test1202-cluster/no1_protein_all.fasta" inDB
##mmseqs cluster inDB outDB tmp
### 3. 转换聚类结果为 TSV
##mmseqs createtsv inDB inDB outDB cluster_results.tsv
#
## 4. 提取代表序列
##mmseqs extractclusters inDB outDB rep_seq_db
##mmseqs convert2fasta rep_seq_db rep_sequences.fasta
#
#
#
##
##mmseqs easy-cluster "/storage/caishangLab/fangxiunan/ASAP/test0912-clustering/sequences.fasta" \
##clusterRes tmp --min-seq-id 0.1 -c 0.1 --cov-mode 1
#
#
#
## 1. 使用Pfam扫描结构域（MSCRAMM特化参数）
##conda activate hmmer-env   
##hmmscan --tblout initial_domains.txt "/storage/caishangLab/fangxiunan/database/PFAM/Pfam-A.hmm" "/storage/caishangLab/fangxiunan/ASAP/test0912-clustering/sequences.fasta"
##hmmscan --domtblout pro_domains.txt --cut_ga  "/storage/caishangLab/fangxiunan/database/PFAM/Pfam-A.hmm"  "/storage/caishangLab/fangxiunan/ASAP/test0912-clustering/sequences.fasta"
#
#
