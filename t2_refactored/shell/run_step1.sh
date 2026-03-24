#!/bin/bash
# run_step1.sh
# SLURM 作业脚本：Foldseek 聚类
#
# 资源需求：100G 内存，20 cores
# conda activate foldseek-env
#
# 步骤：
#   1. 从 all_ids_merged.txt 生成 PDB ID 列表
#   2. 从 AlphaFold DB 提取 PDB 文件
#   3. foldseek easy-cluster，coverage=0.5
#   4. python step1_cluster/cluster_stats.py
#   5. python step1_cluster/build_family_map.py
#
# 前置：step0 完成
# 产出：foldseek-res_cluster.tsv / protein_to_family_labeled.tsv / cluster_stats.tsv
#
# ⚠️ 假说1.1 判定点：看 cluster_stats 结果再决定是否继续
