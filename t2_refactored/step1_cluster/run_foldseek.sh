#!/bin/bash
# run_foldseek.sh
# 用 Foldseek 对合并后的全量蛋白做结构聚类
#
# 前置：
#   - step0 产出的 all_ids_merged.txt
#   - AlphaFold DB 中对应的 PDB 文件
#
# 步骤：
#   1. 从 all_ids_merged.txt 生成 PDB ID 列表（AF-{UniProtID}-F1-model_v6.pdb.gz）
#   2. 从 AlphaFold DB tar 包中提取对应 PDB 文件到 sw_sub/
#   3. 运行 foldseek easy-cluster sw_sub/ output tmp -c 0.5
#
# 产出：
#   - foldseek-res_cluster.tsv: 聚类结果（rep, member）
#   - foldseek-res_rep_seq.fasta: 代表序列
#
# SLURM 参数参考：
#   - 内存 100G
#   - CPU 20 cores
#   - conda activate foldseek-env
