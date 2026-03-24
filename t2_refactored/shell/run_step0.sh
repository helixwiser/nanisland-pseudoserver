#!/bin/bash
# run_step0.sh
# SLURM 作业脚本：数据准备（合并ID + 重叠统计）
#
# 计算量小，普通队列即可
#
# 步骤：
#   1. python step0_prepare/check_overlap.py
#   2. python step0_prepare/merge_ids.py
#
# 前置：确认 ca1_ids.txt / no1_ids.txt 存在
# 产出：all_ids_merged.txt / protein_source_labels.tsv / overlap_stats.tsv
