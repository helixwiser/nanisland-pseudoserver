#!/bin/bash
# run_step2_3.sh
# SLURM 作业脚本：StringDB 子网络提取 + 家族级网络构建
#
# 资源需求：100G 内存（StringDB links 文件大），32 cores
# conda activate metagenomics
#
# 步骤：
#   1. python step2_stringdb/extract_subnet.py
#   2. python step2_stringdb/network_summary.py（可选）
#   3. python step3_network/map_ppi_to_family.py
#   4. python step3_network/count_interactions.py
#   5. python step3_network/poisson_test.py
#   6. python step3_network/build_network.py
#
# 前置：step1 完成 + 假说1.1 判定通过
# 产出：network_ca.csv / network_no.csv
