#!/bin/bash
# run_step4_5.sh
# SLURM 作业脚本：拓扑分析 + 可视化
#
# 资源需求：50G 内存，16 cores
# conda activate metagenomics
#
# 步骤：
#   1. python step4_topology/compute_metrics.py
#   2. python step4_topology/detect_communities.py
#   3. python step4_topology/extract_hubs.py
#   4. python step5_visualize/annotate_uniprot.py
#   5. python step5_visualize/network_pyvis.py
#   6. python step5_visualize/network_static.py
#
# 前置：step2_3 完成
# 产出：topo / hubs / communities / HTML / PNG
