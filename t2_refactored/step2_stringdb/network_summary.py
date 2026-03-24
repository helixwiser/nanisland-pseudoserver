# network_summary.py
# StringDB 子网络基础统计
#
# 输入：
#   - string_uniprot_subnet.csv: STRING ID 子网络
#   - all_ids_merged.txt: 全量 UniProt ID
#
# 逻辑：
#   - 统计子网络中每个 UniProt ID 的边数（degree）
#   - 用 pd.melt + value_counts 向量化实现，替代旧版 for 循环
#
# 产出：
#   - network_summary.csv: 两列（id / edge_count）
#
# 对应旧代码：step1_networksummary.py
# 改进：旧版逐个 ID for 循环，大网络性能差。改为向量化操作。
