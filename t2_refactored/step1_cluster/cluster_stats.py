# cluster_stats.py
# 统计 Foldseek 聚类结果的质量指标
#
# 输入：
#   - foldseek-res_cluster.tsv: Foldseek 聚类结果
#
# 逻辑：
#   - 统计总家族数量
#   - 家族大小分布（min / max / median / mean）
#   - singleton 家族比例（只有1个成员的家族）
#   - 最大的 N 个家族及其成员数
#   - 生成家族大小分布直方图
#
# 产出：
#   - cluster_stats.tsv: 摘要统计表
#   - cluster_size_distribution.png: 分布图
#
# 假说1.1 判定依据：
#   - 家族数量在数百~数千级别 → 合理
#   - singleton 比例过高（>80%）→ 聚类太严，考虑降低 coverage
#   - 少数家族成员数极多（>1000）→ 聚类太松，考虑提高 coverage
