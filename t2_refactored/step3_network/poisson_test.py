# poisson_test.py
# 分组 Poisson 富集检验，判断家族对之间的互作是否显著高于随机期望
#
# 输入：
#   - family_interaction_count_labeled.tsv: 分来源计数表
#   - protein_to_family_labeled.tsv: 用于计算家族 degree
#
# 逻辑：
#   对 ca / no 各做一次：
#   1. 计算每个家族的 degree（在对应组的家族级 PPI 网络中的连接数）
#   2. 配置模型期望 E = (d1 × d2) / total_edges
#      同家族自环：E = d1 × (d1-1) / (2 × total_edges)
#   3. Poisson 检验：P(X >= observed | mu=E)
#   4. FDR 校正（BH 方法）
#
# 产出：
#   - family_enrichment_ca.tsv:
#     family1 / family2 / observed / expected / fold_enrichment / pval / fdr
#   - family_enrichment_no.tsv: 同上
#
# 对应旧代码：step2_poissontest.py
# 改进：
#   - iterrows() → 向量化（pandas groupby + scipy vectorized）
#   - 分 ca/no 两组各做一次，而非混在一起
