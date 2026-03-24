# build_network.py
# 构建 ca / no 两张加权网络
#
# 输入：
#   - family_enrichment_ca.tsv / family_enrichment_no.tsv: Poisson 检验结果
#   - family_interaction_count_labeled.tsv: 分来源计数
#   - family_composition.tsv: 各家族在 ca/no 侧的成员数
#
# 逻辑：
#   1. 筛选 FDR < 0.05 的显著边
#   2. 计算归一化边权：
#      norm_weight_ca(A,B) = count_total_ca / (|A_ca_members| × |B_ca_members|)
#      norm_weight_no 同理
#   3. 输出两张网络
#
# 产出：
#   - network_ca.csv: family1 / family2 / raw_count / norm_weight / fdr
#   - network_no.csv: 同上
#
# 对应旧代码：step2_buildweightednetwork.py
# 改进：
#   - 旧版只输出一张网络，边权为 raw observed count
#   - 新版分 ca/no 两张，边权归一化，消除家族大小偏差
