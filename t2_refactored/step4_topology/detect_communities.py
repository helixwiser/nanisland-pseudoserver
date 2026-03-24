# detect_communities.py
# 社区检测（方法待定）
#
# 输入：
#   - network_ca.csv / network_no.csv: 加权网络
#
# 逻辑：
#   - 待定：Leiden vs Louvain，有数据后对比选择
#   - resolution 参数待调
#   - 对 ca / no 各做一次
#
# 产出：
#   - communities_ca.tsv: family_rep / community_id
#   - communities_no.tsv: 同上
#   - community_stats.tsv: 各社区大小、模块度
#
# 对应旧代码：step3_test_fa.py 中 step11 部分
# 注意：旧代码试了 Louvain / Leiden / Label Propagation / Louvain多层四种，
#       均未收敛到最终选择。此处先预留接口，有数据后再定。
