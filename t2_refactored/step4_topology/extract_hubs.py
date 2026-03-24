# extract_hubs.py
# 提取各组 hub 蛋白家族列表
#
# 输入：
#   - topo_ca.csv / topo_no.csv: 拓扑指标
#   - family_composition.tsv: 家族成员组成（可选，用于注释）
#
# 逻辑：
#   - 按 degree / betweenness / PageRank 各取 top N（N 从 config 读取）
#   - 取三个指标的并集作为候选 hub 列表
#   - 可选：合并 UniProt 功能注释
#
# 产出：
#   - hubs_ca.csv: family_rep / degree / betweenness / pagerank / rank_degree / rank_betweenness / rank_pagerank
#   - hubs_no.csv: 同上
#
# 下游使用：假说1.3 关键蛋白→关键菌映射的输入
