# count_interactions.py
# 按家族对 + 来源分组统计互作条数
#
# 输入：
#   - ppi_in_families_labeled.csv: 家族级 PPI（含 ppi_label）
#
# 逻辑：
#   1. 家族对排序标准化：sorted([family1, family2])，保证 (A,B)=(B,A)
#   2. 按 (family1, family2, ppi_label) 分组计数
#   3. pivot 成宽表，每对家族一行
#
# 产出：
#   - family_interaction_count_labeled.tsv:
#     family1 / family2 / count_baseline / count_ca / count_no / count_cross /
#     count_total_ca(=baseline+ca) / count_total_no(=baseline+no)
#
# 对应旧代码：step2_countppis.py（逻辑扩展，新增分来源计数）
