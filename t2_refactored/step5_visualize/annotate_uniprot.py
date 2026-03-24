# annotate_uniprot.py
# UniProt 批量 ID mapping 注释
#
# 输入：
#   - hubs_ca.csv / hubs_no.csv: hub 蛋白家族列表（或全量家族 rep ID）
#
# 逻辑：
#   - 使用 UniProt 批量 ID mapping API（非逐条 REST 查询）
#   - 提取：protein_name / gene_name / organism / function description
#   - 合并到 hub 列表或拓扑表中
#
# 产出：
#   - uniprot_annotations.tsv: uniprot_id / protein_name / gene_name / organism
#   - hubs_ca_annotated.csv: hub 列表 + 注释列
#   - hubs_no_annotated.csv: 同上
#
# 对应旧代码：step3_tempgenerate_ppi_analysis_report.py
# 改进：
#   - 旧版逐条调 UniProt REST API（requests.get 循环），几千个蛋白极慢
#   - 改用 UniProt 批量 ID mapping API，一次请求搞定
