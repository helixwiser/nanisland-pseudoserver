# merge_ids.py
# 合并 ca + no 的 UniProt ID，去冗余，生成来源标签表
#
# 输入：
#   - CA_IDS: cancer R1 的 UniProt ID 列表
#   - NO_IDS: normal R1 的 UniProt ID 列表
#
# 逻辑：
#   - 读取两组 ID
#   - 取并集去冗余
#   - 为每个蛋白打标签：ca_only / no_only / shared
#
# 产出：
#   - all_ids_merged.txt: 合并后的全量 UniProt ID（一列）
#   - protein_source_labels.tsv: 三列（protein_id / source / origin_file）
#     source 取值：ca_only, no_only, shared
