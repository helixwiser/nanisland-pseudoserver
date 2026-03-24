# extract_subnet.py
# 从 StringDB 提取 UniProt 子网络
#
# 输入：
#   - all_ids_merged.txt: 合并后的全量 UniProt ID
#   - STRINGDB_ALIASES: protein.aliases.v12.0.txt.gz
#   - STRINGDB_LINKS: protein.links.detailed.v12.0.txt.gz
#
# 逻辑：
#   1. 读取全量 UniProt ID
#   2. 在 aliases 中筛选 source=UniProt_AC 的记录，映射 UniProt → STRING ID
#   3. 分块读取 links（chunksize=1M），保留两端都在映射表中的边
#   4. STRING ID 回转 UniProt ID
#
# 产出：
#   - string_uniprot_alias_sub.csv: UniProt ↔ STRING 映射子集
#   - string_uniprot_subnet.csv: STRING ID 子网络
#   - string_uniprot_subnet_renamed.csv: UniProt ID 子网络（下游主要使用）
#
# 注意：
#   - 用全量蛋白查一次，不分 ca/no。后续通过标签拆分。
#   - 对应旧代码：step1_getstringdbsubnet.py（逻辑不变，路径从 config 读取）
