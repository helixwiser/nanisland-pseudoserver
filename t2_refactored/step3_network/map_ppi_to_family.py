# map_ppi_to_family.py
# 蛋白级 PPI 升级为家族级，并标记每条 PPI 的来源
#
# 输入：
#   - string_uniprot_subnet_renamed.csv: UniProt 蛋白级 PPI 网络
#   - protein_to_family_labeled.tsv: 蛋白→家族映射（含 source 标签）
#
# 逻辑：
#   1. 读取蛋白级 PPI，两端蛋白都在映射表中的保留
#   2. protein1 → family1, protein2 → family2
#   3. 去掉自环（family1 == family2）
#   4. 为每条 PPI 标记来源，调用 utils/labels.py 的标记逻辑：
#      - shared × shared → baseline
#      - ca_only × ca_only → ca_contribution
#      - ca_only × shared → ca_contribution
#      - no_only × no_only → no_contribution
#      - no_only × shared → no_contribution
#      - ca_only × no_only → cross（特殊情况）
#
# 产出：
#   - ppi_in_families_labeled.csv:
#     protein1 / protein2 / family1 / family2 / source1 / source2 / ppi_label
#
# 对应旧代码：step2_mapppis_mmseq2.py + step2_mapppis_foldseek.py
# 改进：两个脚本合并；新增来源标记逻辑；Foldseek ID 提取由 build_family_map 处理
