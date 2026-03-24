# labels.py
# 蛋白来源标签逻辑（step0 和 step3 共用）
#
# 核心函数：
#
# assign_protein_source(protein_id, ca_set, no_set) → str
#   - 在 ca_set 且在 no_set → "shared"
#   - 仅在 ca_set → "ca_only"
#   - 仅在 no_set → "no_only"
#
# assign_ppi_label(source1, source2) → str
#   来源标记规则：
#   - shared   × shared   → "baseline"
#   - ca_only  × ca_only  → "ca_contribution"
#   - ca_only  × shared   → "ca_contribution"
#   - no_only  × no_only  → "no_contribution"
#   - no_only  × shared   → "no_contribution"
#   - ca_only  × no_only  → "cross"
#
# 注意：
#   - cross 类型的 PPI 在实际数据中应该较少（两端蛋白分别只在一组出现）
#   - cross 类型不计入 ca 或 no 网络，单独标记备查
