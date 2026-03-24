# check_overlap.py
# 统计 ca 和 no 两组 UniProt ID 的重叠度
#
# 输入：
#   - CA_IDS: cancer R1 的 UniProt ID 列表
#   - NO_IDS: normal R1 的 UniProt ID 列表
#
# 逻辑：
#   - 计算 ca_only / no_only / shared 各多少个
#   - 计算重叠比例（shared / union）
#   - 打印摘要
#
# 产出：
#   - 终端输出：各集合大小、重叠比例
#   - overlap_stats.tsv（可选）：留档备查
#
# 对应旧代码：steppre_check_ids.py（重写，旧版 print 信息与变量名不一致）
