# build_family_map.py
# 构建蛋白→家族映射表，附带来源标签和家族大小
#
# 输入：
#   - foldseek-res_cluster.tsv: Foldseek 聚类结果（rep, member）
#   - protein_source_labels.tsv: step0 产出的来源标签表
#
# 逻辑：
#   1. 从聚类结果提取 member → rep 映射
#   2. Foldseek ID 转 UniProt ID（AF-Q7DFV4-F1-model_v6 → Q7DFV4）
#   3. 合并来源标签（ca_only / no_only / shared）
#   4. 统计每个家族的成员数（总数 / ca侧成员数 / no侧成员数）
#
# 产出：
#   - protein_to_family_labeled.tsv: 四列
#     protein / family_rep / source / family_size
#   - family_composition.tsv: 每个家族的成员组成
#     family_rep / total_members / ca_members / no_members / shared_members
#
# 下游使用：
#   - step3 建网络时按 source 分组计数
#   - step3 归一化边权时用 ca_members / no_members 做分母
