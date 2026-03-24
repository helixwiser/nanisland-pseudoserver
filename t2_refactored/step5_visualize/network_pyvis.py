# network_pyvis.py
# 交互式网络可视化（PyVis）
#
# 输入：
#   - network_ca.csv / network_no.csv: 加权网络
#   - topo_ca.csv / topo_no.csv: 拓扑指标（用于节点大小/颜色）
#   - communities_ca.tsv / communities_no.tsv: 社区划分（用于着色，可选）
#
# 逻辑：
#   - 节点大小 = degree 归一化
#   - 节点颜色 = 社区 ID（或 degree 渐变）
#   - 边宽度 = norm_weight 归一化
#   - hover 显示：family_rep / degree / betweenness / pagerank
#   - 物理引擎：ForceAtlas2Based
#   - top20 节点显示 label
#
# 产出：
#   - network_ca_interactive.html
#   - network_no_interactive.html
#
# 对应旧代码：step2_networkpyvis.py + step3_networkbone_pyvis.py
# 改进：合并为一个脚本，ca/no 各生成一份
