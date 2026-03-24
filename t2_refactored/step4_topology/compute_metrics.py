# compute_metrics.py
# 计算网络拓扑指标
#
# 输入：
#   - network_ca.csv / network_no.csv: 加权网络
#
# 逻辑：
#   对 ca / no 各做一次：
#   1. 构建 NetworkX 加权无向图
#   2. 计算三个拓扑指标：
#      - degree: 连接数（hub 识别）
#      - betweenness centrality: 桥接性（瓶颈节点）
#      - PageRank: 影响力（连接到重要节点的节点）
#   3. 不计算 closeness（与 betweenness 高度重叠，大网络计算慢）
#
# 产出：
#   - topo_ca.csv: family_rep / degree / betweenness / pagerank
#   - topo_no.csv: 同上
#
# 对应旧代码：step3_test_fa.py 中 step5 + step14 部分
# 改进：从 690 行单体脚本中拆出，独立运行
