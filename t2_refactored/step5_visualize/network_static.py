# network_static.py
# 静态网络可视化（matplotlib + ForceAtlas2 / Kamada-Kawai 布局）
#
# 输入：
#   - network_ca.csv / network_no.csv: 加权网络
#   - topo_ca.csv / topo_no.csv: 拓扑指标
#   - communities_ca.tsv / communities_no.tsv: 社区划分
#
# 逻辑：
#   生成 2×2 四面板图（每组一张）：
#   1. 布局 + degree 着色
#   2. 布局 + 社区着色
#   3. hub 路径高亮
#   4. 核心子图 + UniProt 注释标签
#
# 产出：
#   - network_ca_four_panel.png（dpi=300）
#   - network_no_four_panel.png
#
# 对应旧代码：step3_test_fa.py（step12~step20）/ step3_test_kk.py
# 改进：
#   - 从 690 行单体脚本拆出
#   - 布局方法从 config 读取（ForceAtlas2 / Kamada-Kawai）
#   - 清理掉注释中残留的历史版本代码
