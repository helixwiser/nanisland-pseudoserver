# config.py
# 全局路径和参数统一管理
# 所有脚本从这里读取路径和参数，不再硬编码

# === 服务器路径 ===
# BASE_DIR: 项目根目录
# DB_DIR: 数据库目录（StringDB, AlphaFold DB, UniProt 等）
# INPUT_DIR: 上游 pipeline 产出目录（MAF 文件、UniProt ID 列表）
# OUTPUT_DIR: 本 pipeline 输出根目录

# === 输入文件 ===
# CA_IDS: cancer R1 的 UniProt ID 列表（来自 LAST 比对）
# NO_IDS: normal R1 的 UniProt ID 列表（来自 LAST 比对）

# === 数据库路径 ===
# STRINGDB_ALIASES: protein.aliases.v12.0.txt.gz
# STRINGDB_LINKS: protein.links.detailed.v12.0.txt.gz
# ALPHAFOLD_PDB: swissprot_pdb_v6.test.tar

# === Foldseek 参数 ===
# FOLDSEEK_COVERAGE = 0.5

# === 网络构建参数 ===
# FDR_THRESHOLD = 0.05
# POISSON_CHUNK_SIZE = 1_000_000  # StringDB links 分块读取大小

# === 拓扑分析参数 ===
# TOP_N_HUBS = 100  # 提取 hub 数量
# COMMUNITY_METHOD = "leiden"  # 待定：leiden / louvain
# COMMUNITY_RESOLUTION = 0.1  # 待调
