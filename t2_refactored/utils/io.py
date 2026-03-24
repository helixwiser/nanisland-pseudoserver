# io.py
# 统一读写函数
#
# 功能：
#   - read_id_list(path) → set: 读取单列 ID 文件
#   - read_tsv(path) → pd.DataFrame: 读取 TSV
#   - read_csv(path) → pd.DataFrame: 读取 CSV
#   - write_tsv(df, path): 写入 TSV
#   - write_csv(df, path): 写入 CSV
#   - read_stringdb_gz(path, chunksize) → iterator: 分块读取 StringDB 压缩文件
#
# 所有路径从 config.py 读取，这里只处理读写逻辑
# 统一 encoding='utf-8'，统一 lineterminator='\n'
