import pandas as pd
import gzip


# 读取已生成的子网络文件
subnet_df = pd.read_csv('string_uniprot_subnet.csv')


# 读取 UniProt_AC_ca.txt 获取目标 ID 集合
with open('UniProt_ids.txt', 'r') as f:
    uniprot_ac_set = set(line.strip() for line in f if line.strip())


# 读取 protein.aliases.v12.0.txt.gz，获取 STRING ID 到 UniProt AC 的映射
with gzip.open('/storage/caishangLab/fangxiunan/database/stringDB/protein.aliases.v12.0.txt.gz', 'rt') as f:
    aliases_df = pd.read_csv(f, sep='\t', header=0, names=['string_protein_id', 'alias', 'source'])


# 筛选 source 为 'UniProt_AC' 且 alias 在目标列表中的记录
uniprot_ac_map = aliases_df[
    (aliases_df['source'] == 'UniProt_AC') & 
    (aliases_df['alias'].isin(uniprot_ac_set))
][['string_protein_id', 'alias']].drop_duplicates()


# 提取子网络中所有出现的 protein1 和 protein2
network_ids = pd.concat([subnet_df['protein1'], subnet_df['protein2']]).unique()
valid_string_ids = set(network_ids)


# 找到这些 ID 对应的 UniProt AC
mapped_uniprot_in_network = uniprot_ac_map[uniprot_ac_map['string_protein_id'].isin(valid_string_ids)]


# 确保所有在子网络中映射到的 UniProt AC 确实来自输入列表
final_uniprot_ids = mapped_uniprot_in_network['alias'].unique()


# 统计每个 UniProt_AC 的边数（出现在子网络中的次数）
edge_count = []
for uid in final_uniprot_ids:
    string_id = uniprot_ac_map[uniprot_ac_map['alias'] == uid]['string_protein_id'].values
    if len(string_id) > 0:
        sid = string_id[0]
        count = ((subnet_df['protein1'] == sid) | (subnet_df['protein2'] == sid)).sum()
        edge_count.append({'id': uid, 'edge count': count})


# 转为 DataFrame 并保存
summary_df = pd.DataFrame(edge_count)
summary_df.to_csv('UniProt_AC_ca_networksummary.csv', index=False)