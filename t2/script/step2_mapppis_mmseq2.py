# Step 1: map_proteins_to_families.py
# 将每个蛋白映射到其所属的蛋白家族（以representative为key）
import pandas as pd
import sys


# 输入文件路径
cluster_file = "cluster_data/cluster_results.tsv"
output_map = "cluster_data/protein_to_family.tsv"


# 读取聚类结果：第一列是representative，第二列是成员
df = pd.read_csv(cluster_file, sep='\t', header=None, names=['rep', 'member'])

####### FOLDSEEK ID #####
## 提取 UniProt ID（去掉 AlphaFold 前缀和后缀，如 AF-Q7DFV4-F1-model_v6 → Q7DFV4）
#def extract_uniprot_id(s):
#    return s.split('-')[1]
#
#
#df['rep_uniprot'] = df['rep'].apply(extract_uniprot_id)
#df['member_uniprot'] = df['member'].apply(extract_uniprot_id)
####### FOLDSEEK ID #####

###### MMSEQ2 ID #####
df['rep_uniprot'] = df['rep']
df['member_uniprot'] = df['member']

###### MMSEQ2 ID #####
# 构建 member → rep 映射
mapping = df.set_index('member_uniprot')['rep_uniprot'].to_dict()


# 输出映射表
pd.DataFrame(list(mapping.items()), columns=['protein', 'family_rep']).to_csv(output_map, index=False)
print(f"[Step 1] Protein to family mapping saved to {output_map}")

ppi_file = "cluster_data/string_uniprot_subnet_renamed.csv"
mapping_file = "cluster_data/protein_to_family.tsv"
output_ppi_family = "cluster_data/ppi_in_families.csv"


# 读取映射
mapping = pd.read_csv(mapping_file, index_col='protein')['family_rep'].to_dict()


# 读取PPI
ppi = pd.read_csv(ppi_file)


# 仅保留两个蛋白都在映射中的PPI
ppi_filtered = ppi[
    (ppi['protein1'].isin(mapping)) & 
    (ppi['protein2'].isin(mapping))
].copy()


# 转换为家族级别互作
ppi_filtered['family1'] = ppi_filtered['protein1'].map(mapping)
ppi_filtered['family2'] = ppi_filtered['protein2'].map(mapping)


# 排除自环（可选）
ppi_filtered = ppi_filtered[ppi_filtered['family1'] != ppi_filtered['family2']]


# 保存家族间原始互作对
ppi_filtered.to_csv(output_ppi_family, index=False)
print(f"[Step 2] Filtered PPI at family level saved to {output_ppi_family}")