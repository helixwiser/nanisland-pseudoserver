# Step 5: build_weighted_network.py
# 构建加权家族网络（仅保留显著互作）
import pandas as pd


enrich_file = "cluster_data/family_enrichment_poisson.tsv"
output_network = "representation_protein_network.csv"


# 读取富集结果
df = pd.read_csv(enrich_file, sep='\t')


# 筛选显著互作（FDR < 0.05）
significant = df[df['fdr'] < 0.05].copy()


# 构建加权网络：family1, family2, weight = observed
significant['weight'] = significant['observed']


# 输出为边列表
network = significant[['family1', 'family2', 'weight']]
network.to_csv(output_network, index=False)


print(f"[Step 5] Weighted family network saved to {output_network}")
print(f"Total significant edges: {len(network)}")
