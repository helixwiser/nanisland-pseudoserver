# Step 3: count_family_interactions.py
# 统计每对家族之间的互作数量
import pandas as pd
from itertools import combinations


input_file = "cluster_data/ppi_in_families.csv"
output_count = "cluster_data/family_interaction_count.tsv"


# 读取家族间互作
df = pd.read_csv(input_file)


# 标准化家族对顺序（保证 (A,B) 和 (B,A) 被视为同一对）
df['family_pair'] = df.apply(lambda x: tuple(sorted([x['family1'], x['family2']])), axis=1)


# 计数
count_series = df['family_pair'].value_counts()


# 转换为 DataFrame
count_df = pd.DataFrame(count_series.items(), columns=['family_pair', 'observed_interactions'])
count_df[['family1', 'family2']] = pd.DataFrame(count_df['family_pair'].tolist(), index=count_df.index)
count_df = count_df.drop('family_pair', axis=1)[['family1', 'family2', 'observed_interactions']]


count_df.to_csv(output_count, sep='\t', index=False)
print(f"[Step 3] Family interaction counts saved to {output_count}")
