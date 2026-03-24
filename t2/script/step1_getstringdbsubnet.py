# -*- coding: utf-8 -*-
import pandas as pd


UNIPROT_AC_PATH = 'UniProt_ids.txt'  # 例如: '/path/to/UniProt_AC_ca.txt'
ALIASES_GZ_PATH = '/storage/caishangLab/fangxiunan/database/stringDB/protein.aliases.v12.0.txt.gz'  # 例如: '/data/protein.aliases.v12.0.txt.gz'
LINKS_GZ_PATH = '/storage/caishangLab/fangxiunan/database/stringDB/protein.links.detailed.v12.0.txt.gz'  # 例如: '/data/protein.links.detailed.v12.0.txt.gz'



# Output file names
OUTPUT_SUBNET = 'string_uniprot_subnet.csv'
OUTPUT_ALIASES_SUB = 'string_uniprot_alias_sub.csv'
OUTPUT_RENAMED = 'string_uniprot_subnet_renamed.csv'




# Read UniProt_AC list
try:
    uniprot_ac = pd.read_csv(UNIPROT_AC_PATH, header=None, squeeze=True)
except TypeError:
    df = pd.read_csv(UNIPROT_AC_PATH, header=None)
    uniprot_ac = df.iloc[:, 0] if df.shape[1] == 1 else None
    if uniprot_ac is None:
        raise ValueError("Input file must be a single-column text file.")
uniprot_ac_clean = uniprot_ac.str.strip().str.replace(r'[,\s]', '', regex=True)
uniprot_ac_set = set(uniprot_ac_clean.drop_duplicates())




# Load aliases and filter for UniProt_AC entries
aliases = pd.read_csv(
    ALIASES_GZ_PATH,
    sep='\t',
    usecols=['#string_protein_id', 'alias', 'source'],
    dtype=str,
    compression='gzip',
    engine='c'
)
ac_filter = (aliases['source'] == 'UniProt_AC') & aliases['alias'].isin(uniprot_ac_set)
uniprot_map = aliases[ac_filter].drop_duplicates(subset=['#string_protein_id', 'alias'])
string_to_uniprot = uniprot_map.set_index('#string_protein_id')['alias'].to_dict()
valid_string_ids = set(string_to_uniprot.keys())




# Stream process links file
chunksize = 1_000_000
header = pd.read_csv(LINKS_GZ_PATH, sep=r'\s+', nrows=1, engine='c').columns.tolist()




# Write header
with open(OUTPUT_SUBNET, 'w') as f:
    f.write(','.join(header) + '\n')




# Process in chunks
for chunk in pd.read_csv(
    LINKS_GZ_PATH,
    sep=r'\s+',
    compression='gzip',
    dtype=str,
    chunksize=chunksize,
    skiprows=1,
    header=None,
    names=header
):
    mask = chunk['protein1'].isin(valid_string_ids) & chunk['protein2'].isin(valid_string_ids)
    if mask.any():
        chunk[mask].to_csv(OUTPUT_SUBNET, mode='a', header=False, index=False, lineterminator='\n')




# Save matched aliases
uniprot_map.to_csv(OUTPUT_ALIASES_SUB, index=False)




# Replace STRING IDs with UniProt ACs in network
subnet = pd.read_csv(OUTPUT_SUBNET)
subnet['protein1'] = subnet['protein1'].map(string_to_uniprot)
subnet['protein2'] = subnet['protein2'].map(string_to_uniprot)
subnet.dropna(subset=['protein1', 'protein2'], inplace=True)
subnet.to_csv(OUTPUT_RENAMED, index=False)