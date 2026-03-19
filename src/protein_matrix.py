"""
Build protein count matrices from parsed MAF/TSV data.

Reads per-patient TSV files, deduplicates PCR artifacts,
extracts UniProt IDs, and builds a patients x proteins count matrix.
"""

import pandas as pd
import os
from collections import Counter
from typing import List, Optional


def build_count_matrix(
    pid_list: List[str],
    pre_dir: str,
    tissue_type: str,
    read_number: int = 1,
    out_dir: Optional[str] = None,
) -> pd.DataFrame:
    """
    Build a protein count matrix from per-patient TSV files.

    For each patient, reads the parsed TSV, removes PCR duplicates (by
    sequence column 6), and extracts UniProt accession IDs from the
    sp|ACCESSION|ENTRY_SPECIES format in column 3.

    Args:
        pid_list: Patient IDs to include.
        pre_dir: Directory containing *read_uniprot.{read_number}.tsv files.
        tissue_type: 'cancer-skin' or 'normal-skin'.
        read_number: 1 or 2.
        out_dir: If provided, saves CSV to this directory.

    Returns:
        DataFrame with patients as rows and UniProt IDs as columns (integer counts).
    """
    if read_number not in [1, 2]:
        raise ValueError("read_number must be 1 or 2")

    pid_proteins: dict[str, list[str]] = {}

    for pid in pid_list:
        file_path = os.path.join(pre_dir, f"{pid}-{tissue_type}read_uniprot.{read_number}.tsv")
        tsv = pd.read_csv(file_path, sep='\t', header=None)

        # Deduplicate PCR artifacts (same sequence in column 6)
        tsv = tsv.drop_duplicates(subset=6, keep='first')

        # Extract UniProt accession: "sp|Q91YR9|PTGR1_MOUSE" -> "Q91YR9"
        names = tsv[3].tolist()
        pid_proteins[pid] = [name.split('|')[1] for name in names]

    # Pivot to count matrix
    rows = []
    for pid, proteins in pid_proteins.items():
        counts = Counter(proteins)
        counts['__index__'] = pid
        rows.append(counts)

    df = pd.DataFrame(rows).set_index('__index__').fillna(0).astype(int)
    df = df.reindex(sorted(df.columns), axis=1)

    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
        out_file = os.path.join(out_dir, f"{tissue_type}_all_uniprotproID_results.{read_number}.csv")
        df.to_csv(out_file)
        print(f"Saved: {out_file}  ({df.shape[0]} patients x {df.shape[1]} proteins)")

    return df
