#!/usr/bin/env python3
"""
Step 4 (optional): Fetch FASTA sequences for significant proteins.

Reads protein ID lists from the differential results and downloads
sequences from UniProt.
"""

import sys
import os
import yaml

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from src.fasta import load_entries, fetch_fasta


def main(config_path: str = "config/params.yaml"):
    with open(config_path) as f:
        cfg = yaml.safe_load(f)

    diff_dir = cfg['paths']['differential_dir']
    fasta_dir = cfg['paths']['fasta_dir']
    os.makedirs(fasta_dir, exist_ok=True)

    # Fetch for cancer-significant proteins
    cancer_ids_file = os.path.join(diff_dir, "sig_uniprotproid_cancer.txt")
    if os.path.exists(cancer_ids_file):
        ids = load_entries(cancer_ids_file)
        if ids:
            print(f"\nFetching FASTA for {len(ids)} cancer-significant proteins...")
            fetch_fasta(ids,
                        out_fasta=os.path.join(fasta_dir, "sig_cancer_sequences.fasta"),
                        error_log=os.path.join(fasta_dir, "failed_cancer.txt"))

    # Fetch for normal-significant proteins
    normal_ids_file = os.path.join(diff_dir, "sig_uniprotproid_normal.txt")
    if os.path.exists(normal_ids_file):
        ids = load_entries(normal_ids_file)
        if ids:
            print(f"\nFetching FASTA for {len(ids)} normal-significant proteins...")
            fetch_fasta(ids,
                        out_fasta=os.path.join(fasta_dir, "sig_normal_sequences.fasta"),
                        error_log=os.path.join(fasta_dir, "failed_normal.txt"))

    print(f"\nDone. FASTA files in: {fasta_dir}")


if __name__ == '__main__':
    config = sys.argv[1] if len(sys.argv) > 1 else "config/params.yaml"
    main(config)
