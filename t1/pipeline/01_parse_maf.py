#!/usr/bin/env python3
"""
Step 1: Parse MAF files and produce per-patient TSVs.

Reads raw MAF alignment files, extracts protein records,
and writes one TSV per patient/tissue/read.
"""

import sys
import os
import yaml

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from src.maf_parser import find_paired_patients, process_maf_proteins


def main(config_path: str = "config/params.yaml"):
    with open(config_path) as f:
        cfg = yaml.safe_load(f)

    raw_dir = cfg['paths']['raw_maf_dir']
    pre_dir = cfg['paths']['pre_res_dir']

    pid_list = find_paired_patients(raw_dir)
    print(f"Found {len(pid_list)} patients with paired cancer/normal samples\n")

    for tissue in ['cancer-skin', 'normal-skin']:
        print(f"Processing {tissue}...")
        process_maf_proteins(raw_dir, pre_dir, pid_list, tissue)

    print("\nDone. TSV files written to:", pre_dir)


if __name__ == '__main__':
    config = sys.argv[1] if len(sys.argv) > 1 else "config/params.yaml"
    main(config)
