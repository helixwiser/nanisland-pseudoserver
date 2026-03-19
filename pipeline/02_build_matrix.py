#!/usr/bin/env python3
"""
Step 2: Build protein count matrices from parsed TSVs.

Produces patients x proteins count CSVs for each tissue type and read number.
"""

import sys
import os
import yaml

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from src.maf_parser import find_paired_patients
from src.protein_matrix import build_count_matrix


def main(config_path: str = "config/params.yaml"):
    with open(config_path) as f:
        cfg = yaml.safe_load(f)

    raw_dir = cfg['paths']['raw_maf_dir']
    pre_dir = cfg['paths']['pre_res_dir']
    out_dir = cfg['paths']['count_matrices_dir']

    pid_list = find_paired_patients(raw_dir)
    print(f"Building count matrices for {len(pid_list)} patients\n")

    for tissue in ['cancer-skin', 'normal-skin']:
        for read_num in [1, 2]:
            print(f"  {tissue} read {read_num}...")
            build_count_matrix(pid_list, pre_dir, tissue,
                               read_number=read_num, out_dir=out_dir)

    print("\nDone. Count matrices written to:", out_dir)


if __name__ == '__main__':
    config = sys.argv[1] if len(sys.argv) > 1 else "config/params.yaml"
    main(config)
