"""
MAF file parsing and paired patient discovery.

Reads LAST aligner MAF output, extracts protein alignment records,
and identifies patients with both cancer-skin and normal-skin samples.
"""

from typing import Dict, List
import os
import re
from collections import defaultdict


def parse_maf_file(maf_file_path: str) -> List[str]:
    """
    Parse a single MAF file and extract protein alignment records.

    Each MAF block is 4 lines: header (a), query (s), subject (s), blank.
    Returns tab-separated records with score fields, query protein, query seq,
    subject coords, and subject seq.
    """
    lines = []
    try:
        with open(maf_file_path, 'r') as f_in:
            for line in f_in:
                if not line.startswith('#'):
                    lines.append(line.strip())
    except FileNotFoundError:
        print(f"Warning: File not found: {maf_file_path}")
        return []

    results = []
    for i in range(0, len(lines), 4):
        group = lines[i:i + 4]
        if len(group) < 4:
            continue

        line1, line2, line3, line4 = group

        if not (line1.startswith('a') and line2.startswith('s') and
                line3.startswith('s') and line4.strip() == ''):
            continue

        tokens1 = line1.split()
        tokens2 = line2.split()
        tokens3 = line3.split()

        if len(tokens2) >= 7 and len(tokens3) >= 7:
            val00 = tokens1[1]
            val01 = tokens1[2]
            val02 = tokens1[3]
            val1 = tokens2[1]
            val2 = tokens2[6]
            vals3 = '_'.join(tokens3[1:6])
            val4 = tokens3[6]
            results.append('\t'.join([val00, val01, val02, val1, val2, vals3, val4]))

    return results


def find_paired_patients(folder_path: str) -> List[str]:
    """
    Find patient IDs that have both cancer-skin and normal-skin MAF files.

    Scans folder for files matching: {patient_id}-{type}.{number}.out.sprot.maf
    Returns sorted list of patient IDs with paired samples.
    """
    if not os.path.exists(folder_path):
        raise FileNotFoundError(f"Folder {folder_path} does not exist")

    patients = defaultdict(lambda: {'cancer-skin': set(), 'normal-skin': set(), 'oral': set()})
    pattern = r'^([^-]+)-(cancer-skin|normal-skin|oral)\.(\d+)\.out\.sprot\.maf$'

    for filename in os.listdir(folder_path):
        match = re.match(pattern, filename)
        if match:
            patient_id = match.group(1)
            file_type = match.group(2)
            patients[patient_id][file_type].add(match.group(3))

    paired = [
        pid for pid, types in patients.items()
        if types['cancer-skin'] and types['normal-skin']
    ]
    return sorted(paired)


def process_maf_proteins(
    input_dir: str,
    output_dir: str,
    patient_ids: List[str],
    tissue_type: str = 'cancer-skin',
) -> Dict[str, List[str]]:
    """
    Parse MAF files for multiple patients and write per-patient TSVs.

    Processes both .1 and .2 read files for each patient.
    Returns dict mapping patient_id -> list of alignment records.
    """
    os.makedirs(output_dir, exist_ok=True)
    pid_pro_dict: Dict[str, List[str]] = {}

    for pid in patient_ids:
        pid_pro_dict[pid] = []

        for suffix in ['.1', '.2']:
            maf_file = os.path.join(input_dir, f"{pid}-{tissue_type}{suffix}.out.sprot.maf")
            results = parse_maf_file(maf_file)

            output_file = os.path.join(output_dir, f"{pid}-{tissue_type}read_uniprot{suffix}.tsv")
            with open(output_file, 'w') as f:
                for item in results:
                    f.write(item + '\n')

            pid_pro_dict[pid].extend(results)

    return pid_pro_dict
