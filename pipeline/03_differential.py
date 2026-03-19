#!/usr/bin/env python3
"""
Step 3: Differential expression analysis.

Loads count matrices, normalizes, runs paired t-tests,
creates a volcano plot, and saves significant protein lists.
"""

import sys
import os
import yaml
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from src.normalization import normalize, library_size_stats
from src.statistics import get_paired_data, paired_ttest, filter_significant, save_protein_lists, save_results
from src.plotting import volcano_plot


def main(config_path: str = "config/params.yaml"):
    with open(config_path) as f:
        cfg = yaml.safe_load(f)

    paths = cfg['paths']
    params = cfg['analysis']

    read_number = params.get('read_number', 1)
    ca_path = os.path.join(paths['count_matrices_dir'],
                           f"cancer-skin_all_uniprotproID_results.{read_number}.csv")
    no_path = os.path.join(paths['count_matrices_dir'],
                           f"normal-skin_all_uniprotproID_results.{read_number}.csv")
    diff_dir = paths['differential_dir']
    fig_dir = paths['figures_dir']
    os.makedirs(diff_dir, exist_ok=True)
    os.makedirs(fig_dir, exist_ok=True)

    # Load raw counts
    print("Loading count matrices...")
    ca_raw = pd.read_csv(ca_path, index_col=0)
    no_raw = pd.read_csv(no_path, index_col=0)

    # Filter samples
    sample_filter = params.get('sample_filter', 'H')
    if sample_filter:
        ca_raw = ca_raw[ca_raw.index.str.startswith(sample_filter)]
        no_raw = no_raw[no_raw.index.str.startswith(sample_filter)]

    print(f"Cancer: {ca_raw.shape}  |  Normal: {no_raw.shape}")

    # Library size stats
    library_size_stats(ca_raw, "Cancer")
    library_size_stats(no_raw, "Normal")

    # Normalize
    norm_method = params.get('normalization', 'log_cpm')
    print(f"\nNormalization: {norm_method}")
    ca_norm = normalize(ca_raw, method=norm_method)
    no_norm = normalize(no_raw, method=norm_method)

    # Paired samples
    common_samples, common_proteins = get_paired_data(ca_norm, no_norm)
    ca_sub = ca_norm.loc[common_samples, common_proteins]
    no_sub = no_norm.loc[common_samples, common_proteins]
    ca_sub_raw = ca_raw.loc[common_samples, common_proteins]
    no_sub_raw = no_raw.loc[common_samples, common_proteins]

    # Paired t-test
    print("\nRunning paired t-tests...")
    results_df, log2fc_df = paired_ttest(ca_sub, no_sub)

    # Volcano plot
    pval_thresh = params.get('pval_threshold', 0.05)
    lfc_thresh = params.get('log2fc_threshold', 0.1)
    use_adj = params.get('use_adjusted_pval', True)

    volcano_path = os.path.join(fig_dir, "volcano_plot.png")
    volcano_plot(results_df,
                 pval_threshold=pval_thresh,
                 log2fc_threshold=lfc_thresh,
                 use_adjusted_pval=use_adj,
                 output_path=volcano_path)

    # Save all results
    print("\nSaving results...")
    save_results(ca_sub, no_sub, results_df, log2fc_df, diff_dir,
                 ca_raw=ca_sub_raw, no_raw=no_sub_raw)

    # Significant protein lists
    sig_pval = params.get('sig_pval_threshold', 0.1)
    sig_lfc = params.get('sig_log2fc_threshold', 0.1)
    sig_adj = params.get('sig_use_adjusted', False)

    cancer_p, normal_p, nonsig_p = filter_significant(
        results_df, pval_threshold=sig_pval,
        log2fc_threshold=sig_lfc, use_adjusted=sig_adj)

    save_protein_lists(cancer_p, normal_p, nonsig_p, diff_dir)

    print(f"\nDone. Results in: {diff_dir}")
    print(f"       Figures in: {fig_dir}")


if __name__ == '__main__':
    config = sys.argv[1] if len(sys.argv) > 1 else "config/params.yaml"
    main(config)
