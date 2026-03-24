"""
Statistical analysis for paired cancer vs normal protein expression.

Performs per-protein paired t-tests with FDR correction and
generates significance-filtered protein lists.
"""

import pandas as pd
import numpy as np
from scipy.stats import ttest_rel
from statsmodels.stats.multitest import multipletests
from typing import Tuple, List, Optional
import os


def get_paired_data(
    ca_df: pd.DataFrame,
    no_df: pd.DataFrame,
    min_median: Optional[float] = None,
) -> Tuple[pd.Index, List[str]]:
    """
    Identify paired samples and common proteins between two datasets.

    Args:
        ca_df: Cancer expression matrix (samples x proteins).
        no_df: Normal expression matrix (samples x proteins).
        min_median: If set, exclude samples with median expression below this.

    Returns:
        (common_samples, common_proteins)
    """
    if min_median is not None:
        valid_ca = ca_df[ca_df.median(axis=1) >= min_median].index
        valid_no = no_df[no_df.median(axis=1) >= min_median].index
    else:
        valid_ca = ca_df.index
        valid_no = no_df.index

    common_samples = valid_ca.intersection(valid_no)
    common_proteins = sorted(set(ca_df.columns) & set(no_df.columns))

    print(f"Paired samples: {len(common_samples)}  |  Common proteins: {len(common_proteins)}")
    return common_samples, common_proteins


def paired_ttest(
    ca_df: pd.DataFrame,
    no_df: pd.DataFrame,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Per-protein paired t-test between cancer and normal samples.

    Both DataFrames must share the same index (patients) and columns (proteins).

    Returns:
        results_df: One row per protein with stats, p-values, and significance flags.
        log2fc_df:  Per-sample log2 fold changes (patients x proteins).
    """
    proteins = ca_df.columns.tolist()
    results = []
    log2fc_dict = {}

    for i, protein in enumerate(proteins):
        if (i + 1) % 10000 == 0:
            print(f"  Testing protein {i+1}/{len(proteins)}...")

        ca_vals = ca_df[protein].values
        no_vals = no_df[protein].values

        ca_mean = np.mean(ca_vals)
        no_mean = np.mean(no_vals)
        ca_std = np.std(ca_vals, ddof=1)
        no_std = np.std(no_vals, ddof=1)

        # Per-sample log2FC with pseudocount; NaN where both are zero
        zero_mask = (ca_vals == 0) & (no_vals == 0)
        lfc = np.log2((ca_vals + 1) / (no_vals + 1))
        lfc[zero_mask] = np.nan
        log2fc_dict[protein] = lfc

        try:
            if ca_std > 0 or no_std > 0:
                t_stat, p_value = ttest_rel(ca_vals, no_vals)
            else:
                t_stat, p_value = np.nan, np.nan
        except Exception:
            t_stat, p_value = np.nan, np.nan

        results.append({
            'protein': protein,
            'cancer_mean': ca_mean,
            'normal_mean': no_mean,
            'cancer_std': ca_std,
            'normal_std': no_std,
            'mean_difference': ca_mean - no_mean,
            'log2_fold_change': np.nanmean(lfc),
            't_statistic': t_stat,
            'p_value': p_value,
        })

    results_df = pd.DataFrame(results)
    log2fc_df = pd.DataFrame(log2fc_dict, index=ca_df.index)

    # FDR correction
    _, pvals_adj, _, _ = multipletests(results_df['p_value'].fillna(1), method='fdr_bh')
    results_df['p_adj'] = pvals_adj

    # Significance flags
    for col, thresh in [('p01', 0.1), ('p005', 0.05), ('p001', 0.01)]:
        results_df[f'sig_{col}'] = results_df['p_value'] < thresh
        results_df[f'sig_adj_{col}'] = results_df['p_adj'] < thresh

    results_df['neg_log10_pval'] = -np.log10(results_df['p_value'])
    results_df['neg_log10_padj'] = -np.log10(results_df['p_adj'])

    print(f"Completed paired t-test for {len(results_df)} proteins")
    return results_df, log2fc_df


def filter_significant(
    results_df: pd.DataFrame,
    pval_threshold: float = 0.1,
    log2fc_threshold: float = 0.1,
    use_adjusted: bool = False,
) -> Tuple[List[str], List[str], List[str]]:
    """
    Split proteins into cancer-upregulated, normal-upregulated, and non-significant.

    Returns:
        (cancer_proteins, normal_proteins, nonsig_proteins) as lists of protein IDs.
    """
    pval_col = 'p_adj' if use_adjusted else 'p_value'

    sig_mask = results_df[pval_col] < pval_threshold
    up_mask = results_df['log2_fold_change'] > log2fc_threshold
    down_mask = results_df['log2_fold_change'] < -log2fc_threshold

    cancer_proteins = results_df.loc[sig_mask & up_mask, 'protein'].tolist()
    normal_proteins = results_df.loc[sig_mask & down_mask, 'protein'].tolist()
    nonsig_proteins = results_df.loc[~(sig_mask & (up_mask | down_mask)), 'protein'].tolist()

    print(f"Upregulated in cancer: {len(cancer_proteins)}")
    print(f"Upregulated in normal: {len(normal_proteins)}")
    print(f"Non-significant:       {len(nonsig_proteins)}")

    return cancer_proteins, normal_proteins, nonsig_proteins


def save_protein_lists(
    cancer_proteins: List[str],
    normal_proteins: List[str],
    nonsig_proteins: List[str],
    out_dir: str,
) -> None:
    """Write protein ID lists to text files."""
    os.makedirs(out_dir, exist_ok=True)

    for filename, proteins in [
        ("sig_uniprotproid_cancer.txt", cancer_proteins),
        ("sig_uniprotproid_normal.txt", normal_proteins),
        ("nonsig_uniprotproID_common.txt", nonsig_proteins),
    ]:
        path = os.path.join(out_dir, filename)
        with open(path, 'w') as f:
            f.write('\n'.join(proteins) + '\n' if proteins else '')
        print(f"  {filename}: {len(proteins)} proteins")


def save_results(
    ca_df: pd.DataFrame,
    no_df: pd.DataFrame,
    results_df: pd.DataFrame,
    log2fc_df: pd.DataFrame,
    out_dir: str,
    ca_raw: Optional[pd.DataFrame] = None,
    no_raw: Optional[pd.DataFrame] = None,
) -> None:
    """
    Save all analysis outputs to CSV files.

    If ca_raw/no_raw are provided, the count file uses raw counts;
    otherwise it uses the (possibly normalized) input data.
    """
    os.makedirs(out_dir, exist_ok=True)

    ca_labeled = (ca_raw if ca_raw is not None else ca_df).copy()
    no_labeled = (no_raw if no_raw is not None else no_df).copy()
    ca_labeled.index = ca_labeled.index + '-cancer'
    no_labeled.index = no_labeled.index + '-normal'

    df_count = pd.concat([ca_labeled, no_labeled], axis=0)

    df_sample = pd.DataFrame(df_count.index, columns=['sample_name'])
    df_sample['type'] = df_sample['sample_name'].str.split('-').str[1]
    df_sample.set_index('sample_name', inplace=True)

    files = {
        "uniprotproID_count_all.csv": df_count,
        "uniprotproID_sample_all.csv": df_sample,
        "paired_ttest_results.csv": results_df,
        "log2_fold_changes_per_sample.csv": log2fc_df,
    }

    for name, data in files.items():
        path = os.path.join(out_dir, name)
        data.to_csv(path, index=(name != "paired_ttest_results.csv"))
        print(f"  Saved: {name}")
