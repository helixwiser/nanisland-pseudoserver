"""
Visualization functions for metagenomic protein analysis.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from typing import Optional, Tuple

sns.set_style("whitegrid")
plt.rcParams['figure.dpi'] = 100
plt.rcParams['savefig.dpi'] = 300


def volcano_plot(
    results_df: pd.DataFrame,
    pval_threshold: float = 0.05,
    log2fc_threshold: float = 1.0,
    use_adjusted_pval: bool = True,
    figsize: Tuple[int, int] = (12, 8),
    output_path: Optional[str] = None,
) -> Tuple[plt.Figure, plt.Axes]:
    """
    Volcano plot of log2 fold change vs -log10 p-value.

    Points are colored by regulation category:
        Red  = upregulated in cancer
        Blue = downregulated in cancer (upregulated in normal)
        Gold = significant but low fold change
        Grey = not significant
    """
    df = results_df.copy()

    if use_adjusted_pval:
        pval_col, neg_log_col = 'p_adj', 'neg_log10_padj'
        pval_label = 'FDR-adjusted p-value'
    else:
        pval_col, neg_log_col = 'p_value', 'neg_log10_pval'
        pval_label = 'p-value'

    df = df.dropna(subset=['log2_fold_change', neg_log_col])

    df['regulation'] = 'Not significant'
    sig = df[pval_col] < pval_threshold
    df.loc[sig & (df['log2_fold_change'] > log2fc_threshold), 'regulation'] = 'Upregulated'
    df.loc[sig & (df['log2_fold_change'] < -log2fc_threshold), 'regulation'] = 'Downregulated'
    df.loc[sig & (df['log2_fold_change'].abs() <= log2fc_threshold), 'regulation'] = 'Significant (low FC)'

    colors = {
        'Not significant': '#CCCCCC',
        'Significant (low FC)': '#FFD700',
        'Upregulated': '#FF4444',
        'Downregulated': '#4444FF',
    }

    fig, ax = plt.subplots(figsize=figsize)

    for cat in ['Not significant', 'Significant (low FC)', 'Downregulated', 'Upregulated']:
        subset = df[df['regulation'] == cat]
        if len(subset) == 0:
            continue
        ax.scatter(
            subset['log2_fold_change'], subset[neg_log_col],
            c=colors[cat], label=f"{cat} ({len(subset)})",
            alpha=0.6, s=20, edgecolors='none',
        )

    ax.axhline(-np.log10(pval_threshold), color='black', ls='--', lw=1, alpha=0.5,
               label=f'{pval_label} = {pval_threshold}')
    ax.axvline(log2fc_threshold, color='black', ls='--', lw=1, alpha=0.5)
    ax.axvline(-log2fc_threshold, color='black', ls='--', lw=1, alpha=0.5)

    ax.set_xlabel('Log2 Fold Change (Cancer / Normal)', fontsize=12, fontweight='bold')
    ax.set_ylabel(f'-Log10 ({pval_label})', fontsize=12, fontweight='bold')
    ax.set_title('Volcano Plot: Cancer vs Normal Skin\nPaired T-Test Protein Expression',
                 fontsize=14, fontweight='bold', pad=20)
    ax.legend(loc='upper right', frameon=True, fancybox=True, shadow=True)
    ax.grid(True, alpha=0.3, linestyle=':')
    plt.tight_layout()

    if output_path:
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        print(f"Saved volcano plot: {output_path}")

    return fig, ax
