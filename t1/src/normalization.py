"""
Normalization methods for metagenomic protein count data.

Metagenomic data is compositional: different samples have different
sequencing depths (library sizes). Normalization removes this technical
variation so that biological differences can be detected.

Methods:
    relative_abundance  - proportions (sum to 1 per sample)
    cpm                 - counts per million
    log_cpm             - log2(CPM + 1)
    clr                 - centered log-ratio (compositional)
    none                - pass-through (raw counts)
"""

import pandas as pd
import numpy as np


def normalize(df: pd.DataFrame, method: str = 'relative_abundance') -> pd.DataFrame:
    """
    Normalize a samples x proteins count matrix.

    Args:
        df: Raw count DataFrame (samples as rows, proteins as columns).
        method: One of 'relative_abundance', 'cpm', 'log_cpm', 'clr', 'none'.

    Returns:
        Normalized DataFrame with the same shape and index/columns.
    """
    if method == 'none':
        return df

    total_counts = df.sum(axis=1)

    if method == 'relative_abundance':
        return df.div(total_counts, axis=0)

    elif method == 'cpm':
        return df.div(total_counts, axis=0) * 1e6

    elif method == 'log_cpm':
        cpm = df.div(total_counts, axis=0) * 1e6
        return np.log2(cpm + 1)

    elif method == 'clr':
        pseudo = df + 1
        geom_means = np.exp(np.log(pseudo).mean(axis=1))
        return np.log(pseudo.div(geom_means, axis=0))

    else:
        raise ValueError(
            f"Unknown method '{method}'. "
            f"Choose from: relative_abundance, cpm, log_cpm, clr, none"
        )


def library_size_stats(df: pd.DataFrame, label: str = "") -> dict:
    """
    Compute library size statistics for a count matrix.

    Returns dict with min, max, mean, median, cv.
    Prints a warning if CV > 30%.
    """
    lib_sizes = df.sum(axis=1)
    cv = lib_sizes.std() / lib_sizes.mean() if lib_sizes.mean() > 0 else 0

    stats = {
        'label': label,
        'min': lib_sizes.min(),
        'max': lib_sizes.max(),
        'mean': lib_sizes.mean(),
        'median': lib_sizes.median(),
        'cv': cv,
    }

    print(f"{label} library sizes: "
          f"min={stats['min']:.0f}, max={stats['max']:.0f}, "
          f"mean={stats['mean']:.0f}, CV={cv:.1%}")

    if cv > 0.3:
        print(f"  WARNING: High library size variation (CV={cv:.1%} > 30%). "
              f"Normalization is strongly recommended.")

    return stats
