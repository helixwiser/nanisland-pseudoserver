# Aisland

Metagenomic protein expression analysis pipeline — comparing cancer skin vs normal skin.

## Quick start

```bash
# 1. Edit config
vi config/params.yaml

# 2. Run full pipeline
bash pipeline/run_all.sh

# 3. Or run steps individually
python pipeline/01_parse_maf.py
python pipeline/02_build_matrix.py
python pipeline/03_differential.py
python pipeline/04_fetch_fasta.py        # optional
```

## Project structure

```
Aisland/
├── config/
│   └── params.yaml            # All paths + analysis parameters
├── src/                       # Shared library (import from here)
│   ├── maf_parser.py          # MAF file parsing, paired patient discovery
│   ├── protein_matrix.py      # Build count matrices from parsed TSVs
│   ├── normalization.py       # Normalize counts (rel. abundance, CPM, CLR, etc.)
│   ├── statistics.py          # Paired t-tests, FDR, significance filtering
│   ├── plotting.py            # Volcano plots
│   └── fasta.py               # Fetch FASTA from UniProt API
├── pipeline/                  # Runnable scripts (thin wrappers around src/)
│   ├── 01_parse_maf.py        # MAF → per-patient TSVs
│   ├── 02_build_matrix.py     # TSVs → count matrices
│   ├── 03_differential.py     # Count matrices → t-test + volcano + protein lists
│   ├── 04_fetch_fasta.py      # Protein IDs → FASTA sequences
│   └── run_all.sh             # Run everything end-to-end
├── notebooks/                 # Exploration & visualization (import from src/)
├── R/                         # R scripts (DESeq2, enrichment)
├── results/                   # Generated outputs (gitignored)
│   ├── pre_res/               # Parsed TSVs
│   ├── count_matrices/        # Protein count CSVs
│   ├── differential/          # T-test results, significant protein lists
│   ├── figures/               # Plots
│   └── fasta/                 # Downloaded sequences
└── requirements.txt
```

## Pipeline overview

```
raw MAF files
     │
     ▼  01_parse_maf.py
per-patient TSVs (results/pre_res/)
     │
     ▼  02_build_matrix.py
count matrices (results/count_matrices/)
     │
     ▼  03_differential.py
├── paired_ttest_results.csv
├── volcano_plot.png
├── sig_uniprotproid_cancer.txt
├── sig_uniprotproid_normal.txt
└── nonsig_uniprotproID_common.txt
     │
     ▼  04_fetch_fasta.py (optional)
FASTA sequences (results/fasta/)
```

## Configuration

All paths and parameters live in `config/params.yaml`. Key settings:

| Parameter | Default | Description |
|-----------|---------|-------------|
| `normalization` | `log_cpm` | `relative_abundance`, `cpm`, `log_cpm`, `clr`, `none` |
| `sample_filter` | `H` | Only include patient IDs starting with this |
| `pval_threshold` | `0.05` | Volcano plot p-value cutoff |
| `log2fc_threshold` | `0.1` | Volcano plot fold-change cutoff |

## Dependencies

```
pandas
numpy
scipy
statsmodels
matplotlib
seaborn
pyyaml
httpx          # for FASTA fetching (step 04 only)
```

## Using in notebooks

```python
import sys
sys.path.insert(0, '/mnt/resources/FXN/metagenomic/Aisland')

from src.normalization import normalize
from src.statistics import paired_ttest
from src.plotting import volcano_plot
```
