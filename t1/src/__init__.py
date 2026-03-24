"""
Aisland: Metagenomic protein expression analysis pipeline.

Modules:
    maf_parser      - Parse MAF alignment files, find paired patients
    protein_matrix  - Build protein count matrices from parsed data
    normalization   - Normalize count data (relative abundance, CPM, CLR, etc.)
    statistics      - Paired t-tests, FDR correction, significance filtering
    plotting        - Volcano plots, visualization helpers
    fasta           - Fetch FASTA sequences from UniProt
"""
