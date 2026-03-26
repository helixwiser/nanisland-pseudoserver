#!/bin/bash


#SBATCH -p intel-sc3
#SBATCH -q normal
#SBATCH --cpus-per-task=4
#SBATCH --output=/storage/caishangLab/fangxiunan/log/annotate_hubs_%A_%a.log
#SBATCH --error=/storage/caishangLab/fangxiunan/log/annotate_hubs_%A_%a.err
#SBATCH --mem=10G


# Activate conda environment
source /home/caishangLab/fangxiunan/anaconda3/etc/profile.d/conda.sh
conda activate metagenomics
SCRIPT="/storage/caishangLab/fangxiunan/ASAP/test260121-clustercore/annotate_hubnodes_redefined.py"
###################################
# Define paths
HUB_FILE="/storage/caishangLab/fangxiunan/ASAP/test260121-clustercore/ca1_foldseek/fa_layouttest/ppi_hub_nodes_annotated.csv"
SUMMARY_FILE="/storage/caishangLab/fangxiunan/ASAP/test260121-clustercore/ca1_foldseek/ca_skeleton_summary/proteinfamily_redefined_summary.csv"
DETAILS_FILE="/storage/caishangLab/fangxiunan/database/merged_details_unique.tsv"


OUTPUT_FILE="/storage/caishangLab/fangxiunan/ASAP/test260121-clustercore/ca1_foldseek/ca_skeleton_summary/cahubnodes_final_annotated.csv"



# Run Python script
python "$SCRIPT" \
    --hub-file "$HUB_FILE" \
    --summary-file "$SUMMARY_FILE" \
    --details-file "$DETAILS_FILE" \
    --output-file "$OUTPUT_FILE"


echo "Hub node annotation completed and saved to $OUTPUT_FILE"
###################################
# Define paths
HUB_FILE="/storage/caishangLab/fangxiunan/ASAP/test260121-clustercore/no1_foldseek/fa_layouttest/ppi_hub_nodes_annotated.csv"
SUMMARY_FILE="/storage/caishangLab/fangxiunan/ASAP/test260121-clustercore/no1_foldseek/no_skeleton_summary/proteinfamily_redefined_summary.csv"
DETAILS_FILE="/storage/caishangLab/fangxiunan/database/merged_details_unique.tsv"


OUTPUT_FILE="/storage/caishangLab/fangxiunan/ASAP/test260121-clustercore/no1_foldseek/no_skeleton_summary/nohubnodes_final_annotated.csv"



# Run Python script
python "$SCRIPT" \
    --hub-file "$HUB_FILE" \
    --summary-file "$SUMMARY_FILE" \
    --details-file "$DETAILS_FILE" \
    --output-file "$OUTPUT_FILE"


echo "Hub node annotation completed and saved to $OUTPUT_FILE"