#!/bin/bash



# Activate conda environment
source /home/caishangLab/fangxiunan/anaconda3/etc/profile.d/conda.sh
conda activate metagenomics




## Define paths
#SUMMARY_FILE="/storage/caishangLab/fangxiunan/ASAP/test260121-clustercore/ca1_foldseek/ca_skeleton_summary/proteinfamily_redefined_summary.csv"
#GO_FILE="/storage/caishangLab/fangxiunan/ASAP/test1226_go/mergedids_GOstats.csv"
#OUTPUT_DIR="/storage/caishangLab/fangxiunan/ASAP/test260121-clustercore/ca1_foldseek/ca_skeleton_summary"
#SCRIPT="/storage/caishangLab/fangxiunan/ASAP/test260121-clustercore/plot_top10_organism.py"
#
#
#
#
## Create output directory
#mkdir -p "$OUTPUT_DIR"
#
#
#
#
## Run Python script
#python "$SCRIPT" \
#    --summary-file "$SUMMARY_FILE" \
#    --go-file "$GO_FILE" \
#    --output-dir "$OUTPUT_DIR"
#
#
#
#
#echo "Top 10 organism plot generated and saved to $OUTPUT_DIR"
#
#
##########################################
## Define paths
#SUMMARY_FILE="/storage/caishangLab/fangxiunan/ASAP/test260121-clustercore/no1_foldseek/no_skeleton_summary/proteinfamily_redefined_summary.csv"
#GO_FILE="/storage/caishangLab/fangxiunan/ASAP/test1226_go/mergedids_GOstats.csv"
#OUTPUT_DIR="/storage/caishangLab/fangxiunan/ASAP/test260121-clustercore/no1_foldseek/no_skeleton_summary"
#SCRIPT="/storage/caishangLab/fangxiunan/ASAP/test260121-clustercore/plot_top10_organism.py"
#
#
#
#
## Create output directory
#mkdir -p "$OUTPUT_DIR"
#
#
#
#
## Run Python script
#python "$SCRIPT" \
#    --summary-file "$SUMMARY_FILE" \
#    --go-file "$GO_FILE" \
#    --output-dir "$OUTPUT_DIR"
#
#
#
#
#echo "Top 10 organism plot generated and saved to $OUTPUT_DIR"


# Define paths
CA_SUMMARY="/storage/caishangLab/fangxiunan/ASAP/test260121-clustercore/ca1_foldseek/ca_skeleton_summary/proteinfamily_redefined_summary.csv"
NO_SUMMARY="/storage/caishangLab/fangxiunan/ASAP/test260121-clustercore/no1_foldseek/no_skeleton_summary/proteinfamily_redefined_summary.csv"
GO_FILE="/storage/caishangLab/fangxiunan/ASAP/test1226_go/mergedids_GOstats.csv"
OUTPUT_DIR="/storage/caishangLab/fangxiunan/ASAP/test260121-clustercore/combined_plots"
SCRIPT="/storage/caishangLab/fangxiunan/ASAP/test260121-clustercore/plot_top20_sidebyside.py"




# Create output directory
mkdir -p "$OUTPUT_DIR"




# Run Python script
python "$SCRIPT" \
    --ca-summary "$CA_SUMMARY" \
    --no-summary "$NO_SUMMARY" \
    --go-file "$GO_FILE" \
    --output-dir "$OUTPUT_DIR"




echo "Combined CA vs NO top 20 organism plot generated and saved to $OUTPUT_DIR"