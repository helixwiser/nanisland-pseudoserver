#!/bin/bash
#SBATCH -p intel-sc3
#SBATCH -q normal
#SBATCH --cpus-per-task=4
#SBATCH --output=/storage/caishangLab/fangxiunan/log/test_rep_%A_%a.log
#SBATCH --error=/storage/caishangLab/fangxiunan/log/test_rep_%A_%a.err
#SBATCH --mem=10G


DETAIL_FILE="/storage/caishangLab/fangxiunan/database/merged_details_unique.tsv"
SCRIPT="/storage/caishangLab/fangxiunan/ASAP/test260121-clustercore/redefine_representative.py"


################################
process_dir=/storage/caishangLab/fangxiunan/ASAP/test260121-clustercore/ca1_foldseek
INPUT_DIR=$process_dir/ca_skeleton_membercount
OUTPUT_DIR=$process_dir/ca_skeleton_summary

# 创建输出和临时目录
mkdir -p "$OUTPUT_DIR"
rm -rf "$OUTPUT_DIR" && mkdir -p "$OUTPUT_DIR"




# 激活环境并运行
source /home/caishangLab/fangxiunan/anaconda3/etc/profile.d/conda.sh
conda activate metagenomics


python redefine_representative.py \
  --input-pattern "$INPUT_DIR/*_membercount.csv" \
  --detail-file "$DETAIL_FILE" \
  --output-dir "$OUTPUT_DIR" \
  --tissue-type cancer





################################
process_dir=/storage/caishangLab/fangxiunan/ASAP/test260121-clustercore/no1_foldseek
INPUT_DIR=$process_dir/no_skeleton_membercount
OUTPUT_DIR=$process_dir/no_skeleton_summary

# 创建输出和临时目录
mkdir -p "$OUTPUT_DIR"
rm -rf "$OUTPUT_DIR" && mkdir -p "$OUTPUT_DIR"




# 激活环境并运行
source /home/caishangLab/fangxiunan/anaconda3/etc/profile.d/conda.sh
conda activate metagenomics


python redefine_representative.py \
  --input-pattern "$INPUT_DIR/*_membercount.csv" \
  --detail-file "$DETAIL_FILE" \
  --output-dir "$OUTPUT_DIR" \
  --tissue-type normal