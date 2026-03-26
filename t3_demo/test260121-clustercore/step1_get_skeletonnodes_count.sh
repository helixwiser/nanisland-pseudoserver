#!/bin/bash

#SBATCH -p intel-sc3,amd-ep2,intel-sc3-32c
#SBATCH -q normal
#SBATCH --cpus-per-task=16
#SBATCH --output=/storage/caishangLab/fangxiunan/log/exec_%A_%a.log
#SBATCH --error=/storage/caishangLab/fangxiunan/log/exec_%A_%a.err
#SBATCH --mem=50G

# Source conda and activate environment
source /home/caishangLab/fangxiunan/anaconda3/etc/profile.d/conda.sh
conda activate metagenomics
processdir="/storage/caishangLab/fangxiunan/ASAP/test260121-clustercore/ca1_foldseek"
rm -rf $processdir/ca_skeleton_membercount && mkdir -p $processdir/ca_skeleton_membercount
cd $processdir/ca_skeleton_membercount
# Define fixed inputs
SKELETON_EDGE=$processdir/fa_layouttest/skeleton_edges.csv
CLUSTER_FILE=$processdir/cluster_data/foldseek-res_cluster.tsv
INPUT_DIR="/storage/caishangLab/fangxiunan/project-ongoing/metagenomic/maf_processedresult"
SCRIPT="/storage/caishangLab/fangxiunan/ASAP/test260120-extractcore/get_skeletonnodes_count.py"


# Loop over all _uniprot.1.tsv and _uniprot.2.tsv files
for uniprot_file in "$INPUT_DIR"/*_uniprot.[1].tsv; do
    if [[ ! -f "$uniprot_file" ]]; then
        continue  # Skip if no files match
    fi


    # Extract filename base (e.g., H100-cancer-skinread_uniprot.1.tsv -> H100-cancer-skinread)
    base_name=$(basename "$uniprot_file" | sed -E 's/cancer-skinread_uniprot\.[1]\.tsv$//')


    # Define output filename
    output_file="${base_name}_ca_skeleton_membercount.csv"


    # Run Python script
    python "$SCRIPT" \
        --skeleton-edge "$SKELETON_EDGE" \
        --cluster "$CLUSTER_FILE" \
        --uniprot-read "$uniprot_file" \
        --output "$output_file"


    echo "Processed: $uniprot_file -> $output_file"
done


echo "All jobs completed."




