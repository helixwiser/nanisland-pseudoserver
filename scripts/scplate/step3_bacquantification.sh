#!/bin/bash
source /home/caishangLab/fangxiunan/anaconda3/etc/profile.d/conda.sh
conda activate metagenomics

usage() {
    echo "Usage: $0 -d destination_dir "
    exit 1

}

# Parse command-line arguments

while getopts ":d:" opt; do

    case $opt in

        d) destination_dir="$OPTARG"
        ;;
        *) usage

        ;;
    esac

done


process_sample() {
    local file=$1
    local destination_dir=$2
    local filename="${file##*/}"
    local folder_name="${filename%%_unmapped.dedupe.1.fq}"
    echo "Processing with ID: $folder_name"

    
    
    print_duration() {
        local start_seconds=$1
        local end_seconds=$2

        local duration=$((end_seconds - start_seconds))
        local duration_formatted=$(date -u -d @"$duration" +"%H:%M:%S")
        echo "Duration: $duration_formatted"
    }


    cd $destination_dir
    exec 3>> exec_time.log  
    mkdir -p result02_backraken2
    
    
    local unmappedread1=${folder_name}_unmapped.dedupe.1.fq
    local unmappedread2=${folder_name}_unmapped.dedupe.2.fq
    
    
    echo "Processing with ID: $folder_name" >&3
    ######################################
    ####### kraken2 quantification########
    ######################################
    #use confidence = 0.7
    
    #kraken2-inspect --db /storage/caishangLab/caishang/reference/kraken2-ref/ --threads 16 > inspect_report.txt
    local start_seconds=$(date +"%s")
    kraken2 \
        --db /storage/caishangLab/caishang/reference/kraken2-ref/ \
        --threads 32 \
        --confidence 0.7 \
        --report result02_backraken2/$folder_name.k2report \
        --report-minimizer-data \
        --minimum-hit-groups 3 \
        --paired $unmappedread1 $unmappedread2\
        > result02_backraken2/$folder_name.kraken2
    local end_seconds=$(date +"%s")
    
    echo "kraken2 finished!!!" >&3
    print_duration "$start_seconds" "$end_seconds" >&3
    
    local start_seconds=$(date +"%s")
    bracken \
        -d /storage/caishangLab/caishang/reference/kraken2-ref/ \
        -i result02_backraken2/$folder_name.k2report \
        -r 100 \
        -l S \
        -t 32 \
        -o result02_backraken2/$folder_name.bracken \
        -w result02_backraken2/$folder_name.brepo
    local end_seconds=$(date +"%s")
    
    echo "braken finished!!!" >&3
    print_duration "$start_seconds" "$end_seconds" >&3

    local start_seconds=$(date +"%s")
    python "/storage/caishangLab/fangxiunan/tools/KrakenTools/kreport2krona.py" -r result02_backraken2/$folder_name.brepo \
        -o result02_backraken2/$folder_name.krona.txt --no-intermediate-ranks

    ktImportText result02_backraken2/$folder_name.krona.txt \
        -o result02_backraken2/$folder_name.krona.html
    local end_seconds=$(date +"%s")
    
    echo "Plot finished!!!" >&3
    print_duration "$start_seconds" "$end_seconds" >&3    
    
    kraken2read=result02_backraken2/$folder_name.kraken2
    kraken2report=result02_backraken2/$folder_name.k2report
    
    
    #extract all bacteria reads
    local start_seconds=$(date +"%s")
    python "/storage/caishangLab/fangxiunan/tools/KrakenTools/extract_kraken_reads.py" -t 2 \
    --include-children -k $kraken2read -r $kraken2report \
    -s1 $unmappedread1 -s2 $unmappedread2 \
    -o ${folder_name}_kraken2_bac.1.fq -o2 ${folder_name}_kraken2_bac.2.fq
    
    local end_seconds=$(date +"%s")
    echo "Extract bac reads finished!!!" >&3
    print_duration "$start_seconds" "$end_seconds" >&3   
    
    
    local text="Kraken2 finished!!!" >&3
    local border=$(printf '#%.0s' $(seq 1 $((${#text} + 8)))) >&3
    echo "$border" >&3
    echo "# $text #" >&3
    print_duration "$zero_seconds" "$end_seconds" >&3
    echo "$border" >&3
    
    return 0
}

# Export the function to be used by GNU Parallel

export -f process_sample


process_file() {
    local file="$1"
    local destination_dir="$2"

    echo "Processing file: $(basename "$file")"
    if ! process_sample "$file" "$destination_dir"; then

        echo "Failed processing: $(basename "$file" | sed 's/_.*//')" >> $destination_dir/s3_failed_files.log

        echo "Failed processing: $(basename "$file" | sed 's/_.*//') - $(date)" >> $destination_dir/s3_debug.log

        return 1

    else

        echo "Successfully processed: $(basename "$file")" >> $destination_dir/s3_debug.log

    fi

    return 0

}

export -f process_file
# List the files to process

files=("$destination_dir"/*_unmapped.dedupe.1.fq)

for file in "${files[@]}"; do 
    process_file $file $destination_dir >> $destination_dir/s3joblog.txt 2>&1
done