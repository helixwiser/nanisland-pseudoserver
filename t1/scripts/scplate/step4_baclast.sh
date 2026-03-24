#!/bin/bash

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


export destination_dir

process_sample() {
    local file=$1
    local destination_dir=$2
    local folder_name=$(basename "$file" | sed 's/_.*//')
    echo "Processing with ID: $folder_name"
    source /home/caishangLab/fangxiunan/anaconda3/etc/profile.d/conda.sh
    conda activate metagenomics
    seqtk=/storage/caishangLab/caishang/tools/seqtk/seqtk
    
    
    print_duration() {
        local start_seconds=$1
        local end_seconds=$2

        local duration=$((end_seconds - start_seconds))
        local duration_formatted=$(date -u -d @"$duration" +"%H:%M:%S")
        echo "Duration: $duration_formatted"
    }



    mkdir -p $destination_dir/result03_bacdetail

    cd $destination_dir/result03_bacdetail
    exec 3>> exec_time.log  

    
    local bacread1=$destination_dir/${folder_name}_kraken2_bac.1.fq
    local bacread2=$destination_dir/${folder_name}_kraken2_bac.2.fq
    
    
    
    local zero_seconds=$(date +"%s")
    #extract all bacteria reads
    
    
#    local start_seconds=$(date +"%s")
#    
#    refdbtrembl=/storage/caishangLab/caishang/reference/xiunan/protein/uniprot_trembl_bacteria_last
#    
#    
#    $seqtk seq -a $bacread1 > ${folder_name}_kraken2_bac.1.fasta
#    last-train --codon $refdbtrembl ${folder_name}_kraken2_bac.1.fasta > ${folder_name}.1.codon.trembl.train
#    lastal -p ${folder_name}.1.codon.trembl.train -m100 -D1e9 -K1 $refdbtrembl ${folder_name}_kraken2_bac.1.fasta > ${folder_name}.1.out.trembl.maf
#    
#    $seqtk seq -a $bacread2 > ${folder_name}_kraken2_bac.2.fasta
#    last-train --codon $refdbtrembl ${folder_name}_kraken2_bac.2.fasta > ${folder_name}.2.codon.trembl.train
#    lastal -p ${folder_name}.2.codon.trembl.train -m100 -D1e9 -K1 $refdbtrembl ${folder_name}_kraken2_bac.2.fasta > ${folder_name}.2.out.trembl.maf    
#    
#    local end_seconds=$(date +"%s")
#    
#    
#    echo "last uniprot_trembl finished!!!" >&3
#    print_duration "$start_seconds" "$end_seconds" >&3
#    
    
    refdbsprot=/storage/caishangLab/caishang/reference/xiunan/protein/uniprot_sprot_last
    
    local start_seconds=$(date +"%s")
    $seqtk seq -a $bacread1 > ${folder_name}_kraken2_bac.1.fasta
    last-train --codon $refdbsprot ${folder_name}_kraken2_bac.1.fasta > ${folder_name}.1.codon.sprot.train
    lastal -p ${folder_name}.1.codon.sprot.train -m100 -D1e9 -K1 $refdbsprot ${folder_name}_kraken2_bac.1.fasta > ${folder_name}.1.out.sprot.maf
    
    $seqtk seq -a $bacread2 > ${folder_name}_kraken2_bac.2.fasta
    last-train --codon $refdbsprot ${folder_name}_kraken2_bac.2.fasta > ${folder_name}.2.codon.sprot.train
    lastal -p ${folder_name}.2.codon.sprot.train -m100 -D1e9 -K1 $refdbsprot ${folder_name}_kraken2_bac.2.fasta > ${folder_name}.2.out.sprot.maf    
    
    local end_seconds=$(date +"%s")
    
    
    echo "last uniprot_sprot finished!!!" >&3
    print_duration "$start_seconds" "$end_seconds" >&3
    
#########################################################
    local text="Bac reads protein mapping finished!!!" >&3
    local border=$(printf '#%.0s' $(seq 1 $((${#text} + 8)))) >&3
    echo "$border" >&3
    echo "# $text #" >&3
    print_duration "$zero_seconds" "$end_seconds" >&3
    echo "$border" >&3
    
    return 0
}

export -f process_sample


process_file() {
    local file="$1"
    destination_dir="$2"
    echo "Processing file: $(basename "$file")"
    if ! process_sample "$file" "$destination_dir"; then

        echo "Failed processing: $(basename "$file" | sed 's/_.*//')" >> $destination_dir/s4_failed_files.log

        echo "Failed processing: $(basename "$file" | sed 's/_.*//') - $(date)" >> $destination_dir/s4_debug.log

        return 1

    else

        echo "Successfully processed: $(basename "$file")" >> $destination_dir/s4_debug.log

    fi

    return 0

}

export -f process_file
# List the files to process
files=("$destination_dir"/*_kraken2_bac.1.fq)
echo 
# Use GNU Parallel to process each file in parallel

for file in "${files[@]}"; do 
    process_file $file $destination_dir >> $destination_dir/s4joblog.txt 2>&1
done
#usage() {
#    echo "Usage: $0 -s source_dir -d destination_dir "
#    exit 1
#
#}
#
## Parse command-line arguments
#
#while getopts ":s:d:" opt; do
#
#    case $opt in
#
#        s) source_dir="$OPTARG"
#        ;;
#        d) destination_dir="$OPTARG"
#        ;;
#        *) usage
#
#        ;;
#    esac
#
#done

#export source_dir
#
#export destination_dir
#
#process_sample() {
#    local file=$1
#    local source_dir=$2
#    local destination_dir=$3
#    local folder_name=$(basename "$file" | sed 's/_.*//')
#    echo "Processing with ID: $folder_name"
#    source /home/caishangLab/fangxiunan/anaconda3/etc/profile.d/conda.sh
#    conda activate metagenomics
#    seqtk=/storage/caishangLab/caishang/tools/seqtk/seqtk
#    
#    
#    print_duration() {
#        local start_seconds=$1
#
#        local end_seconds=$2
#
#        local duration=$((end_seconds - start_seconds))
#        local duration_formatted=$(date -u -d @"$duration" +"%H:%M:%S")
#        echo "Duration: $duration_formatted"
#    }
#
#
#
#    mkdir -p $destination_dir/result03_bacdetail
#
#    cd $destination_dir/result03_bacdetail
#    exec 3>> exec_time.log  
#
#    
#    local bacread1=$source_dir/${folder_name}_kraken2_bac.1.fq
#    local bacread2=$source_dir/${folder_name}_kraken2_bac.2.fq
#    
#    
#    
#    local zero_seconds=$(date +"%s")
#    #extract all bacteria reads
#    
#    
##    local start_seconds=$(date +"%s")
##    
##    refdbtrembl=/storage/caishangLab/caishang/reference/xiunan/protein/uniprot_trembl_bacteria_last
##    
##    
##    $seqtk seq -a $bacread1 > ${folder_name}_kraken2_bac.1.fasta
##    last-train --codon $refdbtrembl ${folder_name}_kraken2_bac.1.fasta > ${folder_name}.1.codon.trembl.train
##    lastal -p ${folder_name}.1.codon.trembl.train -m100 -D1e9 -K1 $refdbtrembl ${folder_name}_kraken2_bac.1.fasta > ${folder_name}.1.out.trembl.maf
##    
##    $seqtk seq -a $bacread2 > ${folder_name}_kraken2_bac.2.fasta
##    last-train --codon $refdbtrembl ${folder_name}_kraken2_bac.2.fasta > ${folder_name}.2.codon.trembl.train
##    lastal -p ${folder_name}.2.codon.trembl.train -m100 -D1e9 -K1 $refdbtrembl ${folder_name}_kraken2_bac.2.fasta > ${folder_name}.2.out.trembl.maf    
##    
##    local end_seconds=$(date +"%s")
##    
##    
##    echo "last uniprot_trembl finished!!!" >&3
##    print_duration "$start_seconds" "$end_seconds" >&3
##    
#    
#    refdbsprot=/storage/caishangLab/caishang/reference/xiunan/protein/uniprot_sprot_last
#    
#    local start_seconds=$(date +"%s")
#    $seqtk seq -a $bacread1 > ${folder_name}_kraken2_bac.1.fasta
#    last-train --codon $refdbsprot ${folder_name}_kraken2_bac.1.fasta > ${folder_name}.1.codon.sprot.train
#    lastal -p ${folder_name}.1.codon.sprot.train -m100 -D1e9 -K1 $refdbsprot ${folder_name}_kraken2_bac.1.fasta > ${folder_name}.1.out.sprot.maf
#    
#    $seqtk seq -a $bacread2 > ${folder_name}_kraken2_bac.2.fasta
#    last-train --codon $refdbsprot ${folder_name}_kraken2_bac.2.fasta > ${folder_name}.2.codon.sprot.train
#    lastal -p ${folder_name}.2.codon.sprot.train -m100 -D1e9 -K1 $refdbsprot ${folder_name}_kraken2_bac.2.fasta > ${folder_name}.2.out.sprot.maf    
#    
#    local end_seconds=$(date +"%s")
#    
#    
#    echo "last uniprot_sprot finished!!!" >&3
#    print_duration "$start_seconds" "$end_seconds" >&3
#    
##########################################################
#    local text="Bac reads protein mapping finished!!!" >&3
#    local border=$(printf '#%.0s' $(seq 1 $((${#text} + 8)))) >&3
#    echo "$border" >&3
#    echo "# $text #" >&3
#    print_duration "$zero_seconds" "$end_seconds" >&3
#    echo "$border" >&3
#    
#    return 0
#}
#
#export -f process_sample
#
#
#process_file() {
#    local file="$1"
#    source_dir="$2"
#    destination_dir="$3"
#    echo "Processing file: $(basename "$file")"
#    if ! process_sample "$file" "$source_dir" "$destination_dir"; then
#
#        echo "Failed processing: $(basename "$file" | sed 's/_.*//')" >> $destination_dir/s4_failed_files.log
#
#        echo "Failed processing: $(basename "$file" | sed 's/_.*//') - $(date)" >> $destination_dir/s4_debug.log
#
#        return 1
#
#    else
#
#        echo "Successfully processed: $(basename "$file")" >> $destination_dir/s4_debug.log
#
#    fi
#
#    return 0
#
#}
#
#export -f process_file
## List the files to process
#files=("$source_dir"/*_kraken2_bac.1.fq)
#echo 
## Use GNU Parallel to process each file in parallel
#
#for file in "${files[@]}"; do 
#    process_file $file $source_dir $destination_dir >> $destination_dir/s4joblog.txt 2>&1
#done