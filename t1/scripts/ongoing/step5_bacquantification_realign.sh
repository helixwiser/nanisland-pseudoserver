#!/bin/bash
source /opt/anaconda3/etc/profile.d/conda.sh


usage() {
    echo "Usage: $0 -s source_dir -d destination_dir "
    exit 1

}

# Parse command-line arguments

while getopts ":s:d:" opt; do

    case $opt in

        s) source_dir="$OPTARG"
        ;;
        d) destination_dir="$OPTARG"
        ;;
        *) usage

        ;;
    esac

done


process_sample() {
    local file=$1
    local source_dir=$2
    local destination_dir=$3
    local folder_name=$(basename "$file" | sed 's/_.*//')
    echo "Processing with ID: $folder_name"
    echo "$source_dir"
    source /opt/anaconda3/etc/profile.d/conda.sh
    print_duration() {
        local start_seconds=$1

        local end_seconds=$2

        local duration=$((end_seconds - start_seconds))
        local duration_formatted=$(date -u -d @"$duration" +"%H:%M:%S")
        echo "Duration: $duration_formatted"
    }


    local trimdir=$destination_dir/$folder_name

    mkdir -p $trimdir

    cd $trimdir
    exec 3>> exec_time.log  
    mkdir -p $trimdir/result04_pathseq
    
    local bacread1=kraken2_bac.1.fq
    local bacread2=kraken2_bac.2.fq
    local zero_seconds=$(date +"%s")

#   ###################################
#   ######## pathseq anno ########
#   ###################################
    local start_seconds=$(date +"%s")
   java -jar /home/available_tools/picard-2.30.0/picard.jar FastqToSam \
     TMP_DIR=/mnt/resources/FXN/tmp \
     F1=$bacread1 \
     F2=$bacread2 \
     O=kraken2bac.bam \
     SM=hostunmapped001 \
     RG=hu \
     V=Standard
     
    local end_seconds=$(date +"%s")   
    echo "fastqtosam finished!!!" >&3
    print_duration "$start_seconds" "$end_seconds" >&3
    
    local start_seconds=$(date +"%s")    
   /mnt/shared/tools/gatk/gatk PathSeqPipelineSpark --java-options "-Xmx64G" \
    --tmp-dir /mnt/resources/FXN/tmp \
   	--input kraken2bac.bam \
   	--filter-bwa-image "/mnt/shared/reference/pathseq/mm10.fa.img" \
   	--kmer-file "/mnt/shared/reference/pathseq/mm10.hss" \
   	--min-clipped-read-length 70 \
   	--microbe-bwa-image "/mnt/shared/reference/pathseq/pathseq_microbe.fa.img" \
   	--taxonomy-file "/mnt/shared/reference/pathseq/pathseq_taxonomy.db" \
   	--microbe-dict "/mnt/shared/reference/pathseq/pathseq_microbe.dict" \
   	--output bacstat.pathseq.bam \
   	--scores-output $trimdir/result04_pathseq/bacstat.mm.pathseq.txt  
    local end_seconds=$(date +"%s")   
    echo "pathseq finished!!!" >&3
    print_duration "$start_seconds" "$end_seconds" >&3    


    local text="pathseq quantification finished!!!" >&3
    local border=$(printf '#%.0s' $(seq 1 $((${#text} + 8)))) >&3
    echo "$border" >&3
    echo "# $text #" >&3
    print_duration "$zero_seconds" "$end_seconds" >&3
    echo "$border" >&3
}

# Export the function to be used by GNU Parallel

export -f process_sample

# List the files to process

files=$(ls $source_dir/*_R1_001.fastq.gz)
# Use GNU Parallel to process each file in parallel
parallel --jobs 4 process_sample {} $source_dir $destination_dir ::: "${files[@]}"


