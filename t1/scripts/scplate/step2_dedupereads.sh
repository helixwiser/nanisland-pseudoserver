#!/bin/bash
source /home/caishangLab/fangxiunan/anaconda3/etc/profile.d/conda.sh
conda activate metagenomics

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
    
    module load samtools/1.16.1
    
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
    #mkdir -p $trimdir/result02_backraken2
    
    # Get file sizes in bytes

    local size_hm=$(wc -c < mapped_hm.1.fq)
    local size_mm=$(wc -c < mapped_mm.1.fq)

    # Compare sizes and assign host

    local host="human"
    if [ "$size_mm" -gt "$size_hm" ]; then

        host="mouse"
    fi
    local zero_seconds=$(date +"%s")
    
    ######################################
    ####### specify read1 and read2 ######
    ######################################
    if [ "$host" = "human" ]; then
    #-f unmapped When the input is coordinate-sorted, unmapped mates of mapped records and supplementary/secondary alignments are not marked as duplicates.  thus use bbmap to do the dedup
      samtools view -f 4 -b aligned_hm.sorted.rmdup.bam > unmapped.bam
    fi
    
    if [ "$host" = "mouse" ]; then
      samtools view -f 4 -b aligned_mm.sorted.rmdup.bam > unmapped.bam  
    fi    
    
    samtools fastq -1 unmapped.1.fq -2 unmapped.2.fq -0 unmapped.unpaired.fq unmapped.bam
    
    ##The “ac=f” flag disables containment removal.
    local start_seconds=$(date +"%s")
    "/storage/caishangLab/fangxiunan/tools/bbmap/dedupe.sh" ac=f maxedits=2 \
    in=unmapped.1.fq \
    out=${folder_name}_unmapped.dedupe.1.fq \
    overwrite==true

    "/storage/caishangLab/fangxiunan/tools/bbmap/dedupe.sh" ac=f maxedits=2 \
    in=unmapped.2.fq \
    out=${folder_name}_unmapped.dedupe.2.fq \
    overwrite==true   
    local end_seconds=$(date +"%s")
    
    echo "bbmap dedupe finished!!!" >&3
    print_duration "$start_seconds" "$end_seconds" >&3  
    
    
#    local unmappedread1=unmapped.dedupe.1.fq
#    local unmappedread2=unmapped.dedupe.2.fq
#    ######################################
#    ####### kraken2 quantification########
#    ######################################
#    #use confidence = 0.7
#    local start_seconds=$(date +"%s")
#    kraken2 \
#        --db /storage/caishangLab/caishang/reference/kraken2-ref/ \
#        --threads 32 \
#        --confidence 0.7 \
#        --report result02_backraken2/$folder_name.k2report \
#        --report-minimizer-data \
#        --minimum-hit-groups 3 \
#        --paired $unmappedread1 $unmappedread2\
#        > result02_backraken2/$folder_name.kraken2
#    local end_seconds=$(date +"%s")
#    
#    echo "kraken2 finished!!!" >&3
#    print_duration "$start_seconds" "$end_seconds" >&3
#    
#    local start_seconds=$(date +"%s")
#    bracken \
#        -d /storage/caishangLab/caishang/reference/kraken2-ref/ \
#        -i result02_backraken2/$folder_name.k2report \
#        -r 100 \
#        -l S \
#        -t 10 \
#        -o result02_backraken2/$folder_name.bracken \
#        -w result02_backraken2/$folder_name.brepo
#    local end_seconds=$(date +"%s")
#    
#    echo "braken finished!!!" >&3
#    print_duration "$start_seconds" "$end_seconds" >&3
#
#    local start_seconds=$(date +"%s")
#    python "/storage/caishangLab/fangxiunan/tools/KrakenTools/kreport2krona.py" -r result02_backraken2/$folder_name.brepo \
#        -o result02_backraken2/$folder_name.krona.txt --no-intermediate-ranks
#
#    ktImportText result02_backraken2/$folder_name.krona.txt \
#        -o result02_backraken2/$folder_name.krona.html
#    local end_seconds=$(date +"%s")
#    
#    echo "Plot finished!!!" >&3
#    print_duration "$start_seconds" "$end_seconds" >&3    


    local text="Dedupe finished!!!" >&3
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
    local source_dir="$2"
    local destination_dir="$3"

    echo "Processing file: $(basename "$file")"
    if ! process_sample "$file" "$source_dir" "$destination_dir"; then

        echo "Failed processing: $(basename "$file" | sed 's/_.*//')" >> $destination_dir/s2_failed_files.log

        echo "Failed processing: $(basename "$file" | sed 's/_.*//') - $(date)" >> $destination_dir/s2_debug.log

        return 1

    else

        echo "Successfully processed: $(basename "$file")" >> $destination_dir/s2_debug.log

    fi

    return 0

}

export -f process_file
# List the files to process

files=("$source_dir"/*_R1_001.fastq.gz)

for file in "${files[@]}"; do 
    process_file $file $source_dir $destination_dir >> $destination_dir/s2joblog.txt 2>&1
done