#!/bin/bash
source /home/caishangLab/fangxiunan/anaconda3/etc/profile.d/conda.sh

usage() {
    echo "Usage: $0 -s source_dir -d destination_dir "
    exit 1

}




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
    local read1=$(ls "$source_dir"/${folder_name}_*_R1_001.fastq.gz 2>/dev/null)
    local read2=$(ls "$source_dir"/${folder_name}_*_R2_001.fastq.gz 2>/dev/null)
    
    module load fastqc/0.11.9
    fastqc=/soft/bio/fastqc-0.11.9/fastqc
    fastqc="/mnt/shared/tools/fastqc-0.11.9/fastqc"
    trim_galore="/storage/caishangLab/caishang/tools/TrimGalore-0.6.7/trim_galore"
    cutadapt="/home/caishangLab/fangxiunan/anaconda3/envs/fastq_qc/bin/cutadapt"

    # Check if the files were found
 
    if [ -z "$read1" ]; then

        echo "No R1 file found in ${source_dir}/${folder_name}"
        return 1

    else

        echo "R1 file found: $read1"
    fi

    if [ -z "$read2" ]; then

        echo "No R2 file found in ${source_dir}/${folder_name}"
        return 1

    else

        echo "R2 file found: $read2"
    fi

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

        echo "Failed processing: $(basename "$file" | sed 's/_.*//')" >> $destination_dir/s0_failed_files.log

        echo "Failed processing: $(basename "$file" | sed 's/_.*//') - $(date)" >> $destination_dir/s0_debug.log

        return 1

    else

        echo "Successfully processed: $(basename "$file")" >> $destination_dir/0_debug.log

    fi

    return 0

}

export -f process_file
# List the files to process



files=("$source_dir"/*_R1_001.fastq.gz)
for file in "${files[@]}"; do 
    /usr/bin/time -v bash -c "process_file '$file' '$source_dir' '$destination_dir'" >> $destination_dir/s0joblog.txt 2>&1
done
