#!/bin/bash

usage() {
    echo "Usage: $0 -s <source_dir> -n <files_per_folder>"
    exit 1

}


source_dir=""
files_per_folder=""



while getopts ":s:n:" opt; do

    case $opt in

        s) source_dir="$OPTARG"
        ;;
        n) files_per_folder="$OPTARG"
        ;;
        *) usage

        ;;
    esac

done



if [ -z "$source_dir" ] || [ -z "$files_per_folder" ]; then

    usage

fi

output_prefix="./splitted_folder"


files=($(ls $source_dir/*_R1_*.fastq.gz 2>/dev/null | sort))


folder_index=1

file_count=0

cd $source_dir

mkdir -p "${output_prefix}_${folder_index}"

for r1_file in "${files[@]}"; do


    r2_file="${r1_file/_R1_/_R2_}"

    cp "$r1_file" "${output_prefix}_${folder_index}/"
    cp "$r2_file" "${output_prefix}_${folder_index}/"

    file_count=$((file_count + 2))

    if [[ $file_count -ge $files_per_folder ]] && [[ ${#files[@]} -gt $file_count ]]; then

        folder_index=$((folder_index + 1))
        mkdir -p "${output_prefix}_${folder_index}"
        file_count=0

    fi

done

if [[ $file_count -eq 0 ]]; then
    rmdir "${output_prefix}_${folder_index}"
fi

echo "Split done！"


