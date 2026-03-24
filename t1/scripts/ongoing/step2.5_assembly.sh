#!/bin/bash


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

export source_dir

export destination_dir

process_sample() {
    local file=$1
    local source_dir=$2
    local destination_dir=$3
    local folder_name=$(basename "$file" | sed 's/_.*//')
    echo "Processing with ID: $folder_name"
    source /home/caishangLab/fangxiunan/anaconda3/etc/profile.d/conda.sh
    conda activate metagenomics
    seqtk=/storage/caishangLab/caishang/tools/seqtk/seqtk
    fastp=/storage/caishangLab/fangxiunan/tools/fastp
    
    print_duration() {
        local start_seconds=$1

        local end_seconds=$2

        local duration=$((end_seconds - start_seconds))
        local duration_formatted=$(date -u -d @"$duration" +"%H:%M:%S")
        echo "Duration: $duration_formatted"
    }


    mkdir -p $destination_dir/result03_assemble

    cd $destination_dir
    exec 3>> exec_time.log  
    local zero_seconds=$(date +"%s")
    exec > ${folder_name}.assemble.log 2>&1
    
    local start_seconds=$(date +"%s")
    $fastp -i ${folder_name}_unmapped.dedupe.1.fq -o ${folder_name}_unmapped.dedupe.trimmed.1.fq \
    --cut_front 10 --cut_tail 10 --cut_window_size 5 --cut_mean_quality 10 --qualified_quality_phred 15 \
    --unqualified_percent_limit 50 --length_required 30 --trim_poly_g --poly_g_min_len 5 \
    --trim_poly_x --poly_x_min_len 5 -y 
    
    $fastp -i ${folder_name}_unmapped.dedupe.2.fq -o ${folder_name}_unmapped.dedupe.trimmed.2.fq \
    --cut_front 10 --cut_tail 10 --cut_window_size 5 --cut_mean_quality 10 --qualified_quality_phred 15 \
    --unqualified_percent_limit 50 --length_required 30 --trim_poly_g --poly_g_min_len 5 \
    --trim_poly_x --poly_x_min_len 5 -y 
    
    local end_seconds=$(date +"%s")
    echo "trim_poly_g finished!!!" >&3
    print_duration "$start_seconds" "$end_seconds" >&3   
    
    local start_seconds=$(date +"%s")
    megahit -1 ${folder_name}_unmapped.dedupe.trimmed.1.fq \
    -2 ${folder_name}_unmapped.dedupe.trimmed.2.fq \
    -o result03_assemble/${folder_name}
    local end_seconds=$(date +"%s")
    echo "megahit finished!!!" >&3
    print_duration "$start_seconds" "$end_seconds" >&3   
    
        
    local unmappedread=result03_assemble/${folder_name}/final.contigs.fa
    local start_seconds=$(date +"%s")
    kraken2 \
        --db /storage/caishangLab/caishang/reference/kraken2-ref/ \
        --threads 32 \
        --confidence 0.7 \
        --report result03_assemble/${folder_name}.k2report \
        --report-minimizer-data \
        --minimum-hit-groups 3 \
        $unmappedread \
        > result03_assemble/${folder_name}.kraken2


    bracken \
        -d /storage/caishangLab/caishang/reference/kraken2-ref/ \
        -i result03_assemble/${folder_name}.k2report \
        -r 100 \
        -l S \
        -t 32 \
        -o result03_assemble/${folder_name}.bracken \
        -w result03_assemble/${folder_name}.brepo
    
    python "/storage/caishangLab/fangxiunan/tools/KrakenTools/kreport2krona.py" -r result03_assemble/${folder_name}.k2report \
        -o result03_assemble/${folder_name}.krona1.txt --no-intermediate-ranks
    
    
    ktImportText result03_assemble/${folder_name}.krona1.txt \
        -o result03_assemble/${folder_name}.krona1.html
        
        
    python "/storage/caishangLab/fangxiunan/tools/KrakenTools/extract_kraken_reads.py" -t 2 \
    --include-children -k result03_assemble/${folder_name}.kraken2 -r result03_assemble/${folder_name}.k2report \
    -s $unmappedread \
    -o result03_assemble/${folder_name}_kraken2_bac.fq 
    local end_seconds=$(date +"%s")
    echo "kraken2 for contig finished!!!" >&3
    print_duration "$start_seconds" "$end_seconds" >&3  
      
    cd result03_assemble
    bacread=${folder_name}_kraken2_bac.fq 
    refdbsprot=/storage/caishangLab/caishang/reference/xiunan/protein/uniprot_sprot_last
    
    local start_seconds=$(date +"%s")
    $seqtk seq -a $bacread > ${folder_name}_kraken2_bac.fasta
    last-train --codon $refdbsprot ${folder_name}_kraken2_bac.fasta > ${folder_name}.codon.sprot.train
    lastal -p ${folder_name}.codon.sprot.train -m100 -D1e9 -K1 $refdbsprot ${folder_name}_kraken2_bac.fasta > ${folder_name}.out.sprot.maf    
    local end_seconds=$(date +"%s")
    echo "contig protein mapping finished!!!" >&3
    print_duration "$start_seconds" "$end_seconds" >&3     

 
    
#########################################################
    local text="Bac contig reads protein mapping finished!!!" >&3
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
    source_dir="$2"
    destination_dir="$3"
    echo "Processing file: $(basename "$file")"
    if ! process_sample "$file" "$source_dir" "$destination_dir"; then

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
files=("$source_dir"/*_unmapped.dedupe.1.fq)
echo 
# Use GNU Parallel to process each file in parallel

for file in "${files[@]}"; do 
    process_file $file $source_dir $destination_dir >> $destination_dir/s4joblog.txt 2>&1
done