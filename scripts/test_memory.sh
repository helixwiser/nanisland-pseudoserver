#!/bin/bash

#SBATCH -p intel-sc3,amd-ep2,intel-sc3-32c
#SBATCH -q normal
#SBATCH -J splitted
#SBATCH -c 16
#SBATCH --output=/storage/caishangLab/fangxiunan/log/exec_%A_%a.log
#SBATCH --error=/storage/caishangLab/fangxiunan/log/exec_%A_%a.err
#SBATCH --array=0-$(($num_folders - 1))
#SBATCH --mem=64G


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
    local read1=$(ls "$source_dir"/${folder_name}_*_R1_001.fastq.gz 2>/dev/null)
    local read2=$(ls "$source_dir"/${folder_name}_*_R2_001.fastq.gz 2>/dev/null)
    
    module load fastqc/0.11.9
    module load picard/2.25.1
    module load samtools/1.16.1
    module load picard/2.25.1
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

    print_duration() {
        local start_seconds=$1

        local end_seconds=$2

        local duration=$((end_seconds - start_seconds))
        local duration_formatted=$(date -u -d @"$duration" +"%H:%M:%S")
        echo "Duration: $duration_formatted"
    }

    ######################################
    ########### Trim and process #########
    ######################################
    local trimdir=$destination_dir/$folder_name

    mkdir -p $trimdir

    cd $trimdir
    exec 3>> exec_time.log  
    local zero_seconds=$(date +"%s")
    echo "Step: Trimming ATAC"

    ######################################
    ########### specify read1 and read2 #########
    ######################################
    cp $read1 "$folder_name.1.fastq.gz"
    cp $read2 "$folder_name.2.fastq.gz"
    "$trim_galore" --paired -q 20 --path_to_cutadapt $cutadapt -o "$trimdir" "$folder_name.1.fastq.gz" "$folder_name.2.fastq.gz"

    fastqc "$folder_name.1.fastq.gz" -d "$trimdir" -t 20 --noextract -o "$trimdir"
    fastqc "$folder_name.2.fastq.gz" -d "$trimdir" -t 20 --noextract -o "$trimdir"
    
    local readt1=$folder_name.1_val_1.fq.gz
    local readt2=$folder_name.2_val_2.fq.gz

    ######################################
    ########### Alignment host############

    ######################################
    echo "Aligning human reads to the human reference genome (GRCh38)..."
    local start_seconds=$(date +"%s")
    bowtie2 -p 32 --very-fast \
            -X 1000 \
            -x /storage/caishangLab/caishang/reference/xiunan/bowtie2/GRCh38_noalt_as/GRCh38_noalt_as \
            -1 $readt1 \
            -2 $readt2 \
            --un-conc unmapped_hm.fq --al-conc mapped_hm.fq -S aligned_hm.sam

    local end_seconds=$(date +"%s")

    echo "Aligning to human reads finished!!!" >&3
    print_duration "$start_seconds" "$end_seconds" >&3
 
    # Aligning unmapped human reads to the mouse reference genome (mm10)
    echo "Aligning unmapped human reads to the mouse reference genome (mm10)..."
    start_seconds=$(date +"%s")
    bowtie2 -p 32 --very-fast \
            -X 1000 \
            -x /storage/caishangLab/caishang/reference/xiunan/bowtie2/mm10/mm10 \
            -1 unmapped_hm.1.fq \
            -2 unmapped_hm.2.fq \
            --un-conc unmapped_mm.fq --al-conc mapped_mm.fq -S aligned_mm.sam

    end_seconds=$(date +"%s") 
    echo "Aligning to mouse reads finished!!!" >&3
    print_duration "$start_seconds" "$end_seconds" >&3

    # Get file sizes in bytes

    local size_hm=$(wc -c < mapped_hm.1.fq)
    local size_mm=$(wc -c < mapped_mm.1.fq)

    # Compare sizes and assign host

    local host="human"
    if [ "$size_mm" -gt "$size_hm" ]; then

        host="mouse"
    fi

    ###################################
    ################ stats ############

    ###################################

    if [ "$host" = "human" ]; then
        start_seconds=$(date +"%s")
        samtools view -S -b aligned_hm.sam > aligned_hm.bam
        samtools sort -o aligned_hm.sorted.bam aligned_hm.bam
        samtools index aligned_hm.sorted.bam

        samtools stats aligned_hm.sorted.bam | grep ^MAPQ | cut -f 2- > stats.mapq.aligned_hm.log
        samtools flagstats aligned_hm.sorted.bam > stats.reads.aligned_hm.log

        samtools view -Sb -F 4 -q 0 -f 0x2 -b aligned_hm.sorted.bam > mapped_hm.bam
        samtools stats mapped_hm.bam > stats.mapped_hm.log
        samtools coverage mapped_hm.bam > stats.cover.mapped_hm.log
        
        end_seconds=$(date +"%s") 
        echo "samtools stats finished!!!" >&3
        print_duration "$start_seconds" "$end_seconds" >&3
        
        start_seconds=$(date +"%s")
        
        picard MarkDuplicates \
        TMP_DIR=/storage/caishangLab/fangxiunan/tmp \
        INPUT=aligned_hm.sorted.bam \
        OUTPUT=aligned_hm.sorted.rmdup.bam \
        M="aligned_hm_marked_dup_metrics.txt" \
        VALIDATION_STRINGENCY=LENIENT ASSUME_SORTED=true REMOVE_DUPLICATES=TRUE 
        
        end_seconds=$(date +"%s") 
        echo "picard markdup finished!!!" >&3
        print_duration "$start_seconds" "$end_seconds" >&3
        
        start_seconds=$(date +"%s")   
        samtools flagstats aligned_hm.sorted.rmdup.bam > stats.reads.aligned_hm.rmdup.log
        samtools view -Sb -F 4 -q 20 -f 0x2 -b aligned_hm.sorted.rmdup.bam > mapped_hm.rmdup.bam
        samtools index mapped_hm.rmdup.bam
        samtools stats  mapped_hm.rmdup.bam > stats.mapped_hm.rmdup.log
        samtools coverage  mapped_hm.rmdup.bam > stats.cover.mapped_hm.rmdup.log
        samtools flagstats mapped_hm.rmdup.bam > stats.reads.mapped_hm.rmdup.log
        end_seconds=$(date +"%s") 
        echo "samtools (stats mapped dedup) finished!!!" >&3
        print_duration "$start_seconds" "$end_seconds" >&3
        
        
        #bamPEFragmentSize -b mapped_hm.rmdup.bam -o  mapped_hm.hist1000_length.log.pdf  --logScale  --binSize 1000 --maxFragmentLength 2000

        #################################Extract information to a txt file#################################
        echo $folder_name > summary.reads.txt

        grep -i "in total" stats.reads.aligned_hm.log | cut -d ' ' -f 1 >> summary.reads.txt

        echo $host >> summary.reads.txt

        grep -i "primary mapped" stats.reads.aligned_hm.log | cut -d ' ' -f 1 >> summary.reads.txt

        grep -i "primary mapped" stats.reads.mapped_hm.rmdup.log | cut -d ' ' -f 1 >> summary.reads.txt

        head -n 26 stats.cover.mapped_hm.rmdup.log |awk '$1 != "chrM" {sum += $6; count++} END {if (count > 0) print sum / count}'  >> summary.reads.txt

        grep -i "chrM" stats.cover.mapped_hm.log|cut -f 4 >> summary.reads.txt

        grep -i "chrM" stats.cover.mapped_hm.rmdup.log|cut -f 4 >> summary.reads.txt

        grep -i "chrM" stats.cover.mapped_hm.rmdup.log|cut -f 6 >> summary.reads.txt

        paste -d '\t' "/storage/caishangLab/fangxiunan/scripts/scplate/hostalignment_index.txt" summary.reads.txt > ${folder_name}_summary.hostreads.tsv

    fi

    if [ "$host" = "mouse" ]; then
        start_seconds=$(date +"%s")
        samtools view -S -b aligned_mm.sam > aligned_mm.bam
        samtools sort -o aligned_mm.sorted.bam aligned_mm.bam
        samtools index aligned_mm.sorted.bam

        samtools stats aligned_mm.sorted.bam | grep ^MAPQ | cut -f 2- > stats.mapq.aligned_mm.log
        samtools flagstats aligned_mm.sorted.bam > stats.reads.aligned_mm.log

        samtools view -Sb -F 4 -q 0 -f 0x2 -b aligned_mm.sorted.bam > mapped_mm.bam
        samtools stats mapped_mm.bam > stats.mapped_mm.log
        samtools coverage mapped_mm.bam > stats.cover.mapped_mm.log

        end_seconds=$(date +"%s") 
        echo "samtools stats finished!!!" >&3
        print_duration "$start_seconds" "$end_seconds" >&3


        start_seconds=$(date +"%s")
        picard MarkDuplicates \
        TMP_DIR=/storage/caishangLab/fangxiunan/tmp \
        INPUT=aligned_mm.sorted.bam \
        OUTPUT=aligned_mm.sorted.rmdup.bam \
        M="aligned_mm_marked_dup_metrics.txt" \
        VALIDATION_STRINGENCY=LENIENT ASSUME_SORTED=true REMOVE_DUPLICATES=TRUE 
        
        end_seconds=$(date +"%s") 
        echo "picard markdup finished!!!" >&3
        print_duration "$start_seconds" "$end_seconds" >&3
        
        
        start_seconds=$(date +"%s")        
        samtools flagstats aligned_mm.sorted.rmdup.bam > stats.reads.aligned_mm.rmdup.log
        samtools view -Sb -F 4 -q 20 -f 0x2 -b aligned_mm.sorted.rmdup.bam > mapped_mm.rmdup.bam      
        samtools index mapped_mm.rmdup.bam
        samtools stats  mapped_mm.rmdup.bam > stats.mapped_mm.rmdup.log
        samtools coverage  mapped_mm.rmdup.bam > stats.cover.mapped_mm.rmdup.log
        samtools flagstats mapped_mm.rmdup.bam > stats.reads.mapped_mm.rmdup.log
        
        end_seconds=$(date +"%s") 
        echo "samtools (stats mapped dedup) finished!!!" >&3
        print_duration "$start_seconds" "$end_seconds" >&3



        #bamPEFragmentSize -b mapped_mm.rmdup.bam -o  mapped_mm.hist1000_length.log.pdf  --logScale  --binSize 1000 --maxFragmentLength 2000


        #################################Extract information to a txt file#################################
        echo $folder_name > summary.reads.txt

        grep -i "in total" stats.reads.aligned_mm.log | cut -d ' ' -f 1 >> summary.reads.txt

        echo $host >> summary.reads.txt

        grep -i "primary mapped" stats.reads.aligned_mm.log | cut -d ' ' -f 1 >> summary.reads.txt

        grep -i "primary mapped" stats.reads.mapped_mm.rmdup.log | cut -d ' ' -f 1 >> summary.reads.txt

        head -n 23 stats.cover.mapped_mm.log |awk '$1 != "chrM" {sum += $6; count++} END {if (count > 0) print sum / count}' >> summary.reads.txt

        grep -i "chrM" stats.cover.mapped_mm.log|cut -f 4 >> summary.reads.txt

        grep -i "chrM" stats.cover.mapped_mm.rmdup.log|cut -f 4 >> summary.reads.txt

        grep -i "chrM" stats.cover.mapped_mm.rmdup.log|cut -f 6 >> summary.reads.txt

        paste -d '\t' "/storage/caishangLab/fangxiunan/scripts/scplate/hostalignment_index.txt" summary.reads.txt > ${folder_name}_summary.hostreads.tsv

    fi

    mkdir -p $trimdir/result01_hostalignment

    cp ${folder_name}_summary.hostreads.tsv $trimdir/result01_hostalignment/.

    local text="Aligning to host finished!!!" >&3
    local border=$(printf '#%.0s' $(seq 1 $((${#text} + 8)))) >&3
    echo "$border" >&3
    echo "# $text #" >&3
    print_duration "$zero_seconds" "$end_seconds" >&3
    echo "$border" >&3
    
    return 0
}

# Export the function to be used by GNU Parallel
#
#export -f process_sample
#
#
#process_file() {
#    file="$1"
#    source_dir="$2"
#    destination_dir="$3"
#
#    if ! process_sample "$file" "$source_dir" "$destination_dir"; then
#
#        echo "Failed processing: $(basename "$file" | sed 's/_.*//')" >> s1_failed_files.log
#
#    fi
#    
#    return 1
#
#}
#
#export -f process_file
## List the files to process
#
#files=$(ls $source_dir/*_R1_001.fastq.gz)
## Use GNU Parallel to process each file in parallel
#
#
#parallel --jobs 4 process_file {} $source_dir $destination_dir ::: "${files[@]}"

export -f process_sample


process_file() {
    local file="$1"
    local source_dir="$2"
    local destination_dir="$3"

    echo "Processing file: $(basename "$file")"
    if ! process_sample "$file" "$source_dir" "$destination_dir"; then

        echo "Failed processing: $(basename "$file" | sed 's/_.*//')" >> $destination_dir/s1_failed_files.log

        echo "Failed processing: $(basename "$file" | sed 's/_.*//') - $(date)" >> $destination_dir/s1_debug.log

        return 1

    else

        echo "Successfully processed: $(basename "$file")" >> $destination_dir/s1_debug.log

    fi

    return 0

}

export -f process_file
# List the files to process



files=("$source_dir"/*_R1_001.fastq.gz)
for file in "${files[@]}"; do 
    /usr/bin/time -v bash -c "process_file '$file' '$source_dir' '$destination_dir'" >> memory_log.txt 2>&1
done
# Use GNU Parallel to process each file in parallel
#files=$(ls $source_dir/*_R1_001.fastq.gz)
#for file in $files; do
#   process_file "$file" "$source_dir" "$destination_dir" >> s0joblog.txt
#done 

#for file in $files; do 
#   /usr/bin/time -v process_file "$file" "$source_dir" "$destination_dir" >> s0joblog.txt 2>&1
#done

#parallel --tag 'echo {#} {} >> pid.log' --verbose --joblog s0joblog.txt  --jobs 4 process_file {} $source_dir $destination_dir ::: "${files[@]}"

