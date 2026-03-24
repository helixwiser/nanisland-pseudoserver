#!/bin/bash

#SBATCH -p intel-sc3,amd-ep2,intel-sc3-32c
#SBATCH -q normal
#SBATCH --cpus-per-task=20
#SBATCH --output=/storage/caishangLab/fangxiunan/log/exec_%A_%a.log
#SBATCH --error=/storage/caishangLab/fangxiunan/log/exec_%A_%a.err
#SBATCH --mem=200G


source /home/caishangLab/fangxiunan/anaconda3/etc/profile.d/conda.sh
conda activate metawrap-env
bbmap=/storage/caishangLab/fangxiunan/tools/bbmap

#########change here
inputdir=/storage/caishangLab/fangxiunan/project-ongoing/metagenomic/batch-20240221
#inputdir=/storage/caishangLab/fangxiunan/project-rawdata/metagenomic/batch-test
#########change here


type=cancer
resdir=$inputdir/result04_assembly/$type
log_file=merge_H_${type}.log

#####################################################
#1.rename  .1.fq into _1.fq
mkdir -p $resdir
cd $resdir

for file in $inputdir/*-${type}-skin_unmapped.dedupe.1.fq; do
    cp "$file" "${file/.1.fq/_1.fastq}"
done

for file in $inputdir/*-${type}-skin_unmapped.dedupe.2.fq; do
    cp "$file" "${file/.2.fq/_2.fastq}"
done

mkdir -p $resdir/RAW_READS


mv $inputdir/*.fastq $resdir/RAW_READS/.



mkdir $resdir/CLEAN_READS
for F in RAW_READS/*_1.fastq; do 
  echo $F
	R=${F%_*}_2.fastq
  echo $R
	BASE=${F##*/}
	SAMPLE=${BASE%_*}

  $bbmap/repair.sh \
  in1=$F in2=$R out1=$SAMPLE.repair_1.fq out2=$SAMPLE.repair_2.fq outs=$SAMPLE.singletons.fq repair
  mv $SAMPLE.repair_1.fq CLEAN_READS/${SAMPLE}_1.fastq
  mv $SAMPLE.repair_2.fq CLEAN_READS/${SAMPLE}_2.fastq
	#metawrap read_qc --skip-trimming -1 $SAMPLE.repair_1.fq -2 $SAMPLE.repair_2.fq -t 32 -o READ_QC/$SAMPLE & 
done	


cat CLEAN_READS/*_1.fastq > CLEAN_READS/ALL_READS_1.fastq
cat CLEAN_READS/*_2.fastq > CLEAN_READS/ALL_READS_2.fastq

metawrap assembly -1 CLEAN_READS/ALL_READS_1.fastq -2 CLEAN_READS/ALL_READS_2.fastq \
-m 200 -t 96 --metaspades \
-o ASSEMBLY

metawrap binning -a ASSEMBLY/final_assembly.fasta \
-t 96 --metabat2 --maxbin2 --concoct \
-o INITIAL_BINNING \
CLEAN_READS/H*fastq

metawrap bin_refinement -A INITIAL_BINNING/metabat2_bins/ -B INITIAL_BINNING/maxbin2_bins/ -C INITIAL_BINNING/concoct_bins/ \
-t 96 -c 50 -x 10 \
-o BIN_REFINEMENT

metawrap blobology -a ASSEMBLY/final_assembly.fasta \
-t 96 --bins BIN_REFINEMENT/metawrap_50_10_bins \
-o BLOBOLOGY \
CLEAN_READS/H*fastq

metawrap quant_bins -a ASSEMBLY/final_assembly.fasta -b BIN_REFINEMENT/metawrap_50_10_bins \
-o QUANT_BINS \
CLEAN_READS/H*fastq
 

metawrap reassemble_bins -1 CLEAN_READS/ALL_READS_1.fastq -2 CLEAN_READS/ALL_READS_2.fastq -b BIN_REFINEMENT/metawrap_50_10_bins \
-t 96 -m 800 -c 50 -x 10 \
-o BIN_REASSEMBLY

metawrap classify_bins -b BIN_REASSEMBLY/reassembled_bins \
-t 48 \
-o BIN_CLASSIFICATION 


mkdir $resdir/BIN_REASSEMBLY/reassembled_bins_renamedcontig

cd $resdir/BIN_REASSEMBLY
for fa in "$resdir/BIN_REASSEMBLY/reassembled_bins/"*.fa; do

  [ -f "$fa" ] || continue
  origname=$(basename "$fa")
  filename="${origname%.fa}"
  
  $bbmap/rename.sh \
  in=$fa \
  out=reassembled_bins_renamedcontig/$filename.renamedcontig.fa \
  prefix=contig
  
done

cd $resdir
metaWRAP annotate_bins -o FUNCT_ANNOT -t 96 -b BIN_REASSEMBLY/reassembled_bins_renamedcontig/

######################################################################################
type=normal
resdir=$inputdir/result04_assembly/$type
log_file=merge_H_${type}.log

#####################################################
#1.rename  .1.fq into _1.fq
mkdir -p $resdir
cd $resdir

for file in $inputdir/*-${type}-skin_unmapped.dedupe.1.fq; do
    cp "$file" "${file/.1.fq/_1.fastq}"
done

for file in $inputdir/*-${type}-skin_unmapped.dedupe.2.fq; do
    cp "$file" "${file/.2.fq/_2.fastq}"
done

mkdir -p $resdir/RAW_READS


mv $inputdir/*.fastq $resdir/RAW_READS/.



mkdir $resdir/CLEAN_READS
for F in RAW_READS/*_1.fastq; do 
  echo $F
	R=${F%_*}_2.fastq
  echo $R
	BASE=${F##*/}
	SAMPLE=${BASE%_*}

  $bbmap/repair.sh \
  in1=$F in2=$R out1=$SAMPLE.repair_1.fq out2=$SAMPLE.repair_2.fq outs=$SAMPLE.singletons.fq repair
  mv $SAMPLE.repair_1.fq CLEAN_READS/${SAMPLE}_1.fastq
  mv $SAMPLE.repair_2.fq CLEAN_READS/${SAMPLE}_2.fastq
	#metawrap read_qc --skip-trimming -1 $SAMPLE.repair_1.fq -2 $SAMPLE.repair_2.fq -t 32 -o READ_QC/$SAMPLE & 
done	


cat CLEAN_READS/*_1.fastq > CLEAN_READS/ALL_READS_1.fastq
cat CLEAN_READS/*_2.fastq > CLEAN_READS/ALL_READS_2.fastq

metawrap assembly -1 CLEAN_READS/ALL_READS_1.fastq -2 CLEAN_READS/ALL_READS_2.fastq \
-m 200 -t 96 --metaspades \
-o ASSEMBLY

metawrap binning -a ASSEMBLY/final_assembly.fasta \
-t 96 --metabat2 --maxbin2 --concoct \
-o INITIAL_BINNING \
CLEAN_READS/H*fastq

metawrap bin_refinement -A INITIAL_BINNING/metabat2_bins/ -B INITIAL_BINNING/maxbin2_bins/ -C INITIAL_BINNING/concoct_bins/ \
-t 96 -c 50 -x 10 \
-o BIN_REFINEMENT

metawrap blobology -a ASSEMBLY/final_assembly.fasta \
-t 96 --bins BIN_REFINEMENT/metawrap_50_10_bins \
-o BLOBOLOGY \
CLEAN_READS/H*fastq

metawrap quant_bins -a ASSEMBLY/final_assembly.fasta -b BIN_REFINEMENT/metawrap_50_10_bins \
-o QUANT_BINS \
CLEAN_READS/H*fastq
 

metawrap reassemble_bins -1 CLEAN_READS/ALL_READS_1.fastq -2 CLEAN_READS/ALL_READS_2.fastq -b BIN_REFINEMENT/metawrap_50_10_bins \
-t 96 -m 800 -c 50 -x 10 \
-o BIN_REASSEMBLY

metawrap classify_bins -b BIN_REASSEMBLY/reassembled_bins \
-t 48 \
-o BIN_CLASSIFICATION 


mkdir $resdir/BIN_REASSEMBLY/reassembled_bins_renamedcontig

cd $resdir/BIN_REASSEMBLY
for fa in "$resdir/BIN_REASSEMBLY/reassembled_bins/"*.fa; do

  [ -f "$fa" ] || continue
  origname=$(basename "$fa")
  filename="${origname%.fa}"
  
  $bbmap/rename.sh \
  in=$fa \
  out=reassembled_bins_renamedcontig/$filename.renamedcontig.fa \
  prefix=contig
  
done

cd $resdir
metaWRAP annotate_bins -o FUNCT_ANNOT -t 96 -b BIN_REASSEMBLY/reassembled_bins_renamedcontig/

