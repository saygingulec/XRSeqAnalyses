#!/bin/bash

# Starting with compressed samples (sample.fastq.gz), performs XR-Seq analysis.

echo This script prepares compressed samples for XR-Seq analysis by creating \
     sorted and filtered bed files.
read -s -p $'Press enter to continue.\n'

# Modules

echo Loading modules

echo Loading Bowtie2
module load bowtie2

echo Loading bbmap
module load bbmap

echo Loading fastx
module load fastx_toolkit/0.0.14

echo Loading samtools
module load samtools

echo Loading bedtools
module load bedtools


# Setting up directories

read -p 'Sample data directory: ' SAMPLE_DIR

read -p 'Bowtie2 index directory: ' BOWTIE2_IND

read -p 'Adaptor reference file: ' ADAPTOR_REF

# Listing samples in an array

samples=()
for i in ${SAMPLE_DIR}/*.fastq.gz; do
  i=${i##*/}
  samples+=("${i%.fastq.gz}")
done

echo Your samples are:
for sample in ${samples[@]}; do
  echo ${sample}
done

read -s -p $'Press enter to continue.\n'


# Cutting the adapter with bbduk

echo Cutting adapters
for sample in ${samples[@]}; do
  echo Cutting adapters from ${sample}
  sbatch --mem=32000 --job-name=${sample} --wrap="bbduk.sh in=${SAMPLE_DIR}/${sample}.fastq.gz out=${SAMPLE_DIR}/${sample}_trimmed.fastq ktrim=r k=21 hdist=2 minlen=1 mink=15 ref=${ADAPTOR_REF}"
done


# Merging identical reads

echo Merging identical reads
for sample in ${samples[@]}; do
  echo Merging reads in ${sample}
  sbatch --dependency=singleton --job-name=${sample} --wrap="fastx_collapser -v -i ${SAMPLE_DIR}/${sample}_trimmed.fastq -o ${SAMPLE_DIR}/${sample}_trimmed.fasta -Q33"
done


# Genome allignment

echo Alligning samples with the reference genome
for sample in ${samples[@]}; do
  echo Alligning ${sample}
  sbatch --dependency=singleton --job-name=${sample} --wrap="bowtie2 -x ${BOWTIE2_IND} -f ${SAMPLE_DIR}/${sample}_trimmed.fasta -S ${SAMPLE_DIR}/${sample}_trimmed.sam"
done


# Convert to sorted BAM

echo Sorting and converting to BAM
for sample in ${samples[@]}; do
  echo Sorting and converting ${sample} to BAM
  sbatch --mem=16000 --dependency=singleton --job-name=${sample} --wrap="samtools sort -o ${SAMPLE_DIR}/${sample}_trimmed_sorted.bam ${SAMPLE_DIR}/${sample}_trimmed.sam"
done


# Convert to BAD

echo Converting allignments from BAM to BAD
for sample in ${samples[@]}; do
  echo Converting ${sample} to BAD
  sbatch --dependency=singleton --job-name=${sample} --wrap="bedtools bamtobed -i ${SAMPLE_DIR}/${sample}_trimmed_sorted.bam > ${SAMPLE_DIR}/${sample}_trimmed_sorted.bed"
done


# Filtering out reads with less than 10 and more than 32 nt

echo Filtering out reads with less than 10 and more than 32 nt
for sample in ${samples[@]}; do
  echo Filtering ${sample}
  sbatch --dependency=singleton --job-name=${sample} --wrap="awk '{if($3-$2>=10 && $3-$2<=32){print}}' ${SAMPLE_DIR}/${sample}_trimmed_sorted.bed > ${sample}_trimmed_sorted_fltrd.bed"
done
