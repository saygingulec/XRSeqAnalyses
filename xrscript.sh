#!/bin/bash

# Starting with compressed samples (sample.fastq.gz), performs XR-Seq analysis.

echo This script performs XR-Seq analysis on compressed samples in the prompted directory.
read -s -p $'Press enter to continue.\n'

# Modules

echo Loading modules

echo Loading Bowtie2
module load bowtie2

echo Loading cutadapt
module load cutadapt

echo Loading samtools
module load samtools

echo Loading bedtools
module load bedtools


# Setting up directories

read -p 'Sample data directory: ' SAMPLE_DIR

read -p 'Bowtie2 index directory: ' BOWTIE2_IND


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


# Unzip without deleting the compressed file

echo Unzipping samples

for sample in ${samples[@]}; do
  echo Unzipping ${sample}.fastq.gz
  gunzip -c "${SAMPLE_DIR}/${sample}.fastq.gz" > "${SAMPLE_DIR}/${sample}.fastq"
done


# Cutting the adapter with cutadapt

echo Cutting adapters
for sample in ${samples[@]}; do
  echo Cutting adapters from ${sample}
  sbatch --job-name=${sample} --wrap="cutadapt -a TGGAATTCTCGGGTGCCAAGGAACTCCAGTNNNNNNACGATCTCGTATGCCGTCTTCTGCTTG -o ${SAMPLE_DIR}/${sample}_cutadapt.fastq ${SAMPLE_DIR}/${sample}.fastq"
done


# Genome allignment

echo Alligning samples with the reference genome
for sample in ${samples[@]}; do
  echo Alligning ${sample}
  sbatch --dependency=singleton --job-name=${sample} --wrap="bowtie2 -p 4 -x ${BOWTIE2_IND} -U ${SAMPLE_DIR}/${sample}_cutadapt.fastq -S ${SAMPLE_DIR}/${sample}_cutadapt.sam"
done


# Convert to BAM

echo Converting allignments to BAM
for sample in ${samples[@]}; do
  echo Converting ${sample} to BAM
  sbatch --dependency=singleton --job-name=${sample} --wrap="samtools view -q 20 -b -o ${SAMPLE_DIR}/${sample}_cutadapt.bam ${SAMPLE_DIR}/${sample}_cutadapt.sam"
done


# Convert to BAD

echo Converting allignments from BAM to BAD
for sample in ${samples[@]}; do
  echo Converting ${sample} to BAD
  sbatch --dependency=singleton --job-name=${sample} --wrap="bedtools bamtobed -i ${SAMPLE_DIR}/${sample}_cutadapt.bam > ${SAMPLE_DIR}/${sample}_cutadapt.bed"
done


# Sorting coordinates by removing duplicates

echo Sorting coordinates and removing duplicates
for sample in ${samples[@]}; do
  echo Sorting ${sample}
  sbatch --dependency=singleton --job-name=${sample} --wrap="sort -u -k1,1 -k2,2n -k3,3n ${SAMPLE_DIR}/${sample}_cutadapt.bed > ${SAMPLE_DIR}/${sample}_cutadapt_sorted.bed"
done 
