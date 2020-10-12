#!/bin/bash

# Starting with compressed samples (sample.fastq.gz), performs XR-Seq analysis.

echo This script prepares compressed samples for XR-Seq analysis by creating \
     sorted and filtered bed files.
read -r -s -p $'Press enter to continue.\n'

# Modules

echo Loading modules

echo Loading Bowtie2
module load bowtie2

echo Loading Cutadapt
module load cutadapt

echo Loading fastx
module load fastx_toolkit/0.0.14

echo Loading samtools
module load samtools

echo Loading bedtools
module load bedtools

echo Loading UCSCtools
module load ucsctools/320


# Setting up directories

read -r -p 'Sample data directory: ' SAMPLE_DIR

read -r -p 'Bowtie2 index directory: ' BOWTIE2_IND

read -r -p 'List of genes divided into bins: ' GENELIST

read -r -p 'Genome.fa: ' GENOME

mkdir results


# Listing samples in an array

SAMPLES=()
for i in "${SAMPLE_DIR}"/*.fastq.gz; do
  i=${i##*/}
  SAMPLES+=("${i%.fastq.gz}")
done

echo Your samples are:
for SAMPLE in "${SAMPLES[@]}"; do
  echo "${SAMPLE}"
done

read -r -s -p $'Press enter to continue.\n'


# Cutting the adapter with Cutadapt

for SAMPLE in "${SAMPLES[@]}"; do
  sbatch --job-name="${SAMPLE}" --wrap="cutadapt -a TGGAATTCTCGGGTGCCAAGGAACTCCAGTNNNNNNACGATCTCGTATGCCGTCTTCTGCTTG --discard-untrimmed -m 1 -o ${SAMPLE}_trimmed.fastq ${SAMPLE_DIR}/${SAMPLE}.fastq.gz"
done


# Merging identical reads

for SAMPLE in "${SAMPLES[@]}"; do
  sbatch --dependency=singleton --job-name="${SAMPLE}" --wrap="fastx_collapser -v -i ${SAMPLE}_trimmed.fastq -o ${SAMPLE}_trimmed.fasta -Q33"
done


# Genome alignment

for SAMPLE in "${SAMPLES[@]}"; do
  sbatch --dependency=singleton --job-name="${SAMPLE}" --wrap="bowtie2 -x ${BOWTIE2_IND} -f ${SAMPLE}_trimmed.fasta -S ${SAMPLE}_trimmed.sam"
done


# Convert to sorted BAM

for SAMPLE in "${SAMPLES[@]}"; do
  sbatch --mem=16000 --dependency=singleton --job-name="${SAMPLE}" --wrap="samtools sort -o ${SAMPLE}_trimmed_sorted.bam ${SAMPLE}_trimmed.sam"
done


# Generate TS and NTS read counts

for SAMPLE in "${SAMPLES[@]}"; do
  sbatch --dependency=singleton --job-name="${SAMPLE}" --wrap="bedtools intersect -c -a ${GENELIST} -b ${SAMPLE}_trimmed_sorted.bam -wa -s -F 0.5 > ${SAMPLE}_NTS.bed"
  sbatch --dependency=singleton --job-name="${SAMPLE}" --wrap="bedtools intersect -c -a ${GENELIST} -b ${SAMPLE}_trimmed_sorted.bam -wa -S -F 0.5 > ${SAMPLE}_TS.bed"
done


# Convert to BAD

for SAMPLE in "${SAMPLES[@]}"; do
  sbatch --dependency=singleton --job-name="${SAMPLE}" --wrap="bedtools bamtobed -i ${SAMPLE}_trimmed_sorted.bam > ${SAMPLE}_trimmed_sorted.bed"
done


# Count total mapped reads

for SAMPLE in "${SAMPLES[@]}"; do
  sbatch --dependency=singleton --job-name="${SAMPLE}" --wrap="grep -c \"^\" ${SAMPLE}_trimmed_sorted.bed > ${SAMPLE}_trimmed_sorted_readCount.txt"
  sbatch --dependency=singleton --job-name="${SAMPLE}" --wrap="wc -l ${SAMPLE}_trimmed_sorted.bed >> results/total_mapped_reads.txt"
done


# Generate BedGraph files

for SAMPLE in "${SAMPLES[@]}"; do
  sbatch --dependency=singleton --job-name="${SAMPLE}" --wrap="bedtools genomecov -i ${SAMPLE}_trimmed_sorted.bed -g ${GENOME}.fai -bg -strand + -scale \$(cat ${SAMPLE}_trimmed_sorted_readCount.txt | awk '{print 10000000/\$1}') > ${SAMPLE}_trimmed_sorted_plus.bedGraph"
  sbatch --dependency=singleton --job-name="${SAMPLE}" --wrap="bedtools genomecov -i ${SAMPLE}_trimmed_sorted.bed -g ${GENOME}.fai -bg -strand - -scale \$(cat ${SAMPLE}_trimmed_sorted_readCount.txt | awk '{print 10000000/\$1}') > ${SAMPLE}_trimmed_sorted_minus.bedGraph"
done


# Generate BigWig files

for SAMPLE in "${SAMPLES[@]}"; do
  sbatch --dependency=singleton --job-name="${SAMPLE}" --wrap="bedGraphToBigWig ${SAMPLE}_trimmed_sorted_plus.bedGraph ${GENOME}.fai results/${SAMPLE}_trimmed_sorted_plus.bw"
  sbatch --dependency=singleton --job-name="${SAMPLE}" --wrap="bedGraphToBigWig ${SAMPLE}_trimmed_sorted_minus.bedGraph ${GENOME}.fai results/${SAMPLE}_trimmed_sorted_minus.bw"
done


# Generate read length distribution table

for SAMPLE in "${SAMPLES[@]}"; do
  sbatch --dependency=singleton --job-name="${SAMPLE}" --wrap="awk '{print \$3-\$2}' ${SAMPLE}_trimmed_sorted.bed | sort -k1,1n | uniq -c | sed 's/\s\s*/ /g' | awk '{print \$2\"\t\"\$1}' > results/${SAMPLE}_read_length_distribution.txt"
done


# Get fasta for monomer analysis

for SAMPLE in "${SAMPLES[@]}"; do
  sbatch --dependency=singleton --job-name="${SAMPLE}" --wrap="bedtools getfasta -s -fi ${GENOME} -bed ${SAMPLE}_trimmed_sorted.bed -fo ${SAMPLE}.fa"
done


# Run Python scripts

for SAMPLE in "${SAMPLES[@]}"; do
  sbatch --dependency=singleton --job-name="${SAMPLE}" --wrap="python3 XR_Seq.py -s ${SAMPLE}"
done