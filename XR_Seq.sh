#!/bin/bash

# Modules

echo Loading modules

echo Loading Bowtie2
module load bowtie2

echo Loading Cutadapt
module load cutadapt/2.9

echo Loading fastx
module load fastx_toolkit/0.0.14

echo Loading samtools
module load samtools

echo Loading bedtools
module load bedtools

echo Loading UCSCtools
module load ucsctools/320

echo Loading Python 3.6.6
module load python/3.6.6


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
  sbatch --mem=16g --dependency=singleton --job-name="${SAMPLE}" --wrap="samtools sort -o ${SAMPLE}_trimmed_sorted.bam ${SAMPLE}_trimmed.sam"
done


# Generate TS and NTS read counts

for SAMPLE in "${SAMPLES[@]}"; do
  sbatch --mem=32g --dependency=singleton --job-name="${SAMPLE}" --wrap="bedtools intersect -c -a ${GENELIST} -b ${SAMPLE}_trimmed_sorted.bam -wa -s -F 0.5 > ${SAMPLE}_NTS.bed"
  sbatch --mem=32g --dependency=singleton --job-name="${SAMPLE}" --wrap="bedtools intersect -c -a ${GENELIST} -b ${SAMPLE}_trimmed_sorted.bam -wa -S -F 0.5 > ${SAMPLE}_TS.bed"
done


# Convert to BED

for SAMPLE in "${SAMPLES[@]}"; do
  sbatch --dependency=singleton --job-name="${SAMPLE}" --wrap="bedtools bamtobed -i ${SAMPLE}_trimmed_sorted.bam > ${SAMPLE}_trimmed_sorted.bed"
done


# Count total mapped reads

for SAMPLE in "${SAMPLES[@]}"; do
  sbatch --dependency=singleton --job-name="${SAMPLE}" --wrap="grep -c \"^\" ${SAMPLE}_trimmed_sorted.bed > ${SAMPLE}_trimmed_sorted_readCount.txt"
  sbatch --dependency=singleton --job-name="${SAMPLE}" --wrap="wc -l ${SAMPLE}_trimmed_sorted.bed >> results/total_mapped_reads.txt"
done


# Get fasta for monomer analysis

for SAMPLE in "${SAMPLES[@]}"; do
  sbatch --dependency=singleton --job-name="${SAMPLE}" --wrap="bedtools getfasta -s -fi ${GENOME} -bed ${SAMPLE}_trimmed_sorted.bed -fo ${SAMPLE}.fa"
done


# Generate BedGraph files

for SAMPLE in "${SAMPLES[@]}"; do
  sbatch --dependency=singleton --job-name="${SAMPLE}" --wrap="bedtools genomecov -i ${SAMPLE}_trimmed_sorted.bed -g ${GENOME}.fai -bg -strand + -scale \$(cat ${SAMPLE}_trimmed_sorted_readCount.txt | awk '{print 10000000/\$1}') > ${SAMPLE}_plus.bedGraph"
  sbatch --dependency=singleton --job-name="${SAMPLE}" --wrap="bedtools genomecov -i ${SAMPLE}_trimmed_sorted.bed -g ${GENOME}.fai -bg -strand - -scale \$(cat ${SAMPLE}_trimmed_sorted_readCount.txt | awk '{print 10000000/\$1}') > ${SAMPLE}_minus.bedGraph"
done


# Sort BedGraph files

for SAMPLE in "${SAMPLES[@]}"; do
  sbatch --dependency=singleton --job-name="${SAMPLE}" --wrap="sort -k1,1 -k2,2n ${SAMPLE}_plus.bedGraph > ${SAMPLE}_plus_sorted.bedGraph"
  sbatch --dependency=singleton --job-name="${SAMPLE}" --wrap="sort -k1,1 -k2,2n ${SAMPLE}_minus.bedGraph > ${SAMPLE}_minus_sorted.bedGraph"
done


# Generate BigWig files

for SAMPLE in "${SAMPLES[@]}"; do
  sbatch --dependency=singleton --job-name="${SAMPLE}" --wrap="bedGraphToBigWig ${SAMPLE}_plus_sorted.bedGraph ${GENOME}.fai results/${SAMPLE}_plus.bw"
  sbatch --dependency=singleton --job-name="${SAMPLE}" --wrap="bedGraphToBigWig ${SAMPLE}_minus_sorted.bedGraph ${GENOME}.fai results/${SAMPLE}_minus.bw"
done


# Generate read length distribution table

for SAMPLE in "${SAMPLES[@]}"; do
  sbatch --dependency=singleton --job-name="${SAMPLE}" --wrap="awk '{print \$3-\$2}' ${SAMPLE}_trimmed_sorted.bed | sort -k1,1n | uniq -c | sed 's/\s\s*/ /g' | awk '{print \$2\"\t\"\$1}' > results/${SAMPLE}_read_length_distribution.txt"
done


# Run Python scripts

for SAMPLE in "${SAMPLES[@]}"; do
  sbatch --dependency=singleton --job-name="${SAMPLE}" --wrap="python3 XR_Seq.py -s ${SAMPLE}"
done
