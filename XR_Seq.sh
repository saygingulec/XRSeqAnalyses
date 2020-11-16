#!/bin/bash


# Get options

helpFunction()
{
   echo "Need help with the parameters?"
   echo -e "\t-h | --help"
   echo -e "\t\tDisplay help."
   echo -e "\t-d | --dir"
   echo -e "\t\tDirectory that contains the samples. ./ if in the current directory."
   echo -e "\t-b | --bowtie2index"
   echo -e "\t\tBowtie2Index. Example: Bowtie2Index/WBcel235"
   echo -e "\t-l | --gene_list"
   echo -e "\t\tGene list. Example: ce11_divided.bed"
   echo -e "\t-g | --genome"
   echo -e "\t\tGenome fasta. Example: ce11.fa"
   echo -e "\t-m | --min_length"
   echo -e "\t\tMinimum read length allowed. Example: 13"
   echo -e "\t-M | --max_length"
   echo -e "\t\tMaximum read length allowed. Example: 27"
   echo -e "\t--mon_min"
   echo -e "\t\tMinimum read length for monomer analysis. Default: 10 or, if given, the value of -m"
   echo -e "\t--mon_max"
   echo -e "\t\tMaximum read length for monomer analysis. Default: 30 or, if given, the value of -M"
   echo -e "\t-p | --pinpoint"
   echo -e "\t\tPinpoint damage sites."
   echo -e "\t--dimers"
   echo -e "\t\tDimers to look for while pinpointing. Example: TC,CT Default: TT"
   echo -e "\t--lower"
   echo -e "\t\tLower boundary for damage location, n bp away from 3' end. Default: 9"
   echo -e "\t--upper"
   echo -e "\t\tUpper boundary for damage location, n bp away from 3' end. Default: 8"
   exit 1 # Exit script after printing help
}

MIN="1"
MAX=""
MON_MIN=""
MON_MAX="$MAX"
PIN=""
DIMERS=""
LOWER=""
UPPER=""

while [[ $# -gt 0 ]]; do
  key="$1"

  case $key in
      -h|--help)
          helpFunction;;
      -d|--dir)
          SAMPLE_DIR="$2"
          shift 2;;
      -b|--bowtie2index)
          BOWTIE2_IND="$2"
          shift 2;;
      -l|--gene_list)
          GENELIST="$2"
          shift 2;;
      -g|--genome)
          GENOME="$2"
          shift 2;;
      -m|--min_length)
          MIN="$2"
          MON_MIN="-m $2"
          shift 2;;
      -M|--max_length)
          MAX="-M $2"
          shift 2;;
      --mon_min)
          MON_MIN="-m $2"
          shift 2;;
      --mon_max)
          MON_MAX="-M $2"
          shift 2;;
      -p|--pinpoint)
          PIN="-p"
          shift;;
      --dimers)
          DIMERS="-d $2"
          shift 2;;
      --lower)
          LOWER="-l $2"
          shift 2;;
      --upper)
          UPPER="-u $2"
          shift 2;;
      *)
          echo "$1 is not recognized. -h for help"
          exit 1;;
  esac
done


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

mkdir results


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


# Cutting the adapter with Cutadapt

for SAMPLE in "${SAMPLES[@]}"; do
  sbatch --job-name="${SAMPLE}" --wrap="cutadapt -a TGGAATTCTCGGGTGCCAAGGAACTCCAGTNNNNNNACGATCTCGTATGCCGTCTTCTGCTTG --discard-untrimmed -m ${MIN} ${MAX} -o ${SAMPLE}_trimmed.fastq ${SAMPLE_DIR}/${SAMPLE}.fastq.gz"
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


# Get fasta for monomer analysis

for SAMPLE in "${SAMPLES[@]}"; do
    sbatch --dependency=singleton --job-name="${SAMPLE}" --wrap="bedtools getfasta -s -fi ${GENOME} -bed ${SAMPLE}_trimmed_sorted.bed -fo ${SAMPLE}.fa"
done


# Pinpoint damage sites

if [ $PIN == "-p" ]; then
  for SAMPLE in "${SAMPLES[@]}"; do
    sbatch --dependency=singleton --job-name="${SAMPLE}" --wrap="python3 XR_Seq.py -s ${SAMPLE} ${PIN} ${DIMERS} ${LOWER} ${UPPER}"
  done
fi


# Generate TS and NTS read counts

if [ $PIN == "-p" ]; then
  for SAMPLE in "${SAMPLES[@]}"; do
    sbatch --mem=32g --dependency=singleton --job-name="${SAMPLE}" --wrap="bedtools intersect -c -a ${GENELIST} -b ${SAMPLE}_trimmed_sorted.bam -wa -s -F 0.5 > ${SAMPLE}_NTS.bed"
    sbatch --mem=32g --dependency=singleton --job-name="${SAMPLE}" --wrap="bedtools intersect -c -a ${GENELIST} -b ${SAMPLE}_trimmed_sorted.bam -wa -S -F 0.5 > ${SAMPLE}_TS.bed"
  done
else
  for SAMPLE in "${SAMPLES[@]}"; do
    sbatch --mem=32g --dependency=singleton --job-name="${SAMPLE}" --wrap="bedtools intersect -c -a ${GENELIST} -b ${SAMPLE}_pinpointed.bed -wa -s -F 0.5 > ${SAMPLE}_NTS.bed"
    sbatch --mem=32g --dependency=singleton --job-name="${SAMPLE}" --wrap="bedtools intersect -c -a ${GENELIST} -b ${SAMPLE}_pinpointed.bed -wa -S -F 0.5 > ${SAMPLE}_TS.bed"
  done
fi

# Convert to BED

for SAMPLE in "${SAMPLES[@]}"; do
  sbatch --dependency=singleton --job-name="${SAMPLE}" --wrap="bedtools bamtobed -i ${SAMPLE}_trimmed_sorted.bam > ${SAMPLE}_trimmed_sorted.bed"
done


# Count total mapped reads

for SAMPLE in "${SAMPLES[@]}"; do
  sbatch --dependency=singleton --job-name="${SAMPLE}" --wrap="grep -c \"^\" ${SAMPLE}_trimmed_sorted.bed > ${SAMPLE}_trimmed_sorted_readCount.txt"
  sbatch --dependency=singleton --job-name="${SAMPLE}" --wrap="wc -l ${SAMPLE}_trimmed_sorted.bed >> results/total_mapped_reads.txt"
done


# Generate BedGraph files

if [ $PIN == "-p" ]; then
  for SAMPLE in "${SAMPLES[@]}"; do
    sbatch --dependency=singleton --job-name="${SAMPLE}" --wrap="bedtools genomecov -i ${SAMPLE}_pinpointed.bed -g ${GENOME}.fai -bg -strand + -scale \$(cat ${SAMPLE}_trimmed_sorted_readCount.txt | awk '{print 10000000/\$1}') > ${SAMPLE}_plus.bedGraph"
    sbatch --dependency=singleton --job-name="${SAMPLE}" --wrap="bedtools genomecov -i ${SAMPLE}_pinpointed.bed -g ${GENOME}.fai -bg -strand - -scale \$(cat ${SAMPLE}_trimmed_sorted_readCount.txt | awk '{print 10000000/\$1}') > ${SAMPLE}_minus.bedGraph"
  done
else
  for SAMPLE in "${SAMPLES[@]}"; do
    sbatch --dependency=singleton --job-name="${SAMPLE}" --wrap="bedtools genomecov -i ${SAMPLE}_trimmed_sorted.bed -g ${GENOME}.fai -bg -strand + -scale \$(cat ${SAMPLE}_trimmed_sorted_readCount.txt | awk '{print 10000000/\$1}') > ${SAMPLE}_plus.bedGraph"
    sbatch --dependency=singleton --job-name="${SAMPLE}" --wrap="bedtools genomecov -i ${SAMPLE}_trimmed_sorted.bed -g ${GENOME}.fai -bg -strand - -scale \$(cat ${SAMPLE}_trimmed_sorted_readCount.txt | awk '{print 10000000/\$1}') > ${SAMPLE}_minus.bedGraph"
  done
fi


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
  sbatch --dependency=singleton --job-name="${SAMPLE}" --wrap="python3 XR_Seq.py -s ${SAMPLE} ${MON_MIN} ${MON_MAX}"
done