#!/bin/bash

# Monomer analysis includes the reads that were filtered out while pinpointing.

# Command log
echo bash "$0" "$@" > command_log.txt

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
   echo -e "\t--div"
   echo -e "\t\tOptional divided gene list. If you want to get both the TS-NTS log2 comparison table and the genome-wide TS-NTS comparison graph you can use this. Example: ce11_divided.bed"
   echo -e "\t-g | --genome"
   echo -e "\t\tGenome fasta. Example: ce11.fa"
   echo -e "\t-m | --min_length"
   echo -e "\t\tMinimum read length allowed. Example: 13"
   echo -e "\t-M | --max_length"
   echo -e "\t\tMaximum read length allowed. Example: 27"
   echo -e "\t-w | --whole_sample"
   echo -e "\t\tUse all of the reads aligned to the genome for RPKM normalization, read length distribution and monomer analysis."
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
      --div)
          DIV_GENELIST="$2"
          shift 2;;
      -g|--genome)
          GENOME="$2"
          shift 2;;
      -m|--min_length)
          MIN="$2"
          MON_MIN="--mon_min $2"
          shift 2;;
      -M|--max_length)
          MAX="$2"
          MON_MAX="--mon_max $2"
          shift 2;;
      -w|--whole_sample)
          WHOLE="true"
          shift 2;;
      --mon_min)
          MON_MIN="--mon_min $2"
          shift 2;;
      --mon_max)
          MON_MAX="--mon_max $2"
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


# Checks

if [[ -n "$UPPER" ]] && [[ -n "$LOWER" ]]; then
  if [ "${UPPER//[!0-9]/}" -gt "${LOWER//[!0-9]/}" ]; then
      echo "Upper and lower boundaries are the number of nucleotides before the 3' end. The correct arguments would be: --upper ${LOWER//[!0-9]/} --lower ${UPPER//[!0-9]/}"
      exit 1
  fi
fi


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

echo Loading fastqc
module load fastqc

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

echo Loading R 4.0.3
module load r/4.0.3


# Quality control
for SAMPLE in "${SAMPLES[@]}"; do
  sbatch --job-name="${SAMPLE}" --output="slurm-%j-${SAMPLE}-fastqc.out" --wrap="fastqc ${SAMPLE_DIR}/${SAMPLE}.fastq.gz"
done

# Cutting the adapter with Cutadapt

for SAMPLE in "${SAMPLES[@]}"; do
  sbatch --dependency=singleton --job-name="${SAMPLE}" --output="slurm-%j-${SAMPLE}-cutadapt.out" --wrap="cutadapt -a TGGAATTCTCGGGTGCCAAGGAACTCCAGTNNNNNNACGATCTCGTATGCCGTCTTCTGCTTG --discard-untrimmed -m 1 -o ${SAMPLE}_trimmed.fastq ${SAMPLE_DIR}/${SAMPLE}.fastq.gz"
done


# Merging identical reads

for SAMPLE in "${SAMPLES[@]}"; do
  sbatch --mem=8g -t 120 --dependency=singleton --job-name="${SAMPLE}" --output="slurm-%j-${SAMPLE}-fastx.out" --wrap="fastx_collapser -v -i ${SAMPLE}_trimmed.fastq -o ${SAMPLE}_trimmed.fasta -Q33"
done


# Genome alignment

for SAMPLE in "${SAMPLES[@]}"; do
  sbatch --dependency=singleton --job-name="${SAMPLE}" --output="slurm-%j-${SAMPLE}-bowtie2.out" --wrap="bowtie2 -x ${BOWTIE2_IND} -f ${SAMPLE}_trimmed.fasta -S ${SAMPLE}_trimmed.sam --very-sensitive"
done


# Convert to sorted BAM

for SAMPLE in "${SAMPLES[@]}"; do
  sbatch --mem=16g --dependency=singleton --job-name="${SAMPLE}" --output="slurm-%j-${SAMPLE}-samtools_sort.out" --wrap="samtools sort -o ${SAMPLE}_trimmed_sorted.bam ${SAMPLE}_trimmed.sam"
done


# Convert to BED

for SAMPLE in "${SAMPLES[@]}"; do
  sbatch --dependency=singleton --job-name="${SAMPLE}" --output="slurm-%j-${SAMPLE}-bamtobed.out" --wrap="bedtools bamtobed -i ${SAMPLE}_trimmed_sorted.bam > ${SAMPLE}_trimmed_sorted.bed"
done


# Filter by length

if [[ -n "$MIN" ]] && [[ -n "$MAX" ]]; then
  for SAMPLE in "${SAMPLES[@]}"; do
    sbatch --dependency=singleton --job-name="${SAMPLE}" --output="slurm-%j-${SAMPLE}-length_filter.out" --wrap="awk '{if(\$3-\$2>=$MIN && \$3-\$2<=$MAX){print}}' ${SAMPLE}_trimmed_sorted.bed > ${SAMPLE}_filtered.bed"
  done
elif [[ -n "$MIN" ]]; then
  for SAMPLE in "${SAMPLES[@]}"; do
    sbatch --dependency=singleton --job-name="${SAMPLE}" --output="slurm-%j-${SAMPLE}-length_filter.out" --wrap="awk '{if(\$3-\$2>=$MIN){print}}' ${SAMPLE}_trimmed_sorted.bed > ${SAMPLE}_filtered.bed"
  done
elif [[ -n "$MAX" ]]; then
  for SAMPLE in "${SAMPLES[@]}"; do
    sbatch --dependency=singleton --job-name="${SAMPLE}" --output="slurm-%j-${SAMPLE}-length_filter.out" --wrap="awk '{if(\$3-\$2<=$MAX){print}}' ${SAMPLE}_trimmed_sorted.bed > ${SAMPLE}_filtered.bed"
  done
fi


# Get FASTA from BED

for SAMPLE in "${SAMPLES[@]}"; do
  sbatch --dependency=singleton --job-name="${SAMPLE}" --output="slurm-%j-${SAMPLE}-getfasta.out" --wrap="bedtools getfasta -s -fi ${GENOME} -bed ${SAMPLE}_trimmed_sorted.bed -fo ${SAMPLE}.fa"
done

if [[ -n "$MIN" ]] || [[ -n "$MAX" ]]; then
  for SAMPLE in "${SAMPLES[@]}"; do
    sbatch --dependency=singleton --job-name="${SAMPLE}" --output="slurm-%j-${SAMPLE}-getfasta_filtered.out" --wrap="bedtools getfasta -s -fi ${GENOME} -bed ${SAMPLE}_filtered.bed -fo ${SAMPLE}_filtered.fa"
  done
fi


# Pinpoint damage sites

if [ "$PIN" == "-p" ]; then
  for SAMPLE in "${SAMPLES[@]}"; do
    sbatch --dependency=singleton --job-name="${SAMPLE}" --output="slurm-%j-${SAMPLE}-pinpoint.out" --output="slurm-%j-${SAMPLE}-pinpoint.out" --wrap="python3 XR_Seq.py -s ${SAMPLE} ${PIN} ${DIMERS} ${LOWER} ${UPPER}"
  done
fi


# Count total mapped reads

if [ -n "$WHOLE" ]; then
  for SAMPLE in "${SAMPLES[@]}"; do
    sbatch --dependency=singleton --job-name="${SAMPLE}" --output="slurm-%j-${SAMPLE}-read_count.out" --wrap="grep -c \"^\" ${SAMPLE}_trimmed_sorted.bed > ${SAMPLE}_readCount.txt"
    sbatch --dependency=singleton --job-name="${SAMPLE}" --output="slurm-%j-${SAMPLE}-total_mapped_reads.out" --wrap="wc -l ${SAMPLE}_trimmed_sorted.bed >> results/total_mapped_reads.txt"
  done
elif [ "$PIN" == "-p" ]; then
  for SAMPLE in "${SAMPLES[@]}"; do
    sbatch --dependency=singleton --job-name="${SAMPLE}" --output="slurm-%j-${SAMPLE}-read_count.out" --wrap="grep -c \"^\" ${SAMPLE}_pinpointed.bed > ${SAMPLE}_pinpointed_readCount.txt"
    sbatch --dependency=singleton --job-name="${SAMPLE}" --output="slurm-%j-${SAMPLE}-total_pinpointed_reads.out" --wrap="wc -l ${SAMPLE}_pinpointed.bed >> results/total_pinpointed_reads.txt"
    sbatch --dependency=singleton --job-name="${SAMPLE}" --output="slurm-%j-${SAMPLE}-total_mapped_reads.out" --wrap="wc -l ${SAMPLE}_trimmed_sorted.bed >> results/total_mapped_reads.txt"
  done
elif [[ -n "$MIN" ]] || [[ -n "$MAX" ]]; then
  for SAMPLE in "${SAMPLES[@]}"; do
    sbatch --dependency=singleton --job-name="${SAMPLE}" --output="slurm-%j-${SAMPLE}-read_count.out" --wrap="grep -c \"^\" ${SAMPLE}_filtered.bed > ${SAMPLE}_filtered_readCount.txt"
    sbatch --dependency=singleton --job-name="${SAMPLE}" --output="slurm-%j-${SAMPLE}-total_filtered_reads.out" --wrap="wc -l ${SAMPLE}_filtered.bed >> results/total_filtered_reads.txt"
    sbatch --dependency=singleton --job-name="${SAMPLE}" --output="slurm-%j-${SAMPLE}-total_mapped_reads.out" --wrap="wc -l ${SAMPLE}_trimmed_sorted.bed >> results/total_mapped_reads.txt"
  done
else
  for SAMPLE in "${SAMPLES[@]}"; do
    sbatch --dependency=singleton --job-name="${SAMPLE}" --output="slurm-%j-${SAMPLE}-read_count.out" --wrap="grep -c \"^\" ${SAMPLE}_trimmed_sorted.bed > ${SAMPLE}_readCount.txt"
    sbatch --dependency=singleton --job-name="${SAMPLE}" --output="slurm-%j-${SAMPLE}-total_mapped_reads.out" --wrap="wc -l ${SAMPLE}_trimmed_sorted.bed >> results/total_mapped_reads.txt"
  done
fi


# Generate TS and NTS read counts

if [ "$PIN" == "-p" ]; then
  for SAMPLE in "${SAMPLES[@]}"; do
    sbatch --mem=32g --dependency=singleton --job-name="${SAMPLE}" --output="slurm-%j-${SAMPLE}-intersect_NTS.out" --wrap="bedtools intersect -c -a ${GENELIST} -b ${SAMPLE}_pinpointed.bed -wa -s -F 0.5 > ${SAMPLE}_pinpointed_NTS.bed"
    sbatch --mem=32g --dependency=singleton --job-name="${SAMPLE}" --output="slurm-%j-${SAMPLE}-intersect_TS.out" --wrap="bedtools intersect -c -a ${GENELIST} -b ${SAMPLE}_pinpointed.bed -wa -S -F 0.5 > ${SAMPLE}_pinpointed_TS.bed"
  done
elif [[ -n "$MIN" ]] || [[ -n "$MAX" ]]; then
  for SAMPLE in "${SAMPLES[@]}"; do
    sbatch --mem=32g --dependency=singleton --job-name="${SAMPLE}" --output="slurm-%j-${SAMPLE}-intersect_NTS.out" --wrap="bedtools intersect -c -a ${GENELIST} -b ${SAMPLE}_filtered.bed -wa -s -F 0.5 > ${SAMPLE}_filtered_NTS.bed"
    sbatch --mem=32g --dependency=singleton --job-name="${SAMPLE}" --output="slurm-%j-${SAMPLE}-intersect_TS.out" --wrap="bedtools intersect -c -a ${GENELIST} -b ${SAMPLE}_filtered.bed -wa -S -F 0.5 > ${SAMPLE}_filtered_TS.bed"
  done
else
  for SAMPLE in "${SAMPLES[@]}"; do
    sbatch --mem=32g --dependency=singleton --job-name="${SAMPLE}" --output="slurm-%j-${SAMPLE}-intersect_NTS.out" --wrap="bedtools intersect -c -a ${GENELIST} -b ${SAMPLE}_trimmed_sorted.bam -wa -s -F 0.5 > ${SAMPLE}_NTS.bed"
    sbatch --mem=32g --dependency=singleton --job-name="${SAMPLE}" --output="slurm-%j-${SAMPLE}-intersect_TS.out" --wrap="bedtools intersect -c -a ${GENELIST} -b ${SAMPLE}_trimmed_sorted.bam -wa -S -F 0.5 > ${SAMPLE}_TS.bed"
  done
fi

if [ -n "$DIV_GENELIST" ]; then
  if [ "$PIN" == "-p" ]; then
    for SAMPLE in "${SAMPLES[@]}"; do
      sbatch --mem=32g --dependency=singleton --job-name="${SAMPLE}" --output="slurm-%j-${SAMPLE}-div_intersect_NTS.out" --wrap="bedtools intersect -c -a ${DIV_GENELIST} -b ${SAMPLE}_pinpointed.bed -wa -s -F 0.5 > ${SAMPLE}_pinpointed_div_NTS.bed"
      sbatch --mem=32g --dependency=singleton --job-name="${SAMPLE}" --output="slurm-%j-${SAMPLE}-div_intersect_TS.out" --wrap="bedtools intersect -c -a ${DIV_GENELIST} -b ${SAMPLE}_pinpointed.bed -wa -S -F 0.5 > ${SAMPLE}_pinpointed_div_TS.bed"
    done
  elif [[ -n "$MIN" ]] || [[ -n "$MAX" ]]; then
    for SAMPLE in "${SAMPLES[@]}"; do
      sbatch --mem=32g --dependency=singleton --job-name="${SAMPLE}" --output="slurm-%j-${SAMPLE}-div_intersect_NTS.out" --wrap="bedtools intersect -c -a ${DIV_GENELIST} -b ${SAMPLE}_filtered.bed -wa -s -F 0.5 > ${SAMPLE}_filtered_div_NTS.bed"
      sbatch --mem=32g --dependency=singleton --job-name="${SAMPLE}" --output="slurm-%j-${SAMPLE}-div_intersect_TS.out" --wrap="bedtools intersect -c -a ${DIV_GENELIST} -b ${SAMPLE}_filtered.bed -wa -S -F 0.5 > ${SAMPLE}_filtered_div_TS.bed"
    done
  else
    for SAMPLE in "${SAMPLES[@]}"; do
      sbatch --mem=32g --dependency=singleton --job-name="${SAMPLE}" --output="slurm-%j-${SAMPLE}-div_intersect_NTS.out" --wrap="bedtools intersect -c -a ${DIV_GENELIST} -b ${SAMPLE}_trimmed_sorted.bam -wa -s -F 0.5 > ${SAMPLE}_div_NTS.bed"
      sbatch --mem=32g --dependency=singleton --job-name="${SAMPLE}" --output="slurm-%j-${SAMPLE}-div_intersect_TS.out" --wrap="bedtools intersect -c -a ${DIV_GENELIST} -b ${SAMPLE}_trimmed_sorted.bam -wa -S -F 0.5 > ${SAMPLE}_div_TS.bed"
    done
  fi
fi


# Generate BedGraph files

if [ "$PIN" == "-p" ]; then
  for SAMPLE in "${SAMPLES[@]}"; do
    sbatch --dependency=singleton --job-name="${SAMPLE}" --output="slurm-%j-${SAMPLE}-bedgraph_plus.out" --wrap="bedtools genomecov -i ${SAMPLE}_pinpointed.bed -g ${GENOME}.fai -bg -strand + -scale \$(cat ${SAMPLE}_pinpointed_readCount.txt | awk '{print 10000000/\$1}') > ${SAMPLE}_plus.bedGraph"
    sbatch --dependency=singleton --job-name="${SAMPLE}" --output="slurm-%j-${SAMPLE}-bedgraph_minus.out" --wrap="bedtools genomecov -i ${SAMPLE}_pinpointed.bed -g ${GENOME}.fai -bg -strand - -scale \$(cat ${SAMPLE}_pinpointed_readCount.txt | awk '{print 10000000/\$1}') > ${SAMPLE}_minus.bedGraph"
  done
elif [[ -n "$MIN" ]] || [[ -n "$MAX" ]]; then
  for SAMPLE in "${SAMPLES[@]}"; do
    sbatch --dependency=singleton --job-name="${SAMPLE}" --output="slurm-%j-${SAMPLE}-bedgraph_plus.out" --wrap="bedtools genomecov -i ${SAMPLE}_filtered.bed -g ${GENOME}.fai -bg -strand + -scale \$(cat ${SAMPLE}_filtered_readCount.txt | awk '{print 10000000/\$1}') > ${SAMPLE}_plus.bedGraph"
    sbatch --dependency=singleton --job-name="${SAMPLE}" --output="slurm-%j-${SAMPLE}-bedgraph_minus.out" --wrap="bedtools genomecov -i ${SAMPLE}_filtered.bed -g ${GENOME}.fai -bg -strand - -scale \$(cat ${SAMPLE}_filtered_readCount.txt | awk '{print 10000000/\$1}') > ${SAMPLE}_minus.bedGraph"
  done
else
  for SAMPLE in "${SAMPLES[@]}"; do
    sbatch --dependency=singleton --job-name="${SAMPLE}" --output="slurm-%j-${SAMPLE}-bedgraph_plus.out" --wrap="bedtools genomecov -i ${SAMPLE}_trimmed_sorted.bed -g ${GENOME}.fai -bg -strand + -scale \$(cat ${SAMPLE}_readCount.txt | awk '{print 10000000/\$1}') > ${SAMPLE}_plus.bedGraph"
    sbatch --dependency=singleton --job-name="${SAMPLE}" --output="slurm-%j-${SAMPLE}-bedgraph_minus.out" --wrap="bedtools genomecov -i ${SAMPLE}_trimmed_sorted.bed -g ${GENOME}.fai -bg -strand - -scale \$(cat ${SAMPLE}_readCount.txt | awk '{print 10000000/\$1}') > ${SAMPLE}_minus.bedGraph"
  done
fi


# Sort BedGraph files

for SAMPLE in "${SAMPLES[@]}"; do
  sbatch --dependency=singleton --job-name="${SAMPLE}" --output="slurm-%j-${SAMPLE}-bedgraph_plus_sort.out" --wrap="sort -k1,1 -k2,2n ${SAMPLE}_plus.bedGraph > ${SAMPLE}_plus_sorted.bedGraph"
  sbatch --dependency=singleton --job-name="${SAMPLE}" --output="slurm-%j-${SAMPLE}-bedgraph_minus_sort.out" --wrap="sort -k1,1 -k2,2n ${SAMPLE}_minus.bedGraph > ${SAMPLE}_minus_sorted.bedGraph"
done


# Generate BigWig files

for SAMPLE in "${SAMPLES[@]}"; do
  sbatch --dependency=singleton --job-name="${SAMPLE}" --output="slurm-%j-${SAMPLE}-bedGraphToBigWig_plus.out" --wrap="bedGraphToBigWig ${SAMPLE}_plus_sorted.bedGraph ${GENOME}.fai results/${SAMPLE}_plus.bw"
  sbatch --dependency=singleton --job-name="${SAMPLE}" --output="slurm-%j-${SAMPLE}-bedGraphToBigWig_minus.out" --wrap="bedGraphToBigWig ${SAMPLE}_minus_sorted.bedGraph ${GENOME}.fai results/${SAMPLE}_minus.bw"
done


# Generate read length distribution table

if [ "$PIN" == "-p" ]; then
  for SAMPLE in "${SAMPLES[@]}"; do
    sbatch --dependency=singleton --job-name="${SAMPLE}" --output="slurm-%j-${SAMPLE}-lengt_dist_pinpointed.out" --wrap="awk '{print length(\$4)}' ${SAMPLE}_pinpointed.bed | sort -k1,1n | uniq -c | sed 's/\s\s*/ /g' | awk '{print \$2\"\t\"\$1}' > results/${SAMPLE}_read_length_distribution_after_pinpointing.txt"
  done
elif [[ -n "$MIN" ]] || [[ -n "$MAX" ]]; then
  for SAMPLE in "${SAMPLES[@]}"; do
    sbatch --dependency=singleton --job-name="${SAMPLE}" --output="slurm-%j-${SAMPLE}-lengt_dist_filtered.out" --wrap="awk '{print \$3-\$2}' ${SAMPLE}_filtered.bed | sort -k1,1n | uniq -c | sed 's/\s\s*/ /g' | awk '{print \$2\"\t\"\$1}' > results/${SAMPLE}_read_length_distribution_after_filtering.txt"
  done
fi

for SAMPLE in "${SAMPLES[@]}"; do
  sbatch --dependency=singleton --job-name="${SAMPLE}" --output="slurm-%j-${SAMPLE}-lengt_dist.out" --wrap="awk '{print \$3-\$2}' ${SAMPLE}_trimmed_sorted.bed | sort -k1,1n | uniq -c | sed 's/\s\s*/ /g' | awk '{print \$2\"\t\"\$1}' > results/${SAMPLE}_read_length_distribution.txt"
done


# Run Python scripts

for SAMPLE in "${SAMPLES[@]}"; do
  sbatch -t 1-00:00:00 --dependency=singleton --job-name="${SAMPLE}" --output="slurm-%j-${SAMPLE}-python.out" --wrap="python3 XR_Seq.py -s ${SAMPLE} ${MON_MIN} ${MON_MAX}"
done


# Make graphs with R

for SAMPLE in "${SAMPLES[@]}"; do
  sbatch --dependency=singleton --job-name="${SAMPLE}" --output="slurm-%j-${SAMPLE}-R.out" --wrap="Rscript --vanilla XR_Seq.R ${SAMPLE}"
done
