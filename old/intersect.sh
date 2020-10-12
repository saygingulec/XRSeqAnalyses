#!/bin/bash

# Maps reads to the transcripts divided into bins strandedly.
# Have the BAM files and the gene list in the current directory.

echo Loading bedtools
module load bedtools

read -p 'Gene list: ' genelist

samples=()
for i in ./*.bam; do
  i=${i##*/}
  samples+=("${i}")
done

echo Your samples are:
for sample in ${samples[@]}; do
  echo ${sample}
done

read -s -p $'Press enter to continue.\n'

for sample in ${samples[@]}; do
  echo Intersecting ${sample}
  bedtools intersect -c -a ${genelist} -b ${sample} -wa -S -F 0.5 > ${sample%.bam}_TS.txt
  echo TS done
  bedtools intersect -c -a ${genelist} -b ${sample} -wa -s -F 0.5 > ${sample%.bam}_NTS.txt
  echo NTS done
done
