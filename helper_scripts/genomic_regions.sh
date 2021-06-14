PREFIX=dm6_ucsc
ASSEMBLY=dm6
GENOME=dm6.chrom.sizes

module load bedtools

# Intron
bedtools subtract -a ${PREFIX}_introns.bed -b ${PREFIX}_exons.bed -s > ${ASSEMBLY}_introns.bed

# 5' UTR
bedtools subtract -a ${PREFIX}_5utr.bed -b ${PREFIX}_coding.bed -s > ${ASSEMBLY}_5utr.bed

# 3' UTR
bedtools subtract -a ${PREFIX}_3utr.bed -b ${PREFIX}_coding.bed -s > ${ASSEMBLY}_3utr-coding.bed
bedtools subtract -a ${ASSEMBLY}_3utr-coding.bed -b ${ASSEMBLY}_5utr.bed -s > ${ASSEMBLY}_3utr.bed

# Core Promoter
cat ${ASSEMBLY}_introns.bed ${ASSEMBLY}_5utr.bed ${PREFIX}_coding.bed ${ASSEMBLY}_3utr.bed > ${ASSEMBLY}_gb.bed
bedtools sort -i ${ASSEMBLY}_gb.bed > ${ASSEMBLY}_genebody.bed
bedtools subtract -a ${PREFIX}_cp.bed -b ${ASSEMBLY}_genebody.bed -s > ${ASSEMBLY}_cp.bed

# Distal Promoter
bedtools subtract -a ${PREFIX}_dp.bed -b ${ASSEMBLY}_genebody.bed -s > ${ASSEMBLY}_predp.bed
bedtools subtract -a ${ASSEMBLY}_predp.bed -b ${ASSEMBLY}_cp.bed -s > ${ASSEMBLY}_dp.bed

# 1kb Downstream
cat ${ASSEMBLY}_genebody.bed ${ASSEMBLY}_cp.bed ${ASSEMBLY}_dp.bed > ${ASSEMBLY}_gp.bed
bedtools sort -i ${ASSEMBLY}_gp.bed > ${ASSEMBLY}_genebody_prom.bed
bedtools subtract -a ${PREFIX}_1kbdown.bed -b ${ASSEMBLY}_genebody_prom.bed -s > ${ASSEMBLY}_1kbdown.bed

# Intergenic
cat ${ASSEMBLY}_genebody_prom.bed ${ASSEMBLY}_1kbdown.bed > ${ASSEMBLY}_g.bed
bedtools sort -i ${ASSEMBLY}_g.bed > ${ASSEMBLY}_genic.bed
bedtools complement -L -i ${ASSEMBLY}_genic.bed -g ${GENOME} > ${ASSEMBLY}_intergenic.bed

# Cleanup
rm ${ASSEMBLY}_3utr-coding.bed
rm ${ASSEMBLY}_predp.bed
rm ${ASSEMBLY}_genebody_prom.bed
rm ${ASSEMBLY}_genic.bed
rm ${ASSEMBLY}_gb.bed
rm ${ASSEMBLY}_gp.bed
rm ${ASSEMBLY}_g.bed

# Lengths

for i in dm6_*; do
  printf "%s" "${i}\t" >> region_lengths.txt
  awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}' "${i}" >> region_lengths.txt
done