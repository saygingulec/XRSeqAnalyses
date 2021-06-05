GENELIST_NAME=hg38_genes
GENE_LONGER_THAN=5000
DISTANCE_GREATER_THAN=5000

bedtools closest -a ${GENELIST_NAME}.bed -b ${GENELIST_NAME}.bed -d -N > ${GENELIST_NAME}_closest.bed
awk '$13 >= DISTANCE_GREATER_THAN' ${GENELIST_NAME}_closest.bed > ${GENELIST_NAME}_far.bed
awk '$3 - $2 > GENE_LONGER_THAN' ${GENELIST_NAME}_far.bed > ${GENELIST_NAME}_far_and_long.bed
cut -f 1-6 ${GENELIST_NAME}_far_and_long.bed > ${GENELIST_NAME}_selected.bed
