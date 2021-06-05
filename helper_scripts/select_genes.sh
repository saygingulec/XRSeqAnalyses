GENELIST_NAME=hg38_genes
declare -i GENE_LONGER_THAN=5000
declare -i DISTANCE_GREATER_THAN=5000

bedtools closest -a ${GENELIST_NAME}.bed -b ${GENELIST_NAME}.bed -d -N -t first > ${GENELIST_NAME}_closest.bed
awk -v dist="$DISTANCE_GREATER_THAN" '$13 >= dist' ${GENELIST_NAME}_closest.bed > ${GENELIST_NAME}_far.bed 
awk -v gene_length="$GENE_LONGER_THAN" '$3 - $2 > gene_length' ${GENELIST_NAME}_far.bed > ${GENELIST_NAME}_far_and_long.bed
cut -f 1-6 ${GENELIST_NAME}_far_and_long.bed > ${GENELIST_NAME}_selected.bed

