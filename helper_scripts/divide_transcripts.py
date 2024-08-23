"""
Divides transcripts into 100 (tx_bin_count) bins and add 2 kb (up_length) of 
up- and downstream sequences divided into 25 (up_bin_count) bins of length 80 bp.

The input and output follows the BED format. The score column is used to denote
the bin index.

input row
chr2L	82420	87387	net_NM_080081.3	0	-

output rows from one transcript
...
chr2L	252360	252373	Rpp30_NM_001272879.1	26	-
chr2L	252346	252360	Rpp30_NM_001272879.1	27	-
chr2L	252333	252346	Rpp30_NM_001272879.1	28	-
...
"""

from statistics import mean

# Parameters
input_gene_list = 'ce11_chosen.bed'
output_divided_gene_list = 'ce11_150bin.bed'
up_length = 2000
up_bin_count = 25
tx_bin_count = 100

up_bin_length = int(up_length / up_bin_count)
total_number_of_bins = tx_bin_count + up_bin_count * 2

with open(input_gene_list) as f:
    for line in f:
        col = line.strip().split()
        chrom = col[0]
        tx_start = int(col[1])
        tx_end = int(col[2])
        tx_name = col[3]
        strand = col[5]
        tx_length = tx_end - tx_start

        # upstream
        up_starts = []
        up_start = tx_start - up_length
        up_starts.append(up_start)
        for i in range(up_bin_count - 1):
            up_starts.append(up_starts[-1] + up_bin_length)
        up_ends = up_starts[1:]
        up_ends.append(tx_start)

        # transcript
        ideal_bin_size = tx_length / tx_bin_count
        bin_size = int(ideal_bin_size)
        tx_starts = [tx_start]
        bin_sizes = [ideal_bin_size]
        for i in range(1, tx_bin_count):
            if mean(bin_sizes) < ideal_bin_size:
                tx_starts.append(tx_starts[-1] + int(bin_size) + 1)
                bin_sizes.append(tx_starts[-1] - tx_starts[-2])
            else:
                tx_starts.append(tx_starts[-1] + int(bin_size))
                bin_sizes.append(tx_starts[-1] - tx_starts[-2])
        tx_ends = tx_starts[1:]
        tx_ends.append(tx_end)

        # downstream
        down_starts = [tx_end + x * up_bin_length for x in range(up_bin_count)]
        down_ends = [tx_end + (x+1) * up_bin_length for x in range(up_bin_count)]

        starts = up_starts + tx_starts + down_starts
        ends = up_ends + tx_ends + down_ends

        if strand == '-':
            starts.reverse()
            ends.reverse()

        with open(output_divided_gene_list, 'a') as d:
            for i in range(total_number_of_bins):
                d.write('{chrom}\t{tx_start}\t{tx_end}\t{tx_name}\t{bin_no}\t{strand}\n'.format(
                    chrom=chrom,
                    tx_start=starts[i],
                    tx_end=ends[i],
                    bin_no=(i+1),
                    tx_name=tx_name,
                    strand=strand,
                ))
