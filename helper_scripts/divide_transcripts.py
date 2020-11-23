from statistics import mean

up_length = 2000
up_bin_count = 25
tx_bin_count = 100

up_bin_length = int(up_length / up_bin_count)
total_number_of_bins = tx_bin_count + up_bin_count * 2

with open('ce11_chosen.bed') as f:
    for line in f:
        col = line.strip().split()
        chrom = col[0]
        tx_start = int(col[1])
        tx_end = int(col[2])
        strand = col[3]
        tx_name = col[5]
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

        with open('ce11_150bin.bed', 'a') as d:
            for i in range(total_number_of_bins):
                d.write('{chrom}\t{tx_start}\t{tx_end}\t{tx_name}\t{bin_no}\t{strand}\n'.format(
                    chrom=chrom,
                    tx_start=starts[i],
                    tx_end=ends[i],
                    bin_no=(i+1),
                    tx_name=tx_name,
                    strand=strand,
                ))
