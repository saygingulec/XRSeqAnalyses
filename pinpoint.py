from Bio import SeqIO

fasta_in = 'CelWTL1_6-4_1h_ATCACG_S3_L008_R1_001_trimmed_sorted_fltrd.fa'
bed_out = 'CelWTL1_6-4_1h_ATCACG_S3_L008_R1_001_pinpointed.bed'
dimers = ['TT', 'TC', 'CT']  # these need to be in uppercase
position_interval = [5, 10]  # 5-10 bp away from 3' end

with open(fasta_in) as f:
    with open(bed_out, 'a') as d:
        for record in SeqIO.parse(f, 'fasta'):

            chrom = record.id.split(':')[0]
            start = int(record.id.split('-')[0].split(':')[1])
            end = int(record.id.split('-')[1].split('(')[0])
            strand = record.id[-2]

            seq = record.seq.lower()
            length = len(seq)
            uppper_boundary = length - position_interval[0]
            lower_boundary = length - position_interval[1]
            dimer_found = 0
            poses = []
            for pos, nt in enumerate(seq):
                if lower_boundary < pos < uppper_boundary:
                    dimer = nt + seq[pos-1]
                    if dimer.upper() in dimers:
                        seq = seq[:pos-1] + dimer.upper() + seq[pos + 1:]
                        if strand == '+':
                            poses.append(pos)
                            poses.append(pos - 1)
                        elif strand == '-':
                            poses.append(length - pos)
                            poses.append(length - pos + 1)
                        dimer_found += 1
            if dimer_found > 0:
                end = start + max(poses)
                start += min(poses)
                d.write(f'{chrom}\t{start}\t{end}\t{seq}\t{dimer_found}\t{strand}\n')
