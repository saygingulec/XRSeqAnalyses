from math import log2
from sys import version_info
import argparse
from os import listdir
from statistics import mean
from Bio import SeqIO
import requests
import string


class BedLine:
    def __init__(self, line):
        columns = line.strip().split()
        self.chrom = columns[0]
        self.start = int(columns[1])
        self.end = int(columns[2])
        self.name = columns[3]
        try:
            self.score = int(columns[4])
        except ValueError:
            self.score = columns[4]
        self.strand = columns[5]
        try:
            self.count = float(columns[-1])
        except ValueError:
            pass

    def bed_line(self):
        """
        This part won't work properly with versions of Python earlier than 3.6 because by default and class attributes
        aren't ordered until Python 3.6.
        Be careful, Longleaf loads Python 3.5.1 by default, but has 3.6.8 if you don't load anything.
        """
        return "\t".join([str(value) for value in vars(self).values()]) + "\n"

    def sequence(self, assembly):

        query = {
            'genome': assembly,
            'chrom': self.chrom,  # might need to change chrom names
            'start': self.start,
            'end': self.end
        }

        # UCSC gives you the sequence on the + strand, so you need to keep it in mind for genes on - strand.
        r = requests.get('https://api.genome.ucsc.edu/getData/sequence?', params=query)
        seq = r.json()['dna'].upper()

        if self.strand == '+':
            return seq

        return seq.translate(string.maketrans("ACTG", "TGAC"))

    def __len__(self):
        return self.end - self.start


def rpkm(name):
    read_count_file = name + "_readCount.txt"
    with open(read_count_file) as f:
        read_count = int(f.readline().strip())
        nf = 10 ** 9 / read_count

    intersect_list = [name + '_NTS.bed', name + '_TS.bed']

    for intersect in intersect_list:
        print("Normalizing " + intersect)
        with open(intersect) as f:
            out = "results/" + intersect[:-4] + "_rpkm.bed"
            for line in f:
                line = BedLine(line)
                line.count = (line.count * nf) / len(line)
                with open(out, 'a') as d:
                    d.write(line.bed_line())


def calc_avg(name):

    rpkm_list = [name + '_NTS_rpkm.bed', name + '_TS_rpkm.bed']

    for rpkm_file in rpkm_list:
        print('Averaging ' + rpkm_file)
        with open(rpkm_file) as f:
            bins = {}
            for line in f:
                line = BedLine(line)
                # Score column is used to show which bin it is.
                if line.score not in bins:
                    bins[line.score] = []
                bins[line.score].append(line.count)

        out = rpkm_file[:-4] + '_rdg_averages.txt'
        with open(out, 'a') as f:
            f.write("Bin\tAverage Read Count\n")
            for key, value in bins.items():
                f.write(f'{key}\t{mean(value)}\n')


def monomer_analysis(name, desired_lengths=None):
    fasta = name + '.fa'
    if desired_lengths is None:
        desired_lengths = list(range(10, 31))

    dict_of_monomers = {
        length: {
            i + 1: {
                'A': 0,
                'T': 0,
                'G': 0,
                'C': 0
            } for i in range(length)
        } for length in desired_lengths
    }

    # Count the number of monomers
    for entry in SeqIO.parse(fasta, 'fasta'):
        seq = entry.seq.upper()
        if 'N' in seq:
            continue
        length = len(seq)
        if length in desired_lengths:
            for pos, nt in enumerate(seq):
                pos += 1  # position is 1 based, not 0 based
                dict_of_monomers[length][pos][nt] += 1

    # Convert numbers to ratios
    for length in dict_of_monomers:
        for pos in dict_of_monomers[length]:
            count_sum = sum(dict_of_monomers[length][pos].values())
            for nt, count in dict_of_monomers[length][pos].items():
                try:
                    dict_of_monomers[length][pos][nt] = count / count_sum
                except ZeroDivisionError:
                    dict_of_monomers[length][pos][nt] = 0

    # Convert to R data frame and write
    r_df_name = 'results/' + fasta[:-3] + '_monomer_R_df.txt'
    with open(r_df_name, 'a') as f:
        for length in dict_of_monomers:
            for pos in dict_of_monomers[length]:
                for nt, ratio in dict_of_monomers[length][pos].items():
                    f.write(f'{length}\t{pos}\t{nt}\t{ratio}\n')

    # Convert to human readable table and write
    hr_table = 'results/' + fasta[:-3] + '_monomer_hr.txt'
    with open(hr_table, 'a') as f:
        f.write('\tA\tT\tG\tC\n1\t')
        for length in dict_of_monomers:
            for pos in dict_of_monomers[length]:
                for nt, ratio in dict_of_monomers[length][pos].items():
                    if pos == length:
                        if nt == 'C':
                            if length != desired_lengths[-1]:
                                f.write(f'{ratio}\n\n\tA\tT\tG\tC\n1\t')
                            else:
                                f.write(f'{ratio}\n')
                        else:
                            f.write(f'{ratio}\t')
                    else:
                        if nt == 'C':
                            f.write(f'{ratio}\n{pos+1}\t')
                        else:
                            f.write(f'{ratio}\t')


def pinpoint(name, lower_boundary=9, upper_boundary=8, dimer_list=None):
    if name + '_filtered.fa' in listdir():
        fasta = name + '_filtered.fa'
    else:
        fasta = name + '.fa'
    out = name + '_pinpointed.bed'
    position_interval = [lower_boundary, upper_boundary]  # 8-9 bp away from 3' end, inclusive
    if dimer_list is None:
        if "6-4" in fasta:
            dimer_list = ['TT', 'TC', 'CT', 'CC']  # these need to be in uppercase
        else:
            dimer_list = ['TT']  # these need to be in uppercase

    with open(fasta) as f, open(out, 'a') as d:
        for record in SeqIO.parse(f, 'fasta'):

            chrom = record.id.split(':')[0]
            start = int(record.id.split('-')[0].split(':')[1])
            # end = int(record.id.split('-')[1].split('(')[0])
            strand = record.id[-2]

            seq = record.seq.lower()
            length = len(seq)
            lower_boundary = length - position_interval[0]
            upper_boundary = length - position_interval[1] + 1
            dimer_found = 0
            poses = []
            for pos, nt in enumerate(seq):
                if lower_boundary < pos <= upper_boundary:
                    dimer = seq[pos - 1] + nt
                    if dimer.upper() in dimer_list:
                        seq = seq[:pos - 1] + dimer.upper() + seq[pos + 1:]
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


def gene_body_total_reads(name):
    out = name + '_TS_NTS_rpkm.txt'
    nts_bins = name + '_NTS_rpkm.bed'
    ts_bins = name + '_TS_rpkm.bed'
    genes = {}
    with open(nts_bins) as f:
        for line in f:
            line = BedLine(line)
            genes[line.name] = {}
            genes[line.name]['NTS'] = float(line.count)
    with open(ts_bins) as f:
        for line in f:
            line = BedLine(line)
            genes[line.name]['TS'] = float(line.count)
    for gene in genes:
        ts = genes[gene]['TS']
        nts = genes[gene]['NTS']
        if ts == 0:
            ts = 1
        if nts == 0:
            nts = 1
        ratio = ts / nts
        lg2 = log2(ratio)
        genes[gene]['ratio'] = ratio
        genes[gene]['log2'] = lg2
    with open(out, 'a') as f:
        f.write('Gene\tTS\tNTS\tRatio\tlog2\n')
        for gene, numbers in sorted(genes.items(), key=lambda x: x[1]['log2'], reverse=True):
            f.write(f'{gene}\t{numbers["TS"]}\t{numbers["NTS"]}\t{numbers["ratio"]}\t{numbers["log2"]}\n')


if version_info[1] < 6:
    raise Exception('Use at least Python 3.6!')

parser = argparse.ArgumentParser(description='Performs the Python half of the XR-Seq analysis.')
parser.add_argument('-s', '--sample_name', help='Name of the sample without extensions.')
parser.add_argument('-m', '--min_length', type=int, help='Minimum allowed read length.')
parser.add_argument('-M', '--max_length', type=int, help='Maximum allowed read length.')
parser.add_argument('--mon_min', type=int, default=10, help='Minimum oligomer length for monomer analysis.')
parser.add_argument('--mon_max', type=int, default=30, help='Maximum oligomer length for monomer analysis.')
parser.add_argument('-p', '--pinpoint', action='store_true', help='Pinpoint damage sites')
parser.add_argument('-d', '--dimers', help='Dimers to look for while pinpointing. Example: TC,CT Default: TT')
parser.add_argument('-l', '--lower', default=9, type=int, help="Lower boundary for damage location, n bp away from 3' "
                                                               "end. Default: 9")
parser.add_argument('-u', '--upper', default=8, type=int, help="Upper boundary for damage location, n bp away from 3' "
                                                               "end. Default: 8")
parser.add_argument('--ma', action='store_true', help='Perform monomer analysis.')
args = parser.parse_args()

# Collect sample names
if args.sample_name:
    samples = [args.sample_name]
else:
    samples = [i[:-19] for i in listdir() if '_trimmed_sorted.bed' in i]

# Change defaults
if args.dimers:
    dimers = [dimer.upper() for dimer in args.dimers.split(',')]
else:
    dimers = ["TT"]

mon_lenghts = list(range(args.mon_min, args.mon_max + 1))

# Only perform monomer analysis
if args.ma:
    for sample in samples:
        print("Performing monomer analysis on " + sample)
        monomer_analysis(sample, mon_lenghts)
        if sample + '_pinpointed' in listdir():
            monomer_analysis(sample + '_pinpointed', mon_lenghts)
        elif sample + '_filtered' in listdir():
            monomer_analysis(sample + '_filtered', mon_lenghts)
    quit()

# Pinpoint
if args.pinpoint:
    for sample in samples:
        print("Pinpointing " + sample)
        pinpoint(sample, args.lower, args.upper, args.dimers)
    quit()

# Calculate RPKM and averages for TCR graph, and perform monomer analysis
for sample in samples:
    if str(sample + "_pinpointed.bed") in listdir():
        print("Normalizing " + sample)
        rpkm(sample + "_pinpointed")
        with open("results/" + sample + "_pinpointed_NTS_rpkm.bed") as score_test:
            if BedLine(score_test.readline()).score > 0:
                print("Calculating bin averages of " + sample)
                calc_avg("results/" + sample + "_pinpointed")
            else:
                print("Comparing TS-NTS")
                gene_body_total_reads("results/" + sample + "_pinpointed")
    elif str(sample + "_filtered.bed") in listdir():
        print("Normalizing " + sample)
        rpkm(sample + "_filtered")
        with open("results/" + sample + "_filtered_NTS_rpkm.bed") as score_test:
            if BedLine(score_test.readline()).score > 0:
                print("Calculating bin averages of " + sample)
                calc_avg("results/" + sample + "_filtered")
            else:
                print("Comparing TS-NTS")
                gene_body_total_reads("results/" + sample + "_filtered")
    else:
        print("Normalizing " + sample)
        rpkm(sample)
        with open("results/" + sample + "_NTS_rpkm.bed") as score_test:
            if BedLine(score_test.readline()).score > 0:
                print("Calculating bin averages of " + sample)
                calc_avg("results/" + sample)
            else:
                print("Comparing TS-NTS")
                gene_body_total_reads("results/" + sample)

    print("Performing monomer analysis on " + sample)
    monomer_analysis(sample, mon_lenghts)
    if sample + '_pinpointed' in listdir():
        monomer_analysis(sample + '_pinpointed', mon_lenghts)
    elif sample + '_filtered' in listdir():
        monomer_analysis(sample + '_filtered', mon_lenghts)
