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
    read_count_file = name + "_trimmed_sorted_readCount.txt"
    with open(read_count_file) as f:
        read_count = int(f.readline().strip())
        nf = 10 ** 9 / read_count

    if args.pinpoint:
        intersect_list = [name + 'pinpointed_NTS.bed', name + 'pinpointed_TS.bed']
    else:
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

    if args.pinpoint:
        rpkm_list = [name + 'pinpointed_NTS_rpkm.bed', name + 'pinpointed_TS_rpkm.bed']
    else:
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

        out = 'results/' + rpkm_file[:-4] + '_rdg_averages.txt'
        with open(out, 'a') as f:
            f.write("Bin\tAverage Read Count\n")
            for key, value in bins.items():
                f.write(f'{key}\t{mean(value)}\n')


def monomer_analysis(name, desired_lengths=None):
    fasta = name + '.fa'
    if desired_lengths is None:
        desired_lengths = list(range(10, 31))

    print('Performing monomer analysis')

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


def pinpoint(name, lower_boundary=8, upper_boundary=9, dimer_list=None):
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
            uppper_boundary = length - position_interval[0] + 1
            lower_boundary = length - position_interval[1]
            dimer_found = 0
            poses = []
            for pos, nt in enumerate(seq):
                if lower_boundary < pos < uppper_boundary:
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


if version_info[1] < 6:
    raise Exception('Use at least Python 3.6!')

parser = argparse.ArgumentParser(description='Performs the Python half of the XR-Seq analysis.')
parser.add_argument('-s', '--sample_name', help='name of the sample without extensions')
parser.add_argument('-m', '--min_length', default=10, help='Minimum oligomer length for monomer analysis.')
parser.add_argument('-M', '--max_length', default=30, help='Maximum oligomer length for monomer analysis.')
parser.add_argument('-p', '--pinpoint', action='store_true', help='Pinpoint damage sites')
parser.add_argument('-d', '--dimers', help='Dimers to look for while pinpointing. Example: TC,CT Default: TT')
parser.add_argument('-l', '--lower', default=9, type=int, help="Lower boundary for damage location, n bp away from 3' "
                                                               "end. Default: 9")
parser.add_argument('-u', '--upper', default=8, type=int, help="Upper boundary for damage location, n bp away from 3' "
                                                               "end. Default: 8")
args = parser.parse_args()

if args.sample_name:
    samples = [args.sample_name]
else:
    samples = [i[:-19] for i in listdir() if '_trimmed_sorted.bed' in i]

if args.dimers:
    dimers = [dimer.upper() for dimer in args.dimers.split(',')]
else:
    dimers = ["TT"]

mon_lenghts = list(range(args.min_length, args.max_length + 1))

for sample in samples:
    if args.pinpoint:
        pinpoint(sample, args.lower, args.upper, args.dimers)
    rpkm(sample)
    with open("results/" + sample + "_NTS_rpkm.bed") as score_test:
        if type(BedLine(score_test.readline()).score) == int:
            calc_avg(sample)
    monomer_analysis(sample, mon_lenghts)
