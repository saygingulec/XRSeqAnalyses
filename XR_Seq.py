from sys import version_info
import argparse
from os import listdir
from statistics import mean
from Bio import SeqIO
import requests
import string

if version_info[1] < 6:
    raise Exception('Use at least Python 3.6!')

parser = argparse.ArgumentParser(description='Performs the Python half of the XR-Seq analysis.')
parser.add_argument('-s', '--sample_name', help='name of the sample without extensions')
args = parser.parse_args()

if args.sample_name:
    samples = [args.sample_name]
else:
    samples = [i[:-19] for i in listdir() if '_trimmed_sorted.bed' in i]


class BedLine:
    def __init__(self, line):
        columns = line.strip().split()
        self.chrom = columns[0]
        self.start = int(columns[1])
        self.end = int(columns[2])
        self.name = columns[3]
        self.score = int(columns[4])
        self.strand = columns[5]
        try:
            self.count = float(columns[6])
        except IndexError:
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


def rpkm(sample_list):
    intersect_list = []
    for name in sample_list:
        intersect_list.append(name + '_NTS.bed')
        intersect_list.append(name + '_TS.bed')

    normalization_factors = {}
    with open('results/total_mapped_reads.txt') as f:
        for line in f:
            normalization_factors[line.strip().split()[1][:-19]] = 10 ** 9 / int(line.strip().split()[0])

    for intersect in intersect_list:
        print("Normalizing " + intersect)
        with open(intersect) as f:
            sample_name = "_".join(intersect.split('_')[:-1])
            nf = normalization_factors[sample_name]
            out = intersect[:-4] + "_rpkm.bed"
            for line in f:
                line = BedLine(line)
                line.count = (line.count * nf) / len(line)
                with open(out, 'a') as d:
                    d.write(line.bed_line())


def calc_avg(sample_list):
    rpkm_list = []
    for name in sample_list:
        rpkm_list.append(name + '_NTS_rpkm.bed')
        rpkm_list.append(name + '_TS_rpkm.bed')

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


def monomer_analysis(sample_list, desired_lengths=None):
    if desired_lengths is None:
        desired_lengths = list(range(10, 31))

    print('Performing monomer analysis')

    fasta_list = []
    for name in sample_list:
        fasta_list.append(name + '.fa')

    for fasta in fasta_list:
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


rpkm(samples)
calc_avg(samples)
monomer_analysis(samples)
