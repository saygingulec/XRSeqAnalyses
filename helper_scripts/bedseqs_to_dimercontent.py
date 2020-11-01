"""
Coding strand:          AAATTTAAATTT
Template strand:        TTTAAATTTAAA
Transcribed strand:     TTTAAATTTAAA
Non-transcribed strand: AAATTTAAATTT

We have sequences of the coding strand.
"""

nts_dimers = ['TT', 'TC', 'CT', 'CC']
ts_dimers = ['AA', 'AG', 'GA', 'GG']  # look for AA on coding strand to find TT on transcribed strand


class BedLine:
    def __init__(self, line):
        columns = line.strip().split()
        self.chrom = columns[0]
        self.start = int(columns[1])
        self.end = int(columns[2])
        self.name = columns[3]
        self.score = int(columns[4])
        self.strand = columns[5]
        self.seq = columns[6]

    def bed_line(self):
        """
        This part won't work properly with versions of Python earlier than 3.6 because by default and class attributes
        aren't ordered until Python 3.6.
        Be careful, Longleaf loads Python 3.5.1 by default, but has 3.6.8 if you don't load anything.
        """
        return "\t".join([str(value) for value in vars(self).values()]) + "\n"

    def __len__(self):
        return self.end - self.start


with open('ce11_150bin_seqs.bed') as f,\
        open('ce11_150bin_TS_tpc_content.bed', 'a') as t,\
        open('ce11_150bin_NTS_tpc_content.bed', 'a') as n:

    for line in f:
        line = BedLine(line)
        ts_dimer_count = 0
        nts_dimer_count = 0
        for pos, nt in enumerate(line.seq.upper()):
            try:
                if nt + line.seq[pos + 1] in ts_dimers:
                    ts_dimer_count += 1
                if nt + line.seq[pos + 1] in nts_dimers:
                    nts_dimer_count += 1
            except IndexError:
                pass
        line.seq = ts_dimer_count / len(line)
        t.write(line.bed_line())
        line.seq = nts_dimer_count / len(line)
        n.write(line.bed_line())
