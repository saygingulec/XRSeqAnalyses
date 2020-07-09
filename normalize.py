"""
Normalizes the read counts using RPBM.
"""

sample_list = [
    "Experiment1_Sample1_trimmed_sorted_NTS.txt",
    "Experiment1_Sample2_trimmed_sorted_TS.txt",
    "Experiment2_Sample1_trimmed_sorted_NTS.txt",
    "Experiment2_Sample2_trimmed_sorted_TS.txt"
]

read_counts = [1000000]

normalization_factors = []
for i in read_counts:
    normalization_factors.append(1000000 / i)
    normalization_factors.append(1000000 / i)

for i in range(len(sample_list)):
    print("Normalizing " + sample_list[i])
    with open(sample_list[i]) as f:
        normalized = sample_list[i][:-4] + "_normalized.txt"
        for line in f:
            line = line.strip().split()
            read_count = int(line[7])
            bin_start = int(line[1])
            bin_end = int(line[2])
            bin_length = bin_end - bin_start
            line[7] = str((read_count * normalization_factors[i]) / bin_length)
            line = "\t".join(line) + "\n"
            with open(normalized, "a") as d:
                d.write(line)
