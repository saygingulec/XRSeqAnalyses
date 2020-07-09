from statistics import mean

number_of_bins = 150
sample_list = [
    "Experiment1_Sample1_trimmed_sorted_NTS.txt",
    "Experiment1_Sample2_trimmed_sorted_TS.txt",
    "Experiment2_Sample1_trimmed_sorted_NTS.txt",
    "Experiment2_Sample2_trimmed_sorted_TS.txt"
]

for i in sample_list:
    with open(i) as f:

        bins = {}
        for j in range(number_of_bins):
            bins[j+1] = []

        line_no = 0
        for line in f:
            line_no += 1
            bins[line_no].append(float(line.split()[7]))
            if line_no == number_of_bins:
                line_no = 0

    avg_out = i[:-4] + "_rdg_averages.txt"
    with open(avg_out, "a") as f:
        f.write("Bin\tAverage Read Count\n")
        for key, value in bins.items():
            f.write(str(key) + "\t" + str(mean(value)) + "\n")
