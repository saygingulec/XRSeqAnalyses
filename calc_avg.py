from statistics import mean

sample_list = [
    'CelRB1801L16-41h_TTAGGC_S1_L008_R1_001_NTS_rpbm.txt',
    'CelRB1801L16-41h_TTAGGC_S1_L008_R1_001_TS_rpbm.txt',
    'CelRB1801L1CPD1h_CGATGT_S2_L008_R1_001_NTS_rpbm.txt',
    'CelRB1801L1CPD1h_CGATGT_S2_L008_R1_001_TS_rpbm.txt',
    'CelWTL1_6-4_1h_ATCACG_S3_L008_R1_001_NTS_rpbm.txt',
    'CelWTL1_6-4_1h_ATCACG_S3_L008_R1_001_TS_rpbm.txt'
]

for i in sample_list:
    with open(i) as f:

        bins = {}
        for j in range(20):
            bins[j+1] = []

        line_no = 0
        for line in f:
            line_no += 1
            bins[line_no].append(float(line.strip().split()[6]))
            if line_no == 20:
                line_no = 0

    avg_out = i[:-4] + "_rdg_averages.txt"
    with open(avg_out, "a") as f:
        f.write("Bin\tAverage Read Count\n")
        for key, value in bins.items():
            f.write(str(key) + "\t" + str(mean(value)) + "\n")

