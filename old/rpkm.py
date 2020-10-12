sample_list = [
    'CelRB1801L16-41h_TTAGGC_S1_L008_R1_001_NTS.txt',
    'CelRB1801L16-41h_TTAGGC_S1_L008_R1_001_TS.txt',
    'CelRB1801L1CPD1h_CGATGT_S2_L008_R1_001_NTS.txt',
    'CelRB1801L1CPD1h_CGATGT_S2_L008_R1_001_TS.txt',
    'CelWTL1_6-4_1h_ATCACG_S3_L008_R1_001_NTS.txt',
    'CelWTL1_6-4_1h_ATCACG_S3_L008_R1_001_TS.txt'
]

read_counts = [24572379, 22517652, 25996586]

normalization_factors = []
for i in read_counts:
    normalization_factors.append(1000000000 / i)
    normalization_factors.append(1000000000 / i)

for i in range(len(sample_list)):
    print("Normalizing " + sample_list[i])
    with open(sample_list[i]) as f:
        normalized = sample_list[i][:-4] + "_rpkm.txt"
        for line in f:
            line = line.strip().split()
            tx_start = int(line[1])
            tx_end = int(line[2])
            tx_length = tx_end - tx_start
            if tx_length < 0:
                print(line)
            read_count = float(line[6])
            line[6] = str((read_count * normalization_factors[i]) / tx_length)
            line = "\t".join(line) + "\n"
            with open(normalized, "a") as d:
                d.write(line)

