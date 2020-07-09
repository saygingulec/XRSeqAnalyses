"""
Divides transcripts into 100 bins and add 2 kb of up- and downstream sequences divided into 80 bp long bins.

Used the code in the block comments instead if you want to add up- and downstream sequences divided into 20 bins
that are sized similarly to the transcript bins.

input row
chr2L	7528	9484	+	CG11023	NM_001272857.1	1617

output rows from one transcript
...
chr2L	7528	7547	+	CG11023	NM_001272857.1	1617
chr2L	7547	7567	+	CG11023	NM_001272857.1	1617
...
"""

transcript_list = "transcript_list.bed"
output = "divided_transcripts.bed"

with open(transcript_list) as f:
    for line in f:
        line = line.split()
        start = int(line[1])
        upstart = start - 2000
        end = int(line[2])
        downend = end + 2000
        txlen = end - start

        starts = []
        ends = []

        # make 25 bins of 80 bp long upstream
        for i in range(25):
        starts.append(str(upstart + i*80))

        # divide the genes into 100 bins
        prev = start
        for i in range(100):

            if txlen/100 - int(txlen/100) >= 0.5:   # the last bin might be too short if you don't account for this

                # round up for even numbered bins
                if i%2 == 0:
                    starts.append(str(prev))
                    llen = int(txlen/100)
                    prev = prev + llen

                # round down for odd numbered bins
                else:
                    starts.append(str(prev))
                    llen = int(txlen/100) + 1
                    prev = prev + llen

            else:
                starts.append(str(prev))
                llen = int(txlen / 100)
                prev = prev + llen

        """
        # insert 20 bins of upstream sized similarly to the transcript bins
        prev = start - int(txlen / 100)
        for i in range(20):

            if txlen / 100 - int(txlen / 100) >= 0.5:  # the last bin might be too short if you don't account for this

                # round up for even numbered bins
                if i % 2 == 0:
                    starts.insert(0, str(prev))
                    llen = int(txlen / 100)
                    prev = prev - llen

                # round down for odd numbered bins
                else:
                    starts.insert(0, str(prev))
                    llen = int(txlen / 100) + 1
                    prev = prev - llen

            else:
                starts.insert(0, str(prev))
                llen = int(txlen / 100)
                prev = prev - llen

        # insert 20 bins of downstream sized similarly to the transcript bins
        prev = end
        for i in range(20):

            if txlen / 100 - int(txlen / 100) >= 0.5:  # the last bin might be too short if you don't account for this

                # round up for even numbered bins
                if i % 2 == 0:
                    starts.append(str(prev))
                    llen = int(txlen / 100)
                    prev = prev + llen

                # round down for odd numbered bins
                else:
                    starts.append(str(prev))
                    llen = int(txlen / 100) + 1
                    prev = prev + llen

            else:
                starts.append(str(prev))
                llen = int(txlen / 100)
                prev = prev + llen
        """

        # make 25 bins of 80 bp long downstream
        for i in range(25):
            starts.append(str(end + i * 80))

        ends = starts[1:]
        # endend = int(starts[-1]) + int(txlen / 100)
        # ends.append(str(endend)) for similarly sized up-down stream
        ends.append(str(downend))

        with open(output, "a") as d:
            for i in range(len(starts)):
                line[1] = starts[i]
                line[2] = ends[i]
                d.write("\t".join(line) + "\n")
