# XRSeqAnalyses

These two scripts allow you to get data for monomer analysis, read length distribution, bigWig and TS-NTS comparison graphs, and pinpoint the damage sites. Each sample is run concurrently and indepent of one another.

To do the analysis download XR_Seq.sh and XR_Seq.py. You will also need the locations of the following: Bowtie2Index for alignment, whole genome fasta for monomer analysis, your samples and your genelist.

Help:

```
Need help with the parameters?
        -h | --help
                Display help.
        -d | --dir
                Directory that contains the samples. ./ if in the current directory.
        -b | --bowtie2index
                Bowtie2Index. Example: Bowtie2Index/WBcel235
        -l | --gene_list
                Gene list. Example: ce11_divided.bed
        -g | --genome
                Genome fasta. Example: ce11.fa
        -m | --min_length
                Minimum read length allowed. Example: 13
        -M | --max_length
                Maximum read length allowed. Example: 27
        --mon_min
                Minimum read length for monomer analysis. Default: 10 or, if given, the value of -m
        --mon_max
                Maximum read length for monomer analysis. Default: 30 or, if given, the value of -M
        -p | --pinpoint
                Pinpoint damage sites.
        --dimers
                Dimers to look for while pinpointing. Example: TC,CT Default: TT
        --lower
                Lower boundary for damage location, n bp away from 3' end. Default: 9
        --upper
                Upper boundary for damage location, n bp away from 3' end. Default: 8
```

Example run:

```
[sayging@longleaf-login1 test]$ bash XR_Seq.sh -d ./ -b bowtie2Index/Mycobacterium_smegmatis_str_mc2_155 -l smeg_XRseq.bed -g Mycobacterium_smegmatis_str_mc2_155.ASM1500v1.dna.chromosome.Chromosome.fa -m 13 -M 13 -p
Your samples are:
MsmegUvrD5mnRep2_CGTACG_S11_L007_R1_001
MsmegWT5minRep2_GTTTCG_S12_L007_R1_001
Press enter to continue.
Loading modules
Loading Bowtie2
     Genome indexes for use with BOWTIE2 are available
       in /proj/seq/data .
     bwa/bowtie/bowtie2 indexes are located in "Sequence" directory
       under the top level genome directory with name in CAPS, e.g. MM9_UCSC
Loading Cutadapt
Loading fastx
Loading samtools
Loading bedtools
Loading UCSCtools
Loading Python 3.6.6
Submitted batch job 7328659
Submitted batch job 7328660
Submitted batch job 7328661
Submitted batch job 7328662
Submitted batch job 7328663
Submitted batch job 7328664
Submitted batch job 7328665
Submitted batch job 7328666
Submitted batch job 7328667
Submitted batch job 7328668
Submitted batch job 7328669
Submitted batch job 7328670
Submitted batch job 7328671
Submitted batch job 7328672
Submitted batch job 7328673
Submitted batch job 7328674
Submitted batch job 7328675
Submitted batch job 7328676
Submitted batch job 7328678
Submitted batch job 7328679
Submitted batch job 7328680
Submitted batch job 7328681
Submitted batch job 7328682
Submitted batch job 7328683
Submitted batch job 7328684
Submitted batch job 7328685
Submitted batch job 7328686
Submitted batch job 7328687
Submitted batch job 7328689
Submitted batch job 7328690
Submitted batch job 7328691
Submitted batch job 7328692
Submitted batch job 7328693
Submitted batch job 7328694
Submitted batch job 7328695
Submitted batch job 7328696

```
