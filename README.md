# XRSeqAnalyses

These two scripts allow you to get data for monomer analysis, read length distribution, bigWig and TS-NTS comparison graphs. Each sample is run in parallel to save time.

To do the analysis download XR-Seq.sh and XR-Seq.py. Then you need the Bowtie2Index for alignment, a whole genome fasta for monomer analysis, your samples and your choice of intervals for read counts to be reported.

Example run:

```
This script prepares compressed samples for XR-Seq analysis by creating sorted and filtered bed files.
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
Sample data directory: ./
Bowtie2 index directory: Bowtie2Index/WBcel235
List of genes divided into bins: ce11_divided.bed
Genome.fa: Caenorhabditis_elegans.WBcel235.dna.toplevel.fa
Your samples are:
CelWTL1_6-4_1h_ATCACG_S3_L008_R1_001
Press enter to continue.
Submitted batch job 5110002
Submitted batch job 5110003
Submitted batch job 5110005
Submitted batch job 5110006
Submitted batch job 5110008
Submitted batch job 5110009
Submitted batch job 5110011
Submitted batch job 5110012
Submitted batch job 5110013
Submitted batch job 5110015
Submitted batch job 5110016
Submitted batch job 5110018
Submitted batch job 5110019
Submitted batch job 5110020
Submitted batch job 5110022
Submitted batch job 5110023
```
