# HAWK
Hitting associations with k-mers

## Installation

To install HAWK run (X.Y.Z is the version)

```
tar xf hawk-X.Y.Z-beta.tar
cd hawk-X.Y.Z-beta
make
```

## Prerequisites

JELLYFISH (modified version available in supplements)

EIGENSTRAT (modified version available in supplements)

R (with foreach and doParallel packages)

ABYSS 

GNU `sort` with parallel support

## Counting k-mers

The first step in the pipeline is to count k-mers in each sample, find 
total number of k-mers per sample, discard k-mers that appear once in samples and sort
the k-mers. The k-mer file contains one line per k-mer present and each 
line contains an integer representing the k-mer and its count separated 
by a space. The integer representation is given by using 0 for 'A', 
1 for 'C', 2 for 'G' and 3 for 'T'.

k-mer counting can be done using a modified version of the tool JELLYFISH
provided in the 'supplements' folder with HAWK. All of the steps mentioned 
above can be performed by installing this version of JELLYFISH and then 
running the script 'countKmers' in supplements with necessary modifications.

The version provided assumes reads from each sample is in a separate directory 
and prefixes of all directories containing reads is `Reads`. For example reads from
`sample1`, `sample 2`, etc. could be in directories named `Reads_sample1`, `Reads_sample2`, etc.
It also assumes that the read files are gzipped and have extensions `fastq.gz`. If the read files 
are not gzipped please change the `zcat *.fastq.gz` to `cat *.fastq` in `line 21` in `countKmers`. 


This will write the names of sorted k-mer count files in 'sorted_files.txt' 
and total k-mer count in samples in 'total_kmer_counts.txt'.

## Running HAWK

Copy 'sorted_files.txt' and 'total_kmer_counts.txt' corresponding to the samples 
into a folder as well as a file named 'gwas_info.txt' containing three columns separated by tabs giving a sample ID, male/female/unknown denoted by M/F/U and Case/Control status of the sample for each sample. For example

```
SRR3050845	U	Control
SRR3050846	U	Case
SRR3050847	U	Control
```

Copy the scripts 'runHawk' and 'runAbyss' into the folder and run

```
./runHawk
```

The k-mers with significant association to case and controls will be in 
'case_kmers.fasta' and 'control_kmers.fasta' which can then be assembled by running

```
./runAbyss
```

The assembled 
sequences will be in 'case_abyss.25_49.fasta' and 'control_abyss.25_49.fasta'
respectively.



