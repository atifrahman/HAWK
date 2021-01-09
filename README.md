# HAWK 1.7.0

Hitting associations with k-mers

## Modifications in version 1.7.0 

- Implemented IRLS instead of Newton-Raphson for logistic regression, and scaled feature values to fix convergence issues 


## Modifications in versions 1.6.0 and 1.5.0 

- Confounder correction step has been reimplemented in C++ to improve speed and flexibility
- A bug related to order of samples in confounder factor correction has been fixed
- Option for Benjamini-Hochberg procedure for multiple testing correction has been added
- Support for Jellyfish 2 added
- Stats lookup has been reimplemented in C++

The modifications are described in detail in the [paper](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0245058)

## Installation

To install HAWK run (X.Y.Z is the version)

```
tar xzf hawk-X.Y.Z-beta.tar.gz
cd hawk-X.Y.Z-beta
make
```

## Prerequisites

Jellyfish/Jellyfish 2 (modified versions available in supplements)

EIGENSTRAT (modified version available in supplements)

ABYSS 

GNU `sort` with parallel support



If you want to run the old version

R (with foreach and doParallel packages)

## Counting k-mers

The first step in the pipeline is to count k-mers in each sample, find 
total number of k-mers per sample, discard k-mers that appear once in samples and sort
the k-mers. The k-mer file contains one line per k-mer present and each 
line contains an integer representing the k-mer and its count separated 
by a space. The integer representation is given by using 0 for 'A', 
1 for 'C', 2 for 'G' and 3 for 'T'.

k-mer counting can be done using a modified version of the tool Jellyfish/Jellyfish 2
provided in the 'supplements' folder with HAWK. All of the steps mentioned 
above can be performed by installing this version of Jellyfish/Jellyfish 2 and then 
running the script `countKmers_jf1` or `countKmers_jf2` in supplements with necessary modifications for Jellyfish or Jellyfish2 respectively.

Both versions provided assume reads from each sample is in a separate directory 
and prefixes of all directories containing reads is `Reads`. For example reads from
`sample1`, `sample 2`, etc. could be in directories named `Reads_sample1`, `Reads_sample2`, etc.
It also assumes that the read files are gzipped and have extensions `fastq.gz`. If the read files 
are not gzipped please change the `zcat *.fastq.gz` to `cat *.fastq` in `line 21` in `countKmers`. 

This will write the names of sorted k-mer count files in 'sorted_files.txt' 
and total k-mer count in samples in 'total_kmer_counts.txt'.

To use the most recent version of Jellyfish 2, download the source files from this [link](https://github.com/gmarcais/Jellyfish/releases/tag/v2.2.10) and apply the patch dump_jf2.patch. To apply the patch go to Jellyfish 2 
source folder (parent directory of sub_commands folder) and apply

```Bash
patch -p1 < /path/to/dump_jf2.patch
```

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

The k-mers with significant association to cases and controls will be in 
'case_kmers.fasta' and 'control_kmers.fasta' 

## Assembling k-mers with significant association

The k-mers with significant association to cases and controls can be assembled by running

```
./runAbyss
```

The assembled sequences will be in 'case_abyss.25_49.fasta' and 'control_abyss.25_49.fasta'
respectively.

## Changing default confounders

By default HAWK uses first two principal components found using EIGENSTRAT, sex of samples and sequencing depth in the form of total k-mer counts as confounders. To change the number of principal components to be used as confounders, change the value of the variable `noPC` in the script `runHawk` .

The following variables in the script `runHawk` can be used to control confounding factor correction.

```
noPc      The number of principal components from EIGENSTRAT to be used as confounders. 
          Allowed values are from 1 to 10. Default 2.

useSexConfounder  True or False. Set True if sex information 
                   from gwas_info.txt will be used as confounder.
                   Otherwise set False. Default True

covFile   File path to provide additional covariate from external 
          source. Default is empty (no external covariate).
```

Besides, the variable noThread can be used to set number of threads to be used.

```noPc      The number of principal components from EIGENSTRAT to be used as 			  confounders. Allowed values are from 1 to 10. Default 2.
noThread  # of thread used in correction. To speed up the 
          process increase it upto # of core. Default 1.
```

<strong>Note:</strong> To use sex as confounder, information for all sample must be provided. Even if only one sample has no sex information (i.e. no entry or entry other than entry `M/F`), it will make the program avoid using this information.


#### Additional confounder file format
If there are additional confounders, they need to be specified in a separate file. The file <strong>must</strong> be in the following format:  
* Number of lines in the file must be equal to the number of read samples used to run HAWK (i.e. in gwas_info.txt). 
* In each line, the covariates corresponding to a sample must be specified. The covariates must be numbers separated by tabs or spaces.
* All samples (i.e. for all lines) should contain equal number of entries (i.e. number of columns). No line can have missing entries, other wise program will not use them correctly.

For example, if number of samples is 2 and we like to include 3 additional covariates per sample then the file should be as below:
```
0.1    0.2    0.3
0.4    0.5    0.6
```
First three are for the first sample, second three are for the second sample. 

## Using Benjamini-Hochberg procedure for correcting for multiple testing

HAWK uses Bonferroni correction to correct for multiple testing. However, if the study is underpowered for Bonferroni correction, Benjamini-Hochberg procedure may be used to control the false discovery rate (FDR). To use Benjamini-Hochberg procedure, run ```./runBHCorrection``` 

after `runHawk` is executed successfully.  The outputs will be in `case_kmers_bh_correction.fasta` and  `control_kmers_bh_correction.fasta`.

## Finding sex specific k-mers

To find sex specific k-mers using HAWK, in the file 'gwas_info.txt', the first column will contain sample IDs, second column must contain U's and third column will contain Case/Control status depending on whether the sample is Male/Female. For example if SRR3050845 and 
SRR3050847 are female and SRR3050846 is male, the file will be

```
SRR3050845	U	Control
SRR3050846	U	Case
SRR3050847	U	Control
```

Note that the middle column must be all U. Otherwise, the sexes will be used as confounders and k-mers will not be found correctly.


## Obtaining summary stats

Once the k-mers associated with cases and controls are assembled, average p-values, average counts of constituent k-mers and average number of times they are present can be looked up using the script `runKmerSummary` in the `Supplements` folder. Please copy the script to the directory with the output files, edit HAWK directory, input filename, and whether the sequences are from case or control in the script and run 
```./runKmerSummary```
This will output the following six values to the screen separated by commas for each sequence in the input fasta file.
* Average p-value after adjusting for confounders
* Average p-value without adjusting for confounders
* Total count of k-mers in cases averaged over all constituent k-mers
* Total count of k-mers in controls averaged over all constituent k-mers
* Average number of case samples the constituent k-mers are present in
* Average number of control samples the constituent k-mers are present in

## Downstream analysis

After running HAWK and getting sequences associated with cases and controls, they can be analyzed in a number of ways.
* If no reference genome is available, the assembled sequences can be [Blasted](https://blast.ncbi.nlm.nih.gov/Blast.cgi) to check for hits to sequences in related organisms and subsequently analyzed for presence of genes
* If a reference genome is available, the k-mers can be mapped to the reference and their positions and p-values can be visualized using Manhattan plots. Example scripts to run Bowtie 2 and generate plots are provided in the directory ecoli_analysis. The shell script `runBowtie2` can be used to align k-mers to a reference and the R script `manhattan_plasmid.R` can be modified to generate Manhattan plots using the output file in SAM format generated by Bowtie 2. 


## Citations

If you use HAWK, please cite

* Rahman, Atif, Ingileif Hallgrímsdóttir, Michael Eisen, and Lior Pachter. "Association mapping from sequencing reads using k-mers." Elife 7 (2018): e32920.
* Mehrab Z, Mobin J, Tahmid IA, Rahman A (2021) Efficient association mapping from k-mers—An application in finding sex-specific sequences. PLOS ONE 16(1): e0245058. https://doi.org/10.1371/journal.pone.0245058.
