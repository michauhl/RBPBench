# RBPBench
Evaluate CLIP-seq and other genomic region data using a comprehensive collection of known RBP binding motifs. RBPBench can be used for a variety of
purposes, from RBP motif search in genomic regions, over motif co-occurence analysis, to benchmarking CLIP-seq peak callers.

## Table of contents

- [Introduction](#introduction)
- [Installation](#installation)
    - [Conda](#conda)
    - [Conda package installation](#conda-package-installation)
    - [Manual installation](#manual-installation)
- [Example runs](#example-runs)

## Introduction

CLIP-seq is by far the most common wet-lab protocol to identify the RNA binding sites of an RNA-binding protein (RBP) on a transcriptome-wide scale. Various popular protocol variants exist, such as PAR-CLIP, iCLIP, or eCLIP. In order to identify the binding sites from the mapped CLIP-seq read data, tools termed peak callers are applied on the read data. Regions with enriched read numbers in the CLIP-seq libraries (i.e., peak regions), typically compared to one or more control libraries are subsequently defined as RBP binding sites. Various peak callers exist, applying various techniques in order to define the binding sites.

In order to quantify the performance of a peak caller, the enrichment of known RBP binding motifs in the binding site (or peak region) data can be used. However, there exists no automated tool for this yet. rbp-bench was implemented to fill this gap, providing easy installation, an easy-to-use interface, as well as a comprehensive motif database. rbp-bench is simply installed via conda, and a Galaxy wrapper is also available. Comprehensive statistics and informative plots allow for easy comparisons across multiple RBPs and peak callers.


## Installation

Installation


```
conda create -n rbpbench -c conda-forge -c bioconda
conda activate rbpbench
conda install -c bioconda bedtools
conda install -c bioconda meme
conda install -c bioconda infernal
```

For Kolmogorov-Smirnov test (scipy stats ks_2samp)
```
conda install -c conda-forge scipy
```



### Conda

Conda

### Conda package installation

Conda package installation

### Manual installation

Manual installation

## Example runs

In order to run the examples, we first need to download the human genome sequence:

```
wget https://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/hg38.fa.gz
gunzip hg38.fa.gz
```

We also need some genomic regions to search for motifs. In the examples we mainly use the 
genomic peak regions identified by the eCLIP peak caller CLIPper (IDR method), 
as they are easily available from ENCODE. E.g. for RBP PUM1 (cell type K562):

```
wget https://www.encodeproject.org/files/ENCFF094MQV/@@download/ENCFF094MQV.bed.gz
gunzip -c ENCFF094MQV.bed.gz | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$7"\t"$6}' > PUM1_K562_IDR_peaks.bed
```

Here we reformatted the peak regions BED, to have 6 columns (6-column BED format), 
with column 5 containing the log2 fold change scores of the peak regions.


### Searching for RBP motifs


#### Simple search for single RBPs

Let's first check the CLIPper IDR peak regions of the RBP SLBP for motif occurences. 
For this we first download the peak regions:

```
wget https://www.encodeproject.org/files/ENCFF623WGE/@@download/ENCFF623WGE.bed.gz
gunzip -c ENCFF623WGE.bed.gz | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$7"\t"$6}' > SLBP_K562_IDR_peaks.bed
```

We can use `rbpbench search` to search for RBP motifs by specifying the RBPs we are interested in.
Here we can specify one RBP or several RBPs from the database, as well as all database RBPs, and 
also user-specified ones (see examples further down). 
For a full list of available RBPs for search:

```
rbpbench info
```

Let's first look for SLBP motifs (there is currently one structure motif in the database for SLBP):

```
rbpbench search --in SLBP_K562_IDR_peaks.bed --rbps SLBP --out SLBP_test_search_out --genome hg38.fa
```

Apart from the informative command line output, the search results are stored in two table files:

```
RBP hit stats .tsv:
SLBP_test_search_out/rbp_hit_stats.tsv
Motif hit stats .tsv:
SLBP_test_search_out/motif_hit_stats.tsv
```

The first file (RBP hit stats) contains comprehensive hit statistics for each RBP 
(one row per RBP, see manual below for column descriptions), 
while  the motif hits stats file contains hit statistics for each single motif hit 
(one row per motif hit). 
We can see that out of the 161 input regions, 27 contain a motif hit.


#### Additional search options

There are a number of search options available (check out `rbpbench search -h` and manual below for full descriptions). 
For example, we can extend the given genomic regions, and see how this affects the hit statistics:

```
rbpbench search --in SLBP_K562_IDR_peaks.bed --rbps SLBP --out test_slbp_out --genome hg38.fa --report --ext 10
```

This extends the genomic regions up- and downstream by 10 nt, and results in 51 hits (from 27 hits without extension). 
Moreover, we can do uneven extension to check whether up- or downstream extension has different effects:

```
rbpbench search --in SLBP_K562_IDR_peaks.bed --rbps SLBP --out test_slbp_out --genome hg38.fa --report --ext 20,0
rbpbench search --in SLBP_K562_IDR_peaks.bed --rbps SLBP --out test_slbp_out --genome hg38.fa --report --ext 0,20
```

We note that upstream extension (first command) by 20 nt results in 40 hits, 
while downstream extension (second command) results in 54 hits.
It thus seems that SLBP structure motifs tend to more often reside downstream relative 
to the called peak region (as exemplified in [Uhl et al. 2017](https://doi.org/10.1016/j.ymeth.2017.02.006) Fig. 3).


#### Search with multiple RBPs

```

```



#### Search with all database RBPs

To use all database RBPs for motif search, we just need to specify `--rbps ALL`:

```
rbpbench search --in SLBP_K562_IDR_peaks.bed --rbps ALL --out SLBP_all_search_out --genome hg38.fa --report
```


#### User-provided motifs

Both sequence (MEME XML format) and structure (covariance model .CM) motifs can be supplied by the user 
on top of the database RBPs (using `-rbps USER` option together with `-user-meme-xml` or `--user-cm`).
For example, let's supply a structure motif (SLBP) via `--user-cm`:

```
rbpbench search --in SLBP_K562_IDR_peaks.bed --rbps USER --out SLBP_user_search_out --genome hg38.fa --user-cm test/SLBP_USER.cm  --user-rbp-id SLBP_USER
```

In the same way, we can supply sequence motif(s) (PUM2) via `--user-meme-xml`, and e.g. also combine it with any of the database RBPs (here PUM1 and RBFOX2):

```
rbpbench search --in PUM2_K562_IDR_peaks.bed --rbps USER PUM1 RBFOX2 --out PUM2_user_search_out --genome hg38.fa --user-meme-xml test/PUM2_USER.xml --user-rbp-id PUM2_USER
```


### Batch-processing multiple datasets


```

```




### Comparisons between search results (Benchmarking)


```

```




### Additional functions

#### Plot nucleotide distribution around genomic positions

To plot the nucleotide distribution around genomic positions, we can use 
`rbpbench dist`:

```

```