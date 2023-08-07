# RBPBench
Evaluate CLIP-seq and other genomic region data using a comprehensive collection of known RBP binding motifs. RBPBench can be used for a variety of
purposes, from RBP motif search in genomic regions, over motif co-occurrence analysis, to benchmarking CLIP-seq peak callers.

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

RBPBench was developed and tested on Linux (Ubuntu 22.04 LTS). Currently only Linux operating systems are supported. To install RBPBench, you first need to install the Conda package manager.

### Conda

If you do not have Conda on your system yet, you can e.g. install miniconda, a free + lightweight Conda installer. Get miniconda [here](https://docs.conda.io/en/latest/miniconda.html), choose the latest Miniconda3 Python Linux 64-bit installer and follow the installation instructions. In the end, Conda should be evocable on the command line via (possibly in a different version):

```
$ conda --version
conda 23.7.2
```


### Conda package installation

RBPBench is (soon) available as a Conda package, which makes installation a breeze. 
We simply create a Conda environment and install RBPBench inside the environment:

```
conda create -n rbpbench -c conda-forge -c bioconda
conda activate rbpbench
conda install -c bioconda rbpbench
```

RBPBench should now be available inside the environment:

```
rbpbench -h
```

### Manual installation

Manual installation of RBPBench is only slightly more work. First we create the Conda environment with all necessary dependencies:

```
conda create -n rbpbench -c conda-forge -c bioconda logomaker markdown meme scipy plotly textdistance venn matplotlib-venn infernal bedtools
```

Next we activate the environment, clone the RBPBench repository, and install RBPBench:

```
conda activate rbpbench
git clone https://github.com/michauhl/RBPBench.git
cd RBPBench
python -m pip install . --ignore-installed --no-deps -vv
```

RBPBench should now be available inside the environment:

```
rbpbench -h
```

### RBPBench on Galaxy

RBPBench will soon be available on Galaxy.



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

Let's first check the CLIPper IDR peak regions of the RBP SLBP for motif occurrences. 
For this we first download the peak regions:

```
wget https://www.encodeproject.org/files/ENCFF623WGE/@@download/ENCFF623WGE.bed.gz
gunzip -c ENCFF623WGE.bed.gz | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$7"\t"$6}' > SLBP_K562_IDR_peaks.bed
```

We can use `rbpbench search` to search for RBP motifs by specifying the RBPs we are interested in.
Here we can specify one RBP or several RBPs from the database, as well as all database RBPs, and 
also user-specified ones (see examples further down). 
For a full list of RBPs currently available for search:

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
while the motif hits stats file contains hit statistics for each single motif hit 
(one row per motif hit). 
We can see that out of the 161 input regions, 27 contain a motif hit.


#### Informative statistics

RBPBench produces several informative statistics regarding motif hit (co-)occurrences.
The first statistical test (Wilcoxon rank-sum test) 
**checks whether motif hits are more likely in high-scoring input regions**.
By default BED column 5 is used as region scores (set via `--bed-score-col`), which 
can be e.g. log2 fold changes or -log10 p-values. For example, for the above call

```
rbpbench search --in SLBP_K562_IDR_peaks.bed --rbps SLBP --out SLBP_test_search_out --genome hg38.fa
```

we get a significant enrichment of motif hits (Wilcoxon p-value: 6.87978e-09), 
meaning that input regions with SLBP motif hits also feature significantly higher scores.

The second test addresses **RBP motif co-occurrences** (i.e., between two different RBPs), 
using Fisher's exact test. For this a 2x2 contingency table is constructed between 
each input RBP pair, and the co-occurrence information can also be plotted as an interactive heat map (see Search with multiple RBPs example above). In addition to the co-occurrence statistics, the HTML report also includes a heat map plot of the **correlations between RBPs** (Pearson correlation coefficients). For this genomic input regions are labelled 1 or 0 (RBP motif present or not), resulting in a vector of 1s and 0s for each RBP. Correlations are then calculated by comparing vectors for every pair of RBPs. 


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
to the called peak region (as exemplified in Fig. 3 of [this publication](https://doi.org/10.1016/j.ymeth.2017.02.006)).


#### Search with multiple RBPs

RBPBench can also search for motifs of multiply RBPs in one search run.
For this example we download the PUM1 IDR peak regions:

```
wget https://www.encodeproject.org/files/ENCFF094MQV/@@download/ENCFF094MQV.bed.gz
gunzip -c ENCFF094MQV.bed.gz | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$7"\t"$6}' > PUM1_K562_IDR_peaks.bed
```

When searching with multiple RBPs, motif hit co-occurences and correlations might become of interest. For this we can also output an HTML report file, including hit statistics and interactive plots (`--report`).

```
rbpbench search --in PUM1_K562_IDR_peaks.bed --rbps RBFOX2 PUM1 PUM2 --out test_pum1_out --genome hg38.fa --report
```

The following output files are produced:

```
Co-occurrence p-values for  each RBP pair .tsv:
test_pum1_out/contingency_table_results.tsv
RBP hit stats .tsv:
test_pum1_out/rbp_hit_stats.tsv
Motif hit stats .tsv:
test_pum1_out/motif_hit_stats.tsv
Search report .html:
test_pum1_out/report.rbpbench_search.html
```

The search report `report.rbpbench_search.html` contains a table of motif 
enrichment statistics, as well as the interactive RBP co-occurrences and correlations heat maps.
We can see that PUM1 binding regions with motif hits have significantly higher scores 
(Wilcoxon p-value = 0.00175), and that PUM1 and PUM2 have a significant co-occurrence 
p-value (Fisher's exact test) of 2.48116e-31. In contrast, RBFOX2 co-occurrence p-values
with PUM1 or PUM2 are not significant (0.15234, 0.52303).

Similarly, we can output an HTML report (`motif_plots.rbpbench_search.html`) including the used sequence logos and motif hit statistics (`--plot-motifs`):

```
rbpbench search --in PUM1_K562_IDR_peaks.bed --rbps RBFOX2 PUM1 PUM2 --out test_pum1_out --genome hg38.fa --report --plot-motifs
```

Note that this can take a bit if the whole database (`--rbps ALL`) is used for search (see below).



#### Search with all database RBPs

To use all database RBPs for motif search, we just need to specify `--rbps ALL`:

```
rbpbench search --in SLBP_K562_IDR_peaks.bed --rbps ALL --out SLBP_all_search_out --genome hg38.fa --report
```

This uses all 259 database RBPs (default database: `human_v0.1`, 605 RBP motifs all together). 
Note that SLBP has the lowest Wilcoxon p-value (6.87978e-09) of all 259 RBPs, 
demonstrating that the chosen statistical test is plausible and useful.


#### User-provided motifs

Both sequence (MEME XML format) and structure (covariance model .CM) motifs can be supplied by the user 
on top of the database RBPs (using `-rbps USER` option together with `-user-meme-xml` or `--user-cm`).
For example, let's supply a structure motif (SLBP) via `--user-cm` (the example motif files can be found in the RBPBench repository subfolder `RBPBench/test`):

```
rbpbench search --in SLBP_K562_IDR_peaks.bed --rbps USER --out SLBP_user_search_out --genome hg38.fa --user-cm path_to_test/SLBP_USER.cm  --user-rbp-id SLBP_USER
```

In the same way, we can supply sequence motif(s) (PUM1) via `--user-meme-xml`, and e.g. also combine it with any of the database RBPs (here PUM2 and RBFOX2):

```
rbpbench search --in PUM1_K562_IDR_peaks.bed --rbps USER PUM2 RBFOX2 --out PUM1_user_search_out --genome hg38.fa --user-meme-xml path_to_test/PUM1_USER.xml --user-rbp-id PUM1_USER
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