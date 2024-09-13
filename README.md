# RBPBench

[![GitHub](https://img.shields.io/github/tag/michauhl/RBPBench.svg)](https://github.com/michauhl/RBPBench)
[![Bioconda](https://anaconda.org/bioconda/rbpbench/badges/version.svg)](https://anaconda.org/bioconda/rbpbench)

RBPBench is multi-function tool to evaluate CLIP-seq and other genomic region data 
using a comprehensive collection of known RNA-binding protein (RBP) binding motifs. RBPBench can be used for a variety of
purposes, from RBP motif search (database or user-supplied RBPs) in genomic regions, 
over motif enrichment and co-occurrence analysis, 
to benchmarking CLIP-seq peak caller methods as well as comparisons across cell types and CLIP-seq protocols.

## Table of contents

- [Installation](#installation)
    - [Conda](#conda)
    - [Conda package installation](#conda-package-installation)
    - [Manual installation](#manual-installation)
    - [RBPBench webserver](#rbpbench-webserver)
- [Example runs](#example-runs)
- [Documentation](#documentation)

## Installation

RBPBench was developed and tested on Linux (Ubuntu 22.04 LTS). Currently only Linux operating systems are supported. To install RBPBench, you first need to install the Conda package manager.

### Conda

If you do not have Conda on your system yet, you can e.g. install Miniconda, a free + lightweight Conda installer. Get Miniconda [here](https://docs.conda.io/en/latest/miniconda.html), choose the latest Miniconda3 Python Linux 64-bit installer and follow the installation instructions. In the end, Conda should be evocable on the command line via (possibly in a different version):

```
$ conda --version
conda 24.5.0
```


### Conda package installation

RBPBench is available as a Conda package, which makes installing a breeze. 
We simply create a Conda environment and install RBPBench inside the environment:

```
conda create -n rbpbench -c conda-forge -c bioconda rbpbench
conda activate rbpbench
```

RBPBench should now be available inside the environment:

```
rbpbench -h
```

### Manual installation

Manual installation of RBPBench is only slightly more work. First we create the Conda environment with all necessary dependencies:

```
conda create -n rbpbench -c conda-forge -c bioconda logomaker markdown meme scipy plotly textdistance venn matplotlib-venn infernal bedtools upsetplot scikit-learn goatools python=3.11
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

### RBPBench webserver

RBPBench is also available as a webserver on Galaxy (more infos soon).


## Example runs

More information will be added soon.


## Documentation

This documentation provides further details on RBPBench (version 1.0).

### Program modes

RBPBench is a multi-function tool featuring several modes of operation. 
To get an overview of the currently available modes:

```
$ rbpbench -h
usage: rbpbench [-h] [-v]
                {search,batch,compare,searchseq,searchregex,searchlong,searchrna,searchlongrna,enmo,nemo,streme,tomtom,goa,optex,dist,info}
                ...

Evaluate CLIP-seq and other genomic region data using a comprehensive
collection of known RBP binding motifs (RNA sequence + structure). RBPBench
can be used for a variety of purposes, from RBP motif search in genomic
regions, over motif enrichment and co-occurrence analysis, to benchmarking
CLIP-seq peak callers, as well as comparisons across cell types and CLIP-seq
protocols.

positional arguments:
  {search,batch,compare,searchseq,searchregex,searchlong,searchrna,searchlongrna,enmo,nemo,streme,tomtom,goa,optex,dist,info}
                        Program modes
    search              Search motifs in genomic sites
    batch               Search motifs on > 1 dataset
    compare             Compare different search results
    searchseq           Search motifs in sequences
    searchregex         Search regex in genomic sites or sequences
    searchlong          Search motifs in long genomic regions
    searchrna           Search motifs in spliced transcript sites
    searchlongrna       Search motifs in spliced full transcripts
    enmo                Check for enriched motifs in input sites
    nemo                Check for neighboring motifs in input sites
    streme              Discover motifs in input sites using STREME
    tomtom              Compare motif(s) with database using TOMOTM
    goa                 Run GO enrichment analysis on gene list
    optex               Investigate optimal extension
    dist                Plot nt distribution at genomic positions
    info                Inform about motif database content

options:
  -h, --help            show this help message and exit
  -v, --version         show program's version number and exit

```

##### Single and batch motif search

`rbpbench search` can be used to search for RBP motifs in genomic regions, selecting 
any number of RBPs (database or user-defined (sequence or structure), or regular expressions) 
and to look at motif enrichment and co-occurrences.
`rbpbench batch` is `rbpbench search` extended to multiple input files (BED files containing 
genomic regions), with one RBP to search for each input file. Informative HTML reports are produced 
in most of the modes, e.g. for `rbpbench batch` the input datasets are compared via 
various interactive plots and statistics to learn more about common as well as dataset-specific 
features. 

##### Compare search results

`rbpbench search` and `rbpbench batch` produce the RPB and motif hit statistics files which can 
then be input into `rbpbench compare` in order to benchmark peak callers, or generate comparisons 
between different cell types and CLIP-seq protocols. Again an HTML report is produced with plots 
for each possible comparison found in the input files.


##### Additional motif search modes

Additional more specialised search modes are available as well:
`rbpbench searchseq` allows inputting FASTA sequences to search for motifs.
`rbpbench searchregex` allows regular expression (regex) search in genomic regions or sequences. 
The regular expression can be a simple sequence motif string like `AATAAA`, 
or complex ones like `[CA]CA[CT].{10,25}CGGAC`. Note that regex search is also possible in 
the other search modes, e.g. combined with the database sequence and structure motifs.
`rbpbench searchlong` allows searching long genomic regions for RBP motifs and visualize 
motif prevalences in genomic regions (e.g. whether a motif tends to occur more often in introns, 
5'UTRs, 3'UTRs, etc.).
`rbpbench searchrna` enables RBP motif search in transcript regions, i.e., in spliced transcripts, 
which can be useful for studying e.g. RBP binding on mRNAs.
`rbpbench searchlongrna` similarly allows for RBP motif search in spliced transcripts, but broadens 
the search to all transcripts (or a user-defined set) of the GTF provided transriptome (optionally only on mRNAs),
and e.g. outputs mRNA region preferences of the specified motifs (5'UTR, CDS, 3'UTR coverages).

##### Single motif enrichment search modes

The next two modes, `rbpbench enmo` and `rbpbench nemo`, further zoom in on single motif level 
analysis of motif enrichment and co-occurrence. 
While `rbpbench enmo` looks at enriched motifs over background sequences in general, 
`rbpbench nemo` goes one step further and searches for significantly enriched neighboring motifs. 
I.e., given some genomic or transcript positions, `rbpbench nemo` checks whether there are enriched 
motifs up- or downstream of the input positions, and also tests for statistical preference regarding
up- or downstream occurrence.

##### Checking for similar database motifs

If you have one or more sequence motifs of interest (or even a regex!) and want to know if there are similar motifs 
in the database, you can run `rbpbench tomtom` (including MEME SUITE TOMTOM). 
This mode also informs us whether the reported motif hits are enriched in certain RBP functions (e.g., 
for regex input `--in TTTTTT`, the top 3 enriched functions are: splicing regulation, 
translation regulation, and RNA stability & decay). 

##### More modes

```rbpbench streme``` allows us to discover new motifs using STREME from MEME suite. 
```rbpbench goa``` enables us to run GO term enrichment analysis (GOA) on a set of genes, e.g. obtained 
from the output tables of other modes (like genes covered by input regions from ```rbpbench search```
or ```rbpbench batch```). Note that some of the modes also can run GOA, optimized for the specific mode
(see [GO term analysis](#go-term-analysis)).
Furthermore, we can use `rbpbench dist` to plot the nucleotide distribution at genomic positions, 
and ```rbpbench info``` informs about database RBP motifs and annotated RBP functions.

##### Connecting modes

All motif search modes additionally output motif hits in BED format. This way, the motif hit regions themselves 
can again be used as input regions for the other modes, allowing for easily refined hypothesis checking.


### Inputs


#### Motifs

Motifs for search can be sequence motifs ([MEME motif format](https://meme-suite.org/meme/doc/meme-format.html)), 
structure motifs (covariance models obtained through [Infernal](https://github.com/EddyRivasLab/infernal)),
or regular expressions (e.g., a simple sequence motif string like `AATAAA` or more complex ones like `[CA]CA[CT].{10,25}CGGAC`).
Search motifs are selected via `--rbps`, e.g., to select PUM1 (sequence) and SLBP (structure) motifs `--rbps PUM1 SLBP`.
To select all database motifs, set `--rbps ALL` (internal motif database can be changed via `--motif-db`, 
a custom motif database can be supplied too). Additionally, a regular expression (regex) (e.g. `AATAAA`) can be added to 
search by `--regex AATAAA`. To search only for a regular expression, set `--rbps REGEX --regex AATTA`. Co-occurrence 
and enrichment statistics are calculated on the RBP level in all modes, except `rbpbench enmo` and `rbpbench nemo`, 
which enable enrichment and co-occurrence statistics on single motif level. 
Alternatively, single motifs for search can be selected via `--motifs`.
To list and visualize all selected motifs, use `--plot-motifs` (in some modes automatically output).
RBP motifs can also be filtered by their annotated RBP functions. To only inlcude RBPs with e.g. annotated 3' end processing
function, set `--functions TEP` (for available functions and annotations, run `rbpbench info`).
Moreover, if you have a motif of interest (or a regex) and want to know if there are similar motifs in the database, you can 
use `rbpbench tomtom` (using MEME SUITE TOMTOM). This mode also informs us whether the reported motif hits 
are enriched in certain RBP functions (e.g. for regex input `--in TTTTTT`, the top 3 enriched functions are:
splicing regulation, translation regulation, and RNA stability & decay). 


#### Genomic regions

Genomic input regions have to be provided in BED format. The first 6 columns 
are mandatory and are expected to mean the following: chromosome ID, 
genomic region start position (0-based index),
genomic region end position (1-based index), 
region ID, region score, region strand (plus(+) or minus(-) strand).
Here's a valid input example from the `test` folder:

```
$ head -5 eclip_clipper_idr/PUM1_K562_IDR_peaks.bed 
chr20	62139082	62139128	PUM1_K562_IDR	3.43241573178832	-
chr20	62139146	62139197	PUM1_K562_IDR	3.35874858317445	-
chr7	6156853	6157005	PUM1_K562_IDR	4.85347590189745	+
chr15	82404676	82404753	PUM1_K562_IDR	4.1721908622051	+
chr1	19094999	19095025	PUM1_K562_IDR	5.11052671530143	-
```

Additional columns (> 6) will be ignored, although the region score column 
can be different from column 5 (default). You can specify which column should be used 
as region score via `--bed-score-col` (used for Wilcoxon rank sum test and optionally for 
filtering out regions by `--bed-sc-thr`).
Before motif search, the input regions are filtered and optionally extended via `--ext` 
(e.g. `--ext 20` for up- and downstream extension of 20 nt or `--ext 30,10` for different up- 
and downstream extension).
Regions with invalid chromosome IDs are removed. Furthermore, duplicated regions are merged 
into one region (same chromosome ID + start + end + strand).


#### Transcript regions

Transcript region input is also supported (`rbpbench searchrna`, `rbpbench enmo`, `rbpbench nemo`).
Here the input regions look the following:
```
ENST00000645463	100	150	s1	0	+
ENST00000645463	300	350	s2	0	+
ENST00000561978	450	500	s3	0	+
```

The defined transcript IDs need to be present in the supplied GTF file (`--gtf`).
This allows motif search and subsequent statistics (enrichment, co-occurrence, GO term analysis etc.) 
directly on the transcript regions, which allows us to examine motif binding directly on the 
spliced transcripts.

#### Genomic sequences

Genomic sequences have to be provided in FASTA format. For example, the human 
genome sequence (`hg38.fa`) can be obtained from [here](https://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips).
To download in the command line and unpack:

```
wget https://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/hg38.fa.gz
gunzip hg38.fa.gz
```

### Handling genomic overlaps

As input genomic regions can overlap (also due to extending them via `--ext`), RBPBench 
considers both the called genomic region size as well as the effective genomic 
region size. While the called size ignores overlaps, the effective size 
only counts unique input regions (i.e., merging overlapping parts, only counting
them once). This also holds for the reported motif hits, where **unique** counts 
refer to effective counts. For example, a motif at a specific genomic location 
can appear several times in the input data. The unique count takes care of this 
and counts it only once. This is important e.g. when comparing the results 
of different peak callers. Another interesting statistic (full list of statistics and descriptions [here](#hit-statistics-table-files)) 
is e.g. the number of unique motif hits over 1000 nt of called and effective region size. 
This gives us an idea of how many motifs are included in the regions, normalized over 
the total size of the regions (called or effective size).



### Outputs

Output descriptions




### Statistics

Wilcoxon rank sum test
co-occurrence motif enrichment ...




### Additional functions

#### RBP functions

The selected RBPs can be further filtered by their annotated RBP functions prior to search.
Currently the following RBP functions appear in the database:

```
Function ID                     Function description
EJC                           Exon Junction Complex
MRP                             microRNA processing
MTR                    mitochondrial RNA regulation
PSG                        P-body / stress granules
RBT                    Ribosome & basic translation
RE                                       RNA export
RL                                 RNA localization
RM                                 RNA modification
RRP                                 rRNA processing
RSD                           RNA stability & decay
SNT                     snoRNA / snRNA / telomerase
SR                              Splicing regulation
SS                                      Spliceosome
TEP                               3' end processing
TR                           Translation regulation
TRR                                 tRNA regulation
VRR                            Viral RNA regulation
```

For example, to select all RBPs with annotated 3' end processing and RNA stability & decay functions, 
specify in any search mode `--rbps ALL --functions TEP RSD`.
To see which RBPs are assigned to which functions, run `rbpbench info`. 



#### GO term analysis

GO term (enrichment) analysis (GOA) can also be performed with RBPBench (see `--goa` and related options). 
What target genes are selected for enrichment in GO terms depends on the set mode and command line options.
E.g. in `rbpbench search`, by default GOA is performed on all genes covered by input regions. 
This can be changed via `--goa-cooc-mode`, allowing to select only genes with motif hits from any specified 
RPBs or from all specified RBPs to test different hypotheses.
This differs from `rbpbench batch`, where GOA is performed on genes with input regions from all input datasets
(further restricted via `--goa-only-cooc`, focussing only on genes with motif containing input regions from all datasets).
A third way to run GOA is offered by `rbpbench searchlongrna`, which searches transcripts (e.g. all or selected 
sets of mRNA transcripts) for motifs, and then uses the transcripts (i.e., its underlying genes) with motif hits 
(either by any RBP or all RBPs `--goa-only-cooc`) as target genes for GOA.
We can thus test diverse hypotheses by taking advantage of RBPBench's various program modes.
If you just want to perform GOA on some specific gene list, use `rbpbench goa`.
Some useful options are: `--goa-max-child` (e.g. `--goa-max-child 200` to filter out very broad GO terms),
`--goa-min-depth` (e.g. `--goa-min-depth 5` minimum depth for significant GO terms), 
or `--goa-pval` (set GOA p-value threshold).


### Helper scripts

Various helper scripts are included as well on the command line. 

