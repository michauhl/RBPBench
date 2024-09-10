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
conda create -n rbpbench -c conda-forge -c bioconda logomaker markdown meme scipy plotly textdistance venn matplotlib-venn infernal bedtools upsetplot scikit-learn goatools
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

More information to be added soon.

