# RBPBench

[![GitHub](https://img.shields.io/github/tag/michauhl/RBPBench.svg)](https://github.com/michauhl/RBPBench)
[![Bioconda](https://anaconda.org/bioconda/rbpbench/badges/version.svg)](https://anaconda.org/bioconda/rbpbench)

RBPBench is multi-function tool to evaluate CLIP-seq and other genomic region data 
using a comprehensive collection of known RNA-binding protein (RBP) binding motifs. RBPBench can be used for a variety of
purposes, from RBP motif search (database or user-supplied RBPs) in genomic regions, over motif 
co-occurrence analysis, to benchmarking CLIP-seq peak caller methods as well as comparisons across cell types and CLIP-seq protocols.

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
conda 23.7.2
```


### Conda package installation

RBPBench is available as a Conda package, which makes installation a breeze. 
We simply create a Conda environment and install RBPBench inside the environment:

```
conda create -n rbpbench -c conda-forge -c bioconda rbpbench
conda activate rbpbench
```

If `conda install` is taking too long, one way to speed up the process is to replace the solver of Conda with the Mamba solver.

```
conda install -n base conda-libmamba-solver
conda config --set solver libmamba
```

RBPBench should now be available inside the environment:

```
rbpbench -h
```

### Manual installation

Manual installation of RBPBench is only slightly more work. First we create the Conda environment with all necessary dependencies:

```
conda create -n rbpbench -c conda-forge -c bioconda logomaker markdown meme scipy plotly textdistance venn matplotlib-venn infernal bedtools upsetplot scikit-learn
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

RBPBench is also available as a webserver on Galaxy (more information [here](https://backofenlab.github.io/RBPBench/)).


## Example runs

More information will be added soon.


## Documentation

This documentation provides further details on RBPBench (version 0.91).

### Program modes

RBPBench is a multi-function tool featuring several modes of operation. 
To get an overview of the currently available modes:

```
$ rbpbench -h
usage: rbpbench [-h] [-v]
                {search,batch,searchseq,searchregex,searchlong,searchrna,searchlongrna,optex,info,dist,compare}
                ...

Evaluate CLIP-seq and other genomic region data using a comprehensive
collection of known RBP binding motifs (RNA sequence + structure). RBPBench
can be used for a variety of purposes, from RBP motif search in genomic
regions, over motif co-occurrence analysis, to benchmarking CLIP-seq peak
callers, as well as comparisons across cell types and CLIP-seq protocols.

positional arguments:
  {search,batch,searchseq,searchregex,searchlong,searchrna,searchlongrna,optex,info,dist,compare}
                        Program modes
    search              Search motifs in genomic sites
    batch               Search motifs on > 1 dataset
    searchseq           Search motifs in sequences
    searchregex         Search regex in genomic sites or sequences
    searchlong          Search motifs in long genomic regions
    searchrna           Search motifs in spliced transcript sites
    searchlongrna       Search motifs in spliced full transcripts
    optex               Investigate optimal extension
    info                Print out RBP IDs in database
    dist                Plot nt distribution at genomic positions
    compare             Compare different search results

options:
  -h, --help            show this help message and exit
  -v, --version         show program's version number and exit

```

More information to be added soon.

