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

To run the examples, we change to the cloned repository subfolder `RBPBench\test`
and download the genome FASTA file ([details](#genomic-sequences)) and a 
fitting GTF file from Ensembl ([details](#genomic-annotations)).

### Motif search in single input dataset

For motif search in a single set of genomic input regions (typically CLIP-seq), we can use ``



rbpbench search --in eclip_clipper_idr/SLBP_K562_IDR_peaks.bed --genome hg38.fa --gtf Homo_sapiens.GRCh38.112.gtf.gz --out test_search_slbp_out  --rbps ALL --ext 10 --functions TEP --regex AATAAA --bed-sc-thr 3 --cooc-pval-thr 0.001 --min-motif-dist 5 --add-all-reg-bar
rbpbench search --in eclip_clipper_idr/SLBP_K562_IDR_peaks.bed --genome hg38.fa --gtf Homo_sapiens.GRCh38.112.gtf.gz --out test_search_slbp_out  --rbps ALL --ext 10 --functions TEP --regex AATAAA --bed-sc-thr 3 --cooc-pval-thr 0.05 --min-motif-dist 5 --add-all-reg-bar

What table to see which genes have hnrnpk motif hits ?

Run parameter settings:
test_search_slbp_out/settings.rbpbench_search.out
RBP co-occurrence stats .tsv:
test_search_slbp_out/rbp_cooc_stats.tsv
Filtered input regions .bed:
test_search_slbp_out/in_sites.filtered.bed
Filtered input regions .fa:
test_search_slbp_out/in_sites.filtered.fa
Motif hits .bed:
test_search_slbp_out/motif_hits.rbpbench_search.bed
RBP region occupancies .tsv:
test_search_slbp_out/rbp_region_occupancies.tsv
RBP hit stats .tsv:
test_search_slbp_out/rbp_hit_stats.tsv
Motif hit stats .tsv:
test_search_slbp_out/motif_hit_stats.tsv
Region annotations .tsv:
test_search_slbp_out/region_annotations.tsv
Search report .html:
test_search_slbp_out/report.rbpbench_search.html


region_annotations.tsv

region_id	gene_id	gene_name	transcript_id	region_annotation	transcript_biotype
chr9:127813745-127813823(+)	ENSG00000136877	FPGS	ENST00000373247	3'UTR	protein_coding
chr6:26021705-26021766(+)	ENSG00000278637	H4C1	ENST00000617569	CDS	protein_coding
Add:
RBP_motifs




Then go for enmo/nemo




### Compare search results




```
unzip batch_compare_test.zip
```

```
$ cat batch_compare_test.batch_in.txt
PUM1	clipper_rep1	k562_eclip	batch_compare_test/PUM1.k562_eclip.clipper_rep1.bed
PUM1	clipper_rep2	k562_eclip	batch_compare_test/PUM1.k562_eclip.clipper_rep2.bed
PUM1	clipper_idr	k562_eclip	batch_compare_test/PUM1.k562_eclip.clipper_idr.bed
PUM1	dewseq_w100_s5	k562_eclip	batch_compare_test/PUM1.k562_eclip.dewseq_w100_s5.bed
RBFOX2	clipper_idr	hepg2_eclip	batch_compare_test/RBFOX2.hepg2_eclip.clipper_idr.bed
RBFOX2	clipper_idr	k562_eclip	batch_compare_test/RBFOX2.k562_eclip.clipper_idr.bed


rbpbench batch --bed batch_compare_test.batch_in.txt --genome hg38.fa --out batch_test_out
rbpbench compare --in batch_test_out --out compare_test_out
```




#### Plot nucleotide distribution at genomic positions

We can use `rbpbench dist` to plot the nucleotide distribution at genomic positions.
This can be used e.g. to check for potential nucleotide biases at CLIP-seq crosslink positions.


For illustration, we extract genomic stop codon positions from our ENSEMBL GTF file using 
the [helper script](#helper-scripts) `gtf_extract_tr_feat_bed.py`, 
and run `rbpbench dist` on the obtained stop codons BED file:

```
gtf_extract_tr_feat_bed.py --feat stop_codon --gtf Homo_sapiens.GRCh38.112.gtf.gz --out stop_codons.Homo_sapiens.GRCh38.112.bed --uniq-reg
rbpbench dist --in stop_codons.Homo_sapiens.GRCh38.112.bed --genome hg38.fa --out test_dist_out --ext 5
```

By default the upstream start position of the genomic input regions is used as 
position zero in the plot (change with `--cp-mode` option). The amount of up- and 
downstream context added is defined by `--ext`. 
The generated plot (`test_dist_out/nt_dist_zero_pos.png`) looks the following:

<img src="docs/nt_dist_zero_pos.stop_codons.png" alt="Nucleotide distribution at stop codons"
	title="Nucleotide distribution at stop codons" width="500" />

**Figure:** Nucleotide distribution (position probability matrix) at genomic stop codon positions (human transcript annotations, ENSEMBL GRCh38 release 112).

We can clearly identify the known stop codon triplet sequences (in DNA: TAA, TAG, TGA), starting 
at position 0.


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

`rbpbench search` can be used to search for RBP binding motifs in genomic regions, selecting 
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

If you have one or more sequence motifs of interest (MEME format or regex) and want to know if there are similar motifs 
in the database, you can run `rbpbench tomtom` (incorporating MEME's TOMTOM). 
This mode also informs us whether the reported motif hits are enriched in certain RBP functions (e.g., 
for regex input `--in TTTTTT`, the top 3 enriched functions are: splicing regulation, 
translation regulation, and RNA stability & decay). 

##### More modes

```rbpbench streme``` allows to discover new motifs using STREME from MEME suite. 
```rbpbench goa``` enables us to run GO term enrichment analysis (GOA) on a set of genes, e.g. obtained 
from the output tables of other modes (like genes covered by input regions from ```rbpbench search```
or ```rbpbench batch```). Note that some of the modes also can run GOA, optimized for the specific mode
(see [GO term analysis](#go-term-analysis)).
Furthermore, we can use `rbpbench dist` to plot the nucleotide distribution at genomic positions, 
and ```rbpbench info``` informs about database RBP motifs and annotated RBP functions.

##### Connecting modes

All motif search modes additionally output motif hit regions in BED format. These files for each hit also contain the matched sequence, 
plus the genomic region annotation (if a GTF file was supplied). This way, the motif hit regions themselves 
(or a subset of interest after filtering, e.g. certain matched sequences or only motifs residing in 3'UTRs)
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
which enable enrichment and co-occurrence statistics on single motif level ([details](#co-occurrence-statistics)). 
Alternatively, single motifs for search can be selected via `--motifs`, e.g. `--motifs CSTF1_1 DDX3X_1`.
To list and visualize all selected motifs, use `--plot-motifs` (in some modes automatically output).
RBP motifs can also be filtered by their annotated RBP functions. To only inlcude RBPs with e.g. annotated 3' end processing
function, set `--functions TEP` (for available functions and annotations, run `rbpbench info`).
Moreover, if you have a motif of interest (or a regex), and want to know if there are similar motifs in the database, you can 
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
directly on the transcript regions, which enables us to examine motif binding directly on the 
spliced transcripts.

#### Genomic sequences

Genomic sequences have to be provided in FASTA format. For example, the human 
genome sequence (`hg38.fa`) can be obtained from [here](https://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips).
To download in the command line and unpack:

```
wget https://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/hg38.fa.gz
gunzip hg38.fa.gz
```

#### Genomic annotations

Many of the statistics require a genomic annotations file in GTF format in order to be generated. 
RBPBench was tested mainly using GTF files downloaded from [Ensembl](http://www.ensembl.org/info/data/ftp/index.html), 
but you can also provide a GTF file from  e.g. [GENCODE](https://www.gencodegenes.org/human/). 
For the examples we downloaded the following GTF file:

```
get https://ftp.ensembl.org/pub/release-112/gtf/homo_sapiens/Homo_sapiens.GRCh38.112.gtf.gz
```

#### User-provided motif search

Both sequence (MEME XML format) and structure (covariance model .cm) motifs can be supplied by the user 
on top of the database RBPs (using `-rbps USER` option together with `-user-meme-xml` or `--user-cm`).
This way, the motifs can be used together with the database motifs for search and comparative statistics.
For example, to supply a structure motif (SLBP) via `--user-cm` (the example motif files can be found in the RBPBench repository subfolder `RBPBench/test`):

```
rbpbench search --in SLBP_K562_IDR_peaks.bed --rbps USER --out SLBP_user_search_out --genome hg38.fa --user-cm path_to_test/SLBP_USER.cm  --user-rbp-id SLBP_USER
```

In the same way, we can supply sequence motif(s) (PUM1) via `--user-meme-xml`, and e.g. also combine it with any of the database RBPs (here PUM2 and RBFOX2):

```
rbpbench search --in PUM1_K562_IDR_peaks.bed --rbps USER PUM2 RBFOX2 --out PUM1_user_search_out --genome hg38.fa --user-meme-xml path_to_test/PUM1_USER.xml --user-rbp-id PUM1_USER
```

#### Custom motif database

Apart from the built-in motif database and the option of user-supplied RBP motifs, 
it is also possible to define a custom motif database, which can then be applied just like the built-in motif database. 
The following command line parameters deal with defining a custom motif database:

```
  --custom-db str       Provide custom motif database folder. Alternatively, provide single files via --custom-db-meme-xml, --custom-
                        db-cm, --custom-db-info
  --custom-db-id str    Set ID/name for provided custom motif database via --custom-db (default: "custom")
  --custom-db-meme-xml str
                        Provide custom motif database MEME/DREME XML file containing sequence motifs
  --custom-db-cm str    Provide custom motif database covariance model (.cm) file containing covariance model(s)
  --custom-db-info str  Provide custom motif database info table file containing RBP ID -> motif ID -> motif type assignments
```

The database can be input as a folder (`--custom-db db_folder_path`), which needs to contain an `info.txt` file, 
as well as a sequence motifs file (`seq_motifs.meme`, MEME motif format), and/or
a structure motifs file (`str_motifs.cm`, covariance model format). 
The `info.txt` is a table file containing the RBP ID to motif ID mappings. 
Here is an example of a valid `info.txt` (minimum 3 columns required: RBP_motif_ID, RBP_name, Motif_type) file content:

```
$ cat db_folder_path/info.txt
RBP_motif_ID	RBP_name	Motif_type
A1CF_1	A1CF	meme_xml	human
A1CF_2	A1CF	meme_xml	human
ACIN1_1	ACIN1	meme_xml	human
ACIN1_2	ACIN1	meme_xml	human
ACO1_1	ACO1	meme_xml	human
RF00032	SLBP	cm	human
```

`Motif_type` defines whether a motif is a sequence motif (expected to be found in `seq_motifs.meme`), 
or a structure motif (expected to be found in `str_motifs.cm`).
An ID / name for the custom database can be defined as well via `--custom-db-id`.
Alternatively, the files can be input separately via `--custom-db-info`, 
`--custom-db-meme-xml`, and `--custom-db-cm`.

If you have some short sequences or regular expressions (regexes), you can also quickly create a MEME format motif 
database out of these. Only thing needed is a table file containing regexes and associated RBP/motif IDs:

```
$ cat custom_motifs.tsv
rbp_id	motif_id	regex
RBPX	RBPX_1	AAA[CG]A
RBPX	RBPX_2	AAA[CGT]C
RBPY	RBPY_1	ACGT[ACGT]
```

To create the database folder, we use one of the helper scripts:

```
create_custom_meme_motif_db.py --in custom_motifs.tsv --out custom_db_out
```

Afterwards the folder can be used as described above via `--custom-db` option, and you can specify any RBPs or motifs you want to search for from the database:

```
rbpbench search --in genomic_regions.bed --genome /path/to/hg38.fa --gtf /path/to/Homo_sapiens.GRCh38.112.gtf.gz --rbps ALL --custom-db custom_db_out --out custom_db_search_out --ext 10 --plot-motifs --fimo-pval 0.005
```

Here we used `--plot-motifs` to visualize the motifs and `--rbps ALL`, meaning all motifs from the custom database are included in search. To search only for specific motifs from the specified database, simply add `--motifs RBPX_1 RBPX_2`. Note that RBPBench's default FIMO threshold to filter out motif hits is optimized for 6-nt or longer sequence motifs. If your motifs are shorter, consider setting a more relaxed threshold (e.g., like in the example `--fimo-pval 0.005`). However, keep in mind that a relaxed setting also results in predicting more non-perfect hits, especially for longer motifs.
You can easily check by inputting sequences with motifs you expect to be predicted into `rbpbench searchseq` with `--custom-db` option.


#### FIMO nucleotide frequencies

What p-value a motif hit reported by FIMO gets depends on the set background nucleotide frequencies (FIMO option: `--bfile`). 
By default, frequencies obtained from human ENSEMBL transcripts (excluding introns, A most prominent nucleotide)
are used. You can change the preset via `--fimo-ntf-mode` (other options: transcripts including introns, or uniform distribution). 
Alternatively, you can provide your own background frequencies file via `--fimo-ntf-file`.


### Outputs

RBPBench depending on the mode outputs various useful output files (e.g. motif hits BED, 
region annotations, hit statistics files, HTML reports with statistics and visualizations).

#### Hit statistics table files

RBPBench's motif search modes (e.g. `rbpbench search`, `rbpbench batch`) 
output hit statistics table files which can later be used to 
do comparisons between peak callers or other specified conditions (e.g., cell types, CLIP-seq protocols).
RBPBench  outputs the RBP binding motif hit statistics into two table files 
stored in the results output folder (`--out`):

```
RBP hit stats .tsv:
results_output_folder/rbp_hit_stats.tsv
Motif hit stats .tsv:
results_output_folder/motif_hit_stats.tsv
```

The first file (RBP hit stats) contains comprehensive hit statistics for each RBP 
(one row per RBP), 
while the motif hits stats file contains hit statistics for each single motif hit 
(one row per motif hit).

The RBP hit statistics file `rbp_hit_stats.tsv` contains the following columns:


| Column name | Description |
|:--------------:|:--------------:|
| data_id | Set Data ID (`--data-id`). More [here](#adding-more-information-for-comparisons) |
| method_id | Set method ID (`--method-id`). More [here](#adding-more-information-for-comparisons) |
| run_id | Set run ID (`--run-id`) |
| motif_db | Selected motif database for search (`--motif-db`) |
| rbp_id | RBP ID (i.e., RBP name), e.g. `PUM1` |
| c_regions | Number of input genomic regions used for search (after filtering and extension operations) |
| mean_reg_len | Mean length of genomic regions |
| median_reg_len | Median length of genomic regions |
| median_reg_len | Median length of genomic regions |
| min_reg_len | Minimum length of genomic regions |
| called_reg_size | Called size of genomic regions (including overlaps) |
| effective_reg_size | Effective size of genomic regions (removed overlaps) |
| c_reg_with_hits | Number of regions with motif hits from rbp_id |
| perc_reg_with_hits | Percentage of regions with motif hits from rbp_id |
| c_motif_hits | Total number of motif hits from rbp_id |
| c_uniq_motif_hits | Number of unique motif hits from rbp_id (removed double counts) |
| c_uniq_motif_nts | Number of unique motif nucleotides from rbp_id (removed overlaps) |
| perc_uniq_motif_nts_cal_reg | Percentage of unique motif nucleotides over called region size |
| perc_uniq_motif_nts_eff_reg | Percentage of unique motif nucleotides over effective region size |
| uniq_motif_hits_cal_1000nt | Number of motif hits over 1000 nt of called region size |
| uniq_motif_hits_eff_1000nt | Number of motif hits over 1000 nt of effective region size |
| wc_pval | [Wilcoxon rank-sum test p-value](#informative-statistics) to test whether motif hit regions tend to feature higher scores | 
| seq_motif_ids | Sequence motif IDs. Empty (`-`) if rbp_id has not sequence motifs  | 
| seq_motif_hits | Sequence motif hit counts (count for each motif ID) | 
| str_motif_ids | Structure motif IDs. Empty (`-`) if rbp_id has not structure motifs  | 
| str_motif_hits | Structure motif hit counts (count for each motif ID) | 
| internal_id | Internal ID (unique for each rbp_id run), used for connecting table results  | 


The motif hit statistics file `motif_hit_stats.tsv` contains the following columns:

| Column name | Description |
|:--------------:|:--------------:|
| data_id | Set Data ID (`--data-id`). More [here](#adding-more-information-for-comparisons) |
| method_id | Set method ID (`--method-id`). More [here](#adding-more-information-for-comparisons) |
| run_id | Set run ID (`--run-id`) |
| motif_db | Selected motif database for search (`--motif-db`) |
| region_id | Genomic region ID containing the hit, e.g. `chr1:228458485-228458560(+)` |
| rbp_id | RBP ID (i.e., RBP name), e.g. `SLBP` |
| motif_id | Motif ID |
| chr_id | chromosome ID |
| gen_s | genomic motif hit start (1-based) |
| gen_e | genomic motif hit end (1-based) |
| strand | Chromosome strand (+ or - strand) |
| region_s | region motif hit start (1-based) |
| region_e | region motif hit end (1-based) |
| region_len | Region length |
| uniq_count | Unique motif hit count |
| fimo_score | FIMO score (for sequence motif hits) |
| fimo_pval | FIMO p-value (for sequence motif hits) |
| cms_score | cmsearch score (for structure motif hits) |
| cms_eval | cmsearch e-value (for structure motif hits) |
| matched_seq | matched motif hit sequence |
| internal_id | Internal ID (unique for each rbp_id run), used for connecting table results  | 


#### Motif hits BED files

Motif hits are also output in BED format in each search mode. 
This together with the input regions BED allows for quick studying of motif occurrences 
in a genome viewer (e.g., [IGV](https://software.broadinstitute.org/software/igv/)). 
Additionally, the motif hits BED file can be further filtered by matched sequence (BED column 12) 
or genomic region annotation (BED column 11 if `--gtf` was provided, e.g. filter by `3'UTR` to keep only
3'UTR hits). These motif hit regions can then again be used as search mode input regions, allowing 
for refined hypothesis testing.

There are various settings for controlling motif hit search. E.g., a FIMO p-value threshold for reporting
sequence motif hits (`--fimo-pval`), or a CMSEARCH bit score threshold for reporting structure motif 
hits (`--cmsearch-bs`, also see `--cmsearch-mode`). Furthermore, `--greatest-hits` enabled 
results in only the best hit (lowest p-value, highest bit score) being reported for each input region.



### Statistics

#### Input region score motif enrichment statistics

Given a set of input regions with associated scores (by default BED column 5 is used, change via `--bed-score-col`),
RBPBench for each selected RBP checks whether motif-containing input regions 
have significantly higher scores (using Wilcoxon rank sum test).
This means that a low test p-value for a given RBP indicates that higher-scoring regions are more likely to contain motif hits of the respective RBP. If we assume that the score (e.g. log2 fold change) is somehow correlated with RBP binding affinity,
the test thus can give clues on which RBPs preferentially bind to the provided regions.
Note that the test is only informative if the scores are themselves informative w.r.t. RBP binding (e.g., not all the same), 
or e.g. if not all or too many of the input regions contain motif hits. 

The test results are output in the [output tables](#hit-statistics-table-files), as well as in the HTML reports.
The test can also check for significantly lower scores (change via `--wrs-mode`). Moreover, in `rbpbench batch`, 
the test can be applied for a user-defined regex (via `--regex`). This conveniently allows us to check, 
e.g. over all eCLIP datasets, whether a given regex is significantly enriched in any dataset (results again output to HTML report).


#### Co-occurrence statistics

To test for the co-occurrence of RBP motifs in a set of input regions, 
Fisher's exact test is applied, effectively counting the number of occupied and non-occupied regions 
for each pair of RBPs. This is done on the RBP level (`rbpbench search`, `rbpbench searchrna`),
meaning that all motifs assigned to the same RBP are aggregated in the statistic, 
but can also be done on the single motif level for a more fine-grained analysis (`rbpbench enmo`, `rbpbench nemo`). 
Co-occurrence results are output as interactive heat maps into the respective HTML reports.
Co-occurrence p-value threshold and multiple testing correction method can be adapted via 
`--cooc-pval-thr` and `--cooc-pval-mode` (default: Benjamini Hochberg). 
Moreover, Fisher exact test alternative hypothesis can be changed (`--fisher-mode`) from greater to 
two-sided or less, so instead of reporting
significantly higher co-occurrences, we can e.g. instead report significantly lower co-occurrences 
(i.e., the RBPs or motifs that tend NOT to co-occur together).
Since reported motifs between RBPs can often be very similar, we can further set a minimum mean distance 
between the motif hits required for them to be reported as significant (`--min-motif-dist`). 
This allows us to ignore significant co-occurrences between very similar motifs 
(i.e., they are reported as not significant in the heat map). Another way of influencing co-occurrences
is to filter the input regions by their region score (if the score is meaningful), which can 
be done via `--bed-sc-thr`.

Considering the single motif level co-occurrences (`rbpbench enmo`, `rbpbench nemo`), 
co-occurrence attribution is more strict, since it is only checked for motifs which are significantly enriched 
in the input regions (details [here](#input-region-motif_enrichment-statistics)). 
Furthermore, in addition to `--min-motif-dist`, co-occurrences can be filtered 
by motif pair similarity (only for MEME motif formatted sequence motifs, `--motif-sim-thr`), which 
allows us to focus on co-occurrences with more dissimilar motifs. 


#### Input region motif enrichment statistics

Not to be confused with the region score statistic [above](#input-region-score-motif_enrichment-statistics), 
input region motif enrichment statistics are calculated in the modes `rbpbench enmo`, `rbpbench nemo`, 
using the motif occurrences in the input region dataset and comparing it to their occurrences 
in a background dataset (again via Fisher's exact test). The background set can be generated in 
various ways, either by shuffling the input sequences (`--bg-mode 2`), or by sampling background 
regions from the genome or transcriptome (default: `--bg-mode 1`if input regions are on transcripts).
Various options to modify background set generation exist: 
`--bg-min-size` to sample additional background regions; 
`--bg-mask-bed` to define regions from which no background regions should be sampled;
`--bg-incl-bed` to define regions from which to sample background regions;
`--bg-shuff-factor`, `--bg-shuff-k` to customize shuffling;
and some more options (`--bg-mask-blacklist`, `--bg-ada-sampling`, `--random-seed`, 
check `rbpbench enmo -h` for details).
Similar to the co-occurrence statistics, p-value threshold and multiple testing correction method 
can be adapted (`--enmo-pval-thr`, `--enmo-pval-mode`, `--nemo-pval-thr`, `--nemo-pval-mode`), 
and `--fisher-mode` is also available.


### Additional information

#### Handling genomic overlaps

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


#### Most-prominent transcript selection

In order to obtain genomic region annotations on the transcript level, 
by default one representative transript (i.e., the most prominent transcript (MPT))
is chosen for each gene in the GTF file. MPT selection aims to select the transcript 
with the strongest experimental support. Currently, for a pair of transcripts from 
the same gene, to decide which one is "more prominent":

1. Chose transcript with "basic" tag
2. If both have "basic" tag, chose transcript with higher transcript support level (TSL)
3. If same TSL, chose transcript with "Ensembl_canonical" tag
4. If both have "Ensembl_canonical" tag, chose the longer transcript as current MPT

This is done for all possible pairs, which leaves one transcript as MPT in the end for each gene.
Regions that do not overlap with the selected transcript regions are assigned to "intergenic".
Alternatively, a list of transcript IDs can be supplied `--tr-list`, bypassing the MPT selection.
Which region annotations are to be considered can be further defined via `--tr-types`, and the 
minimum region annotation overlap can be set via `--gtf-feat-min-overlap`.


#### Helper scripts

Various helper scripts are included as well on the command line:

```
bed_print_first_n_pos.py
bed_print_last_n_pos.py
bed_shift_regions.py
create_custom_meme_motif_db.py
gtf_extract_exon_intron_border_bed.py
gtf_extract_exon_intron_region_bed.py
gtf_extract_gene_region_bed.py
gtf_extract_mpt_region_bed.py
gtf_extract_tr_feat_bed.py
gtf_get_gene_region_nt_freqs.py
gtf_get_mpt_nt_freqs.py
gtf_get_mpt_with_introns_nt_freqs.py
```
To can call their help pages to get more infos on what they do and how to use them (e.g., `bed_print_first_n_pos.py -h`).
To get a quick overview: 
`bed_print_first_n_pos.py` prints the first n positions of each region from the provided BED file.
`bed_print_last_n_pos.py` prints the last n positions of each region from the provided BED file.
`bed_shift_regions.py` shifts all regions from the provided BED file a given number of nucleotides up- or downstream.
`create_custom_meme_motif_db.py` as described above generates a custom sequence motifs database which can be used in all search modes via `--custom-db`.
`gtf_extract_exon_intron_border_bed.py` extracts exon-intron border positions from a GTF file and stores them in a BED file, 
which can be used as input e.g. in `rbpbench nemo`.
`gtf_extract_exon_intron_region_bed.py` extracts exon or intron regions from a GTF file and stores them in a BED file.
`gtf_extract_gene_region_bed.py` extracts gene regions from a GTF file and stores them in a BED file.
`gtf_extract_mpt_region_bed.py` extracts most prominent transcript regions from a GTF file and stores them in a BED file. 
Additionally, mRNA regions (5'UTR, CDS, 3'UTR) can be output to a separate BED file.
`gtf_extract_tr_feat_bed.py` extracts transcript feature regions from a GTF file and stores them in a BED file (e.g. stop_codon).
`gtf_get_gene_region_nt_freqs.py` calculates nucleotide frequencies from all gene regions extracted from a GTF and the corresponding genome FASTA.
FIMO can be given this information as a nucleotide frequencies file (see options `--fimo-ntf-file`, `--fimo-ntf-mode`).
`gtf_get_mpt_nt_freqs.py` calculates nucleotide frequencies of from all most prominent transcript (MPT) sequences (introns excluded) 
extracted from a GTF and the corresponding genome FASTA.
`gtf_get_mpt_with_introns_nt_freqs.py` calculates nucleotide frequencies of from all most prominent transcript (MPT) sequences (introns included) 
extracted from a GTF and the corresponding genome FASTA.
