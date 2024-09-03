#!/usr/bin/env python3

import argparse
import os
import re
import subprocess


###############################################################################

def setup_argument_parser():
    """Setup argparse parser."""
    help_description = """
    Python wrapper for RBPBench Galaxy wrapper to work with collections of
    input BED files (i.e. to process them with rbpbench batch).
    """
    # Define argument parser.
    p = argparse.ArgumentParser(add_help=False,
                                prog="batch_table_wrapper_rbpbench.py",
                                description=help_description,
                                formatter_class=argparse.MetavarTypeHelpFormatter)

    # Required arguments.
    p.add_argument("-h", "--help",
                   action="help",
                   help="Print help message")
    p.add_argument("--table",
                   dest="in_table",
                   type=str,
                   metavar='str',
                   required=True,
                   help="Input table file with data ID, method ID, RBP ID and file name (Galaxy element identifier in dataset collection) for each to be processed dataset by rbpbench batch")
    p.add_argument("--paths",
                   dest="in_paths",
                   type=str,
                   metavar='str',
                   nargs='+',
                   required=True,
                   help="List of Galaxy BED file paths (--files path1 path2 .. )")
    p.add_argument("--ids",
                   dest="in_ids",
                   type=str,
                   metavar='str',
                   nargs='+',
                   required=True,
                   help="List of Galaxy element identifiers, equal to the BED dataset names in the dataset collection (--ids id1 id2 .. )")
    p.add_argument("--genome",
                   dest="in_genome",
                   type=str,
                   metavar='str',
                   required=True,
                   help="Genomic sequences file (currently supported formats: FASTA)")
    p.add_argument("--out",
                   dest="out_folder",
                   type=str,
                   metavar='str',
                   required=True,
                   help="Batch results output folder")
    # Optional batch arguments.
    p.add_argument("--ext",
                   dest="ext_up_down",
                   type=str,
                   metavar='str',
                   default="0",
                   help="Up- and downstream extension of --in sites in nucleotides (nt). Set e.g. --ext 30 for 30 nt on both sides, or --ext 20,10 for different up- and downstream extension (default: 0)")
    p.add_argument("--motif-db",
                   dest="motif_db",
                   type=int,
                   default=1,
                   choices=[1, 2, 3],
                   help="Built-in motif database to use (currently there are 3 human databases: 1,2,3. See details in rbpbench (default: 1)")
    p.add_argument("--fimo-ntf-file",
                   dest="fimo_user_ntf_file",
                   type=str,
                   metavar='str',
                   default = False,
                   help = "Provide FIMO nucleotide frequencies (FIMO option: --bfile) file (default: use internal frequencies file, define which with --fimo-ntf-mode)")
    p.add_argument("--fimo-ntf-mode",
                   dest="fimo_ntf_mode",
                   type=int,
                   default=1,
                   choices=[1, 2, 3],
                   help="Set which internal nucleotide frequencies to use for FIMO search. 1: use frequencies from human ENSEMBL transcripts (excluding introns, A most prominent) 2: use frequencies from human ENSEMBL transcripts (including introns, resulting in lower G+C and T most prominent) 3: use uniform frequencies (same for every nucleotide) (default: 1)")
    p.add_argument("--fimo-pval",
                   dest="fimo_pval",
                   type=float,
                   metavar='float',
                   default=0.001,
                   help="FIMO p-value threshold (FIMO option: --thresh) (default: 0.001)")
    p.add_argument("--cmsearch-bs",
                   dest="cmsearch_bs",
                   type=float,
                   metavar='float',
                   default=1.0,
                   help="CMSEARCH bit score threshold (CMSEARCH options: -T --incT). The higher the more strict (default: 1.0)")
    p.add_argument("--cmsearch-mode",
                   dest="cmsearch_mode",
                   type=int,
                   default=1,
                   choices=[1, 2],
                   help="Set CMSEARCH mode to control strictness of filtering. 1: default setting (CMSEARCH option: --default). 2: max setting (CMSEARCH option: --max), i.e., turn all heuristic filters off, slower and more sensitive / more hits) (default: 1)")
    p.add_argument("--greatest-hits",
                   dest="greatest_hits",
                   default = False,
                   action = "store_true",
                   help = "Keep only best FIMO/CMSEARCH motif hits (i.e., hit with lowest p-value / highest bit score for each motif sequence/site combination). By default, report all hits (default: False)")
    p.add_argument("--bed-score-col",
                   dest="bed_score_col",
                   type=int,
                   metavar='int',
                   default=5,
                   help="--in BED score column used for p-value calculations. BED score can be e.g. log2 fold change or -log10 p-value of the region (default: 5)")
    p.add_argument("--bed-sc-thr",
                   dest="bed_sc_thr",
                   type = float,
                   metavar='float',
                   default = None,
                   help = "Minimum site score (by default: --in BED column 5, or set via --bed-score-col) for filtering (assuming higher score == better site) (default: None)")
    p.add_argument("--bed-sc-thr-rev",
                   dest="bed_sc_thr_rev_filter",
                   default = False,
                   action = "store_true",
                   help = "Reverse --bed-sc-thr filtering (i.e. the lower the better, e.g. for p-values) (default: False)")
    p.add_argument("--unstranded",
                   dest="unstranded",
                   default=False,
                   action="store_true",
                   help="Set if --in BED regions are NOT strand-specific, i.e., to look for motifs on both strands of the provided regions. Note that the two strands of a region will still be counted as one region (change with --unstranded-ct) (default: False)")
    p.add_argument("--unstranded-ct",
                   dest="unstranded_ct",
                   default=False,
                   action="store_true",
                   help="Count each --in region twice for RBP hit statistics when --unstranded is enabled. By default, two strands of one region are counted as one region for RBP hit statistics")
    p.add_argument("--meme-no-check",
                   dest="meme_disable_check",
                   default = False,
                   action = "store_true",
                   help = "Disable MEME version check. Make sure --meme-no-pgc is set if MEME version >= 5.5.4 is installed! (default: False)")
    p.add_argument("--meme-no-pgc",
                   dest="meme_no_pgc",
                   default = False,
                   action = "store_true",
                   help = "Manually set MEME's FIMO --no-pgc option (required for MEME version >= 5.5.4). Make sure that MEME >= 5.5.4 is installed! (default: False)")
    # Test modes.
    p.add_argument("--wrs-mode",
                   dest="wrs_mode",
                   type=int,
                   default=1,
                   choices=[1, 2],
                   help="Defines Wilcoxon rank sum test alternative hypothesis for testing whether motif-containing regions have significantly different scores. 1: test for higher (greater) scores, 2: test for lower (less) scores (default: 1)")
    p.add_argument("--fisher-mode",
                   dest="fisher_mode",
                   type=int,
                   default=1,
                   choices=[1, 2, 3],
                   help="Defines Fisher exact test alternative hypothesis for testing co-occurrences of RBP motifs. 1: greater, 2: two-sided, 3: less (default: 1)")
    # Report.
    p.add_argument("--kmer-size",
                   dest="kmer_size",
                   type=int,
                   metavar='int',
                   default=5,
                   help="K-mer size for comparative plots (default: 5)")
    p.add_argument("--no-occ-heatmap",
                   dest="no_occ_heatmap",
                   default = False,
                   action = "store_true",
                   help = "Do not produce gene region occupancy heatmap plot in HTML report (default: False)")
    p.add_argument("--disable-heatmap-cluster-olo",
                   dest="disable_heatmap_cluster_olo",
                   default = False,
                   action = "store_true",
                   help="Disable optimal leave ordering (OLO) for clustering gene region occupancy heatmap. By default, OLO is enabled")
    p.add_argument("--report-header",
                   dest="report_header",
                   default = False,
                   action = "store_true",
                   help = "Add RBPBench to HTML report header. Useful for HTML reports inside Galaxy (default: False)")
    # Regex.
    p.add_argument("--regex",
                   dest="regex",
                   type=str,
                   metavar='str',
                   help="Define regular expression (regex) DNA motif to include in search, e.g. --regex AAACC, --regex 'C[ACGT]AC[AC]', ..")
    p.add_argument("--regex-search-mode",
                   dest="regex_search_mode",
                   type=int,
                   default=1,
                   choices=[1, 2],
                   help="Define regex search mode. 1: when motif hit encountered, continue +1 after motif hit end position, 2: when motif hit encountered, continue +1 of motif hit start position (default: 1)")
    p.add_argument("--max-motif-dist",
                   dest="max_motif_dist",
                   type=int,
                   metavar='int',
                   default=20,
                   help="Set maximum motif distance for regex-RBP co-occurrence statistic in HTML report (default: 20)")
    # GTF.
    p.add_argument("--gtf",
                   dest="in_gtf",
                   type=str,
                   metavar='str',
                   default = False,
                   help = "Input GTF file with genomic annotations to generate genomic regions annotation plots for each input BED file (output to HTML report). By default the most prominent transcripts will be extracted and used for functional annotation. Alternatively, provide a list of expressed transcripts via --tr-list (together with --gtf containing the transcripts). Note that only features on standard chromosomes (1,2,..,X,Y,MT) are currently used for annotation")
    p.add_argument("--tr-list",
                   dest="tr_list",
                   type=str,
                   metavar='str',
                   default = False,
                   help = "Supply file with transcript IDs (one ID per row) to define which transcripts to use from --gtf for genomic regions annotations plots")
    p.add_argument("--tr-types",
                   dest="tr_types_list",
                   type=str,
                   metavar='str',
                   nargs='+',
                   help="List of transcript biotypes to consider in genomic regions annotations plot. By default an internal selection of transcript biotypes is used (in addition to intron, CDS, UTR, intergenic). Note that provided biotype strings need to be in --gtf GTF file")
    p.add_argument("--gtf-feat-min-overlap",
                   dest="gtf_feat_min_overlap",
                   type=float,
                   metavar='float',
                   default=0.1,
                   help="Minimum amount of overlap required for a region to be assigned to a GTF feature (if less or no overlap, region will be assigned to \"intergenic\"). If there is overlap with several features, assign the one with highest overlap (default: 0.1)")
    p.add_argument("--gtf-eib-min-overlap",
                   dest="gtf_eib_min_overlap",
                   type=float,
                   metavar='float',
                   default=0.9,
                   help="Minimum amount input region has to overlap with exon (e), intron (i), i + ei borders to be counted as overlapping with these (note that the amount is reciprocal, i.e., one of the overlapping parts meeting the minimum amount is enough) (default: 0.9)")
    p.add_argument("--gtf-intron-border-len",
                   dest="gtf_intron_border_len",
                   type=int,
                   metavar='int',
                   default=250,
                   help="Set intron border region length (up- + downstream ends) for exon intron overlap statistics (default: 250)")
    # GO enrichment analysis for batch mode.
    p.add_argument("--goa",
                   dest="run_goa",
                   default = False,
                   action = "store_true",
                   help = "Run gene ontology (GO) enrichment analysis on genes occupied by sites in input datasets. Requires --gtf (default: False)")
    p.add_argument("--goa-obo-mode",
                   dest="goa_obo_mode",
                   type=int,
                   default=1,
                   choices=[1, 2, 3],
                   help = "Define how to obtain GO DAG (directed acyclic graph) obo file. 1: download most recent file from internet,  2: use local file,  3: provide file via --goa-obo-file (default: 1)")
    p.add_argument("--goa-obo-file",
                   dest="goa_obo_file",
                   type=str,
                   metavar='str',
                   default = False,
                   help = "Provide GO DAG obo file (default: False)")
    p.add_argument("--goa-gene2go-file",
                   dest="goa_gene2go_file",
                   type=str,
                   metavar='str',
                   default = False,
                   help = "Provide gene ID to GO IDs mapping table (row format: gene_id<tab>go_id1,go_id2). By default, a local file with ENSEMBL gene IDs is used. NOTE that gene IDs need to be compatible with --gtf (default: False)")
    p.add_argument("--goa-pval",
                   dest="goa_pval",
                   type=float,
                   metavar='float',
                   default=0.05,
                   help="GO enrichment analysis p-value threshold (applied on corrected p-value) (default: 0.05)")
    p.add_argument("--goa-only-cooc",
                   dest="goa_only_cooc",
                   default = False,
                   action = "store_true",
                   help = "Only look at genes in GO enrichment analysis which contain motif hits for all input datasets. By default, GO enrichment analysis is performed on the genes covered by sites from all input datasets (default: False)")
    p.add_argument("--goa-bg-gene-list",
                   dest="goa_bg_gene_list",
                   type=str,
                   metavar='str',
                   default = False,
                   help = "Supply file with gene IDs (one ID per row) to use as background gene list for GOA. NOTE that gene IDs need to be compatible with --gtf (default: False)")
    p.add_argument("--goa-max-child",
                   dest="goa_max_child",
                   type=int,
                   metavar='int',
                   default=None,
                   help="Specify maximum number of children for a significant GO term to be reported in HTML table, e.g. --goa-max-child 100. This allows filtering out very broad terms (default: None)")
    p.add_argument("--goa-min-depth",
                   dest="goa_min_depth",
                   type=int,
                   metavar='int',
                   default=None,
                   help="Specify minimum depth number for a significant GO term to be reported in HTML table, e.g. --goa-min-depth 5 (default: None)")
    p.add_argument("--goa-filter-purified",
                   dest="goa_filter_purified",
                   default = False,
                   action = "store_true",
                   help = "Filter out GOA results labeled as purified (i.e., GO terms with significantly lower concentration) in HTML table (default: False)")

    return p


################################################################################

def is_valid_regex(regex):
    """
    Check if regex string is valid regex.

    >>> is_valid_regex(".*")
    True
    >>> is_valid_regex(".*[")
    False
    >>> is_valid_regex("ACGT")
    True

    """

    try:
        re.compile(regex)
        return True
    except re.error:
        return False


################################################################################

def remove_special_chars_from_str(check_str,
                                  reg_ex='[^A-Za-z0-9_-]+'):
    """
    Remove special characters from string.

    reg_ex:
        Regular expression defining what to keep.

    >>> check_str = "{_}[-](_)\V/"
    >>> remove_special_chars_from_str(check_str)
    '_-_V'
    >>> check_str = ""
    >>> remove_special_chars_from_str(check_str)
    ''

    """
    clean_string = re.sub(reg_ex, '', check_str)
    return clean_string


###############################################################################

if __name__ == '__main__':

    parser = setup_argument_parser()
    args = parser.parse_args()

    assert os.path.exists(args.in_table), "--table file \"%s\" not found" % (args.in_file)
    assert os.path.exists(args.in_genome), "--genome file \"%s\" not found" % (args.in_genome)

    c_paths = len(args.in_paths)
    c_ids = len(args.in_ids)
    assert c_paths == c_ids, "given # paths (--paths) != # ids (--ids) (%i != %i). Please provide one ID for each path" % (c_paths, c_ids)

    """
    Check given paths and IDs.

    """

    # Paths.
    paths_dic = {}
    paths_list = []
    for path in args.in_paths:
        assert os.path.exists(path), "--paths %s file not found" % (path)
        if path not in paths_dic:
            paths_dic[path] = 1
        else:
            assert False, "--paths %s given > 1. Please provide unique paths" % (path)
        paths_list.append(path)

    # IDs
    ids_dic = {}
    ids_list = []
    for id in args.in_ids:
        if id not in ids_dic:
            ids_dic[id] = 1
        else:
            assert False, "--ids \"%s\" given > 1. Please provide unique element identifiers (dataset names) inside the dataset collection, in order to unambiguously assign element ID to file path" % (id)
        ids_list.append(id)

    id2path_dic = {}
    for idx, id in enumerate(ids_list):
        path = paths_list[idx]
        id2path_dic[id] = path

    """
    Read in table.

    Column format:
    rbp_id method_id data_id dataset_name

    """

    comb_ids_dic = {}
    id_collect_dic = {}
    id_collect_dic["rbp_id"] = []
    id_collect_dic["method_id"] = []
    id_collect_dic["data_id"] = []
    id_collect_dic["set_name"] = []
    id_collect_dic["path"] = []  # Galaxy file path.

    print("Read in --table ... ")

    with open(args.in_table) as f:
        for line in f:

            if re.search("^#", line):
                continue

            cols = line.strip().split("\t")

            assert len(cols) == 4, "line in --table with # cols != 4 (%i) encountered:%s" % (len(cols), line)

            rbp_id = cols[0]
            method_id = cols[1]
            data_id = cols[2]
            set_name = cols[3]

            if rbp_id == "rbp_id":
                continue

            comb_id = "%s,%s,%s,%s" % (rbp_id, method_id, data_id, set_name)

            if comb_id not in comb_ids_dic:
                comb_ids_dic[comb_id] = 1
            else:
                assert False, "data combination (\"%s\") appears > 1 in --table file. Please provide unique combinations for rbpbench batch calculation" % (comb_id)

            assert set_name in ids_dic, "given dataset name \"%s\" from --table not part of given --ids. Please provide dataset names present in dataset collection" % (set_name)

            id_collect_dic["rbp_id"].append(rbp_id)
            id_collect_dic["method_id"].append(method_id)
            id_collect_dic["data_id"].append(data_id)
            id_collect_dic["set_name"].append(set_name)
            id_collect_dic["path"].append(id2path_dic[set_name])

    f.closed

    assert id_collect_dic["rbp_id"], "nothing read in from --table. Please provide non-empty table in correct format (columns: rbp_id method_id data_id dataset_name)"

    """
    Construct RBPBench batch call.

    """

    batch_call = "rbpbench batch"
    batch_call += " --out %s" % (args.out_folder)
    batch_call += " --genome %s" % (args.in_genome)
    batch_call += " --ext %s" % (args.ext_up_down)
    batch_call += " --motif-db %i" % (args.motif_db)
    if args.fimo_user_ntf_file:
        batch_call += " --fimo-ntf-file %s" % (args.fimo_user_ntf_file)
    batch_call += " --fimo-ntf-mode %i" % (args.fimo_ntf_mode)
    batch_call += " --fimo-pval %s" % (str(args.fimo_pval))
    batch_call += " --cmsearch-bs %s" % (str(args.cmsearch_bs))
    batch_call += " --cmsearch-mode %i" % (args.cmsearch_mode)
    if args.greatest_hits:
        batch_call += " --greatest-hits"
    batch_call += " --bed-score-col %i" % (args.bed_score_col)
    if args.bed_sc_thr is not None:
        batch_call += " --bed-sc-thr %s" % (str(args.bed_sc_thr))
    if args.bed_sc_thr_rev_filter:
        batch_call += " --bed-sc-thr-rev"
    if args.unstranded:
        batch_call += " --unstranded"
    if args.unstranded_ct:
        batch_call += " --unstranded-ct"
    if args.meme_disable_check:
        batch_call += " --meme-no-check"
    if args.meme_no_pgc:
        batch_call += " --meme-no-pgc"

    batch_call += " --kmer-size %i" % (args.kmer_size)
    if args.in_gtf:
        batch_call += " --gtf %s" % (args.in_gtf)
        if args.tr_list:
            batch_call += " --tr-list %s" % (args.tr_list)
        if args.tr_types_list:
            tr_types = (" ").join(args.tr_types_list)
            batch_call += " --tr-types %s" % (tr_types)

    batch_call += " --gtf-feat-min-overlap %s" % (str(args.gtf_feat_min_overlap))
    batch_call += " --gtf-eib-min-overlap %s" % (str(args.gtf_eib_min_overlap))
    batch_call += " --gtf-intron-border-len %i" % (args.gtf_intron_border_len)

    if args.report_header:
        batch_call += " --report-header"
    if args.no_occ_heatmap:
        batch_call += " --no-occ-heatmap"
    if args.disable_heatmap_cluster_olo:
        batch_call += " --disable-heatmap-cluster-olo"

    batch_call += " --fisher-mode %i" % (args.fisher_mode)
    batch_call += " --wrs-mode %i" % (args.wrs_mode)

    if args.regex:
        # Check if given regex is valid.
        assert is_valid_regex(args.regex), "given --regex \"%s\" is not a valid regular expression. Please provide valid expression" % (args.regex)
        # Remove , ; from given regex, to avoid motif_id format conflicts.
        regex = remove_special_chars_from_str(args.regex,
                                              reg_ex="[ ;]")
        
        assert regex, "empty string after removing special chars ( [ ;] ) from --regex. Please provide a valid regex with DNA letters"

        batch_call += " --regex %s" % (regex)
        batch_call += " --regex-search-mode %i" % (args.regex_search_mode)
        batch_call += " --max-motif-dist %i" % (args.max_motif_dist)

    rbp_ids = (" ").join(id_collect_dic["rbp_id"])
    method_ids = (" ").join(id_collect_dic["method_id"])
    data_ids = (" ").join(id_collect_dic["data_id"])
    paths = (" ").join(id_collect_dic["path"])

    batch_call += " --rbp-list %s" % (rbp_ids)
    batch_call += " --method-list %s" % (method_ids)
    batch_call += " --data-list %s" % (data_ids)
    batch_call += " --bed %s" % (paths)

    # GO enrichment analysis.
    if args.run_goa:
        batch_call += " --goa"
        batch_call += " --goa-obo-mode %i" % (args.goa_obo_mode)
        if args.goa_obo_file:
            batch_call += " --goa-obo-file %s" % (args.goa_obo_file)
        if args.goa_gene2go_file:
            batch_call += " --goa-gene2go-file %s" % (args.goa_gene2go_file)
        batch_call += " --goa-pval %s" % (str(args.goa_pval))
        if args.goa_only_cooc:
            batch_call += " --goa-only-cooc"
        if args.goa_bg_gene_list:
            batch_call += " --goa-bg-gene-list %s" % (args.goa_bg_gene_list)
        if args.goa_max_child is not None:
            batch_call += " --goa-max-child %i" % (args.goa_max_child)
        if args.goa_min_depth is not None:
            batch_call += " --goa-min-depth %i" % (args.goa_min_depth)
        if args.goa_filter_purified:
            batch_call += " --goa-filter-purified"

    """
    Execute RBPBench batch call.
    """

    print("")
    print("EXECUTING CALL:\n%s" % (batch_call))
    output = subprocess.getoutput(batch_call)
    print("")
    print("RUN OUTPUT:\n%s" % (output))
    print("")
    print("DONE.")
