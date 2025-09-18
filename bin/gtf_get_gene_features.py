#!/usr/bin/env python3

import argparse
import pyBigWig
from rbpbench import benchlib
import numpy as np
import statistics
import os

"""
Get gene features from GTF file, optionally also including conservation 
data.

"""

###############################################################################

def setup_argument_parser():
    """Setup argparse parser."""
    help_description = """
    Get gene features from GTF file, optionally also including conservation 
    data.

    """
    # Define argument parser.
    p = argparse.ArgumentParser(add_help=False,
                                prog="gtf_get_gene_features.py",
                                description=help_description,
                                formatter_class=argparse.MetavarTypeHelpFormatter)

    # Required arguments.
    p.add_argument("-h", "--help",
                   action="help",
                   help="Print help message")
    p.add_argument("--gtf",
                   dest="in_gtf",
                   type=str,
                   metavar='str',
                   required = True,
                   help = "Input GTF file compatible with transcriptome FASTA file used for read mapping")
    # p.add_argument("--genome",
    #                dest="in_genome",
    #                type=str,
    #                metavar='str',
    #                help = "Genomic sequences file (currently supported formats: FASTA)")
    p.add_argument("--phastcons",
                   dest="pc_bw",
                   type=str,
                   metavar='str',
                   help = "Genomic .bigWig file with phastCons conservation scores")
    p.add_argument("--phylop",
                   dest="pp_bw",
                   type=str,
                   metavar='str',
                   help = "Genomic .bigWig file with phyloP conservation scores")
    p.add_argument("--chr-style",
                   dest="chr_style",
                   type=int,
                   default=1,
                   choices=[1, 2],
                   help="Defines chromosome style in supplied bigWig file, to ensure compatibility with --gtf file 1: chr1,chr2,.. 2: 1,2,.. (default: 1)")
    p.add_argument("--out",
                   dest="out_table",
                   type=str,
                   metavar='str',
                   required = True,
                   help = "Output gene features table file")
    # p.add_argument("--skip-chr-check",
    #                dest="skip_chr_check",
    #                default = False,
    #                action = "store_true",
    #                help="Skip chromosome ID check, i.e. do not filter out entries with non-standard chromosome IDs (default: False)")
    return p


################################################################################

if __name__ == '__main__':

    parser = setup_argument_parser()
    args = parser.parse_args()
    
    # assert os.path.exists(args.in_fasta), "given --fasta file \"%s\" not found" % (args.in_fasta)
    assert os.path.exists(args.in_gtf), "given --gtf file \"%s\" not found" % (args.in_gtf)


    """
    If --genome given, get chrosome IDs style from --genome FASTA file.

    """

    chr_ids_dic = None
    chr_style = args.chr_style  # 1: change to chr1, chr2 to ensure bigWig compatibility.

    # if args.in_genome:

    #     assert args.pc_bw or args.pp_bw, "if --genome given, --phastcons or --phylop must be given as well"
    #     assert os.path.exists(args.in_genome), "given --genome file \"%s\" not found" % (args.in_genome)

    #     print("Get --genome FASTA headers ... ")
    #     chr_ids_dic = nanolib.get_fasta_headers(args.in_genome)

    #     print("Guess chromosome ID style (based on --genome FASTA headers) ... ")
    #     chr_style = nanolib.guess_chr_id_style(chr_ids_dic)
        

    """
    Get gene + transcript infos.

    """

    print("Get gene infos from --gtf ... ")

    # Get gene infos.
    print("Read in gene features from --gtf ... ")
    tr2gid_dic = {}
    tr_types_dic = {}
    add_version_nr = False

    gid2gio_dic = benchlib.gtf_read_in_gene_infos(args.in_gtf,
                                            tr2gid_dic=tr2gid_dic,
                                            tr_types_dic=tr_types_dic,
                                            check_chr_ids_dic=chr_ids_dic,
                                            chr_style=chr_style,
                                            # add_version_nr=add_version_nr,
                                            # skip_chr_check=args.skip_chr_check,
                                            empty_check=False)

    assert gid2gio_dic, "no gene infos read in from --gtf. Please provide a valid/compatible GTF file (e.g. from Ensembl or ENCODE)"
    c_gene_infos = len(gid2gio_dic)
    print("# gene features read in from --gtf:", c_gene_infos)

    # Check exon order (return True if minus strand exon 1 is most downstream, not most upstream, which is the correct way).
    print("Check minus-strand exon order in --gtf ... ")
    correct_min_ex_order = benchlib.gtf_check_exon_order(args.in_gtf)

    if correct_min_ex_order:
        print("Correct order encountered ... ")
    else:
        print("Reverse order encountered ... ")

    # Get transcript infos.
    print("Read in transcript infos from --gtf ... ")
    tid2tio_dic = benchlib.gtf_read_in_transcript_infos(args.in_gtf, 
                                                tr_ids_dic=tr2gid_dic,
                                                correct_min_ex_order=correct_min_ex_order,
                                                chr_style=chr_style,
                                                # add_version_nr=add_version_nr,
                                                # skip_chr_check=args.skip_chr_check,
                                                empty_check=False)

    assert tid2tio_dic, "no transcript infos read in from --gtf. Please provide a valid/compatible GTF file (e.g. from Ensembl or ENCODE)"
    print("# transcript features read in from --gtf:", len(tid2tio_dic))

    # Get mRNA transcripts, with 5'UTR,CDS,3'UTR lengths list.
    print("Get mRNA region lengths ... ")
    tid2regl_dic = benchlib.get_mrna_region_lengths(tid2tio_dic)

    c_mrnas = len(tid2regl_dic)
    print("# of mRNA transcripts: %i" %(c_mrnas))

    # Get most prominent transcripts from gene infos.
    mpt_ids_dic = benchlib.select_mpts_from_gene_infos(gid2gio_dic,
                            basic_tag=False,  # do not be strict (only_tsl=False too).
                            ensembl_canonical_tag=False,
                            prior_basic_tag=True,  # Prioritize basic tag transcript.
                            prior_mane_select=True,  # mane select if set trumps all.
                            prior_lncrna_primary_tag=True,
                            only_tsl=False)
    assert mpt_ids_dic, "most prominent transcript selection from gene infos failed. Please contact developers"
    print("# of transcript IDs (most prominent transcripts): ", len(mpt_ids_dic))

    gid2mpt_dic = {}
    for mpt in mpt_ids_dic:
        gid = mpt_ids_dic[mpt]
        gid2mpt_dic[gid] = mpt

    pc_bw = False
    pp_bw = False
    if args.pc_bw:
        assert os.path.exists(args.pc_bw), "given --phastcons file \"%s\" not found" % (args.pc_bw)
        pc_bw = pyBigWig.open(args.pc_bw)
    if args.pp_bw:
        assert os.path.exists(args.pp_bw), "given --phylop file \"%s\" not found" % (args.pp_bw)
        pp_bw = pyBigWig.open(args.pp_bw)


    with open(args.out_table, "w") as f:

        # Write header.
        f.write("gene_id\tgene_name\tgene_biotype\tgene_len\tc_tr_ids\t"
        "mpt_id\tmpt_tr_length\tmpt_exon_c\tmpt_avg_exon_len\tmpt_median_exon_len\tmpt_std_exon_len\tmpt_max_exon_len\tmpt_min_exon_len\t"
        "mpt_intron_c\tmpt_total_intron_len\tmpt_avg_intron_len\tmpt_median_intron_len\tmpt_std_intron_len\tmpt_max_intron_len\tmpt_min_intron_len\t"
        "mpt_utr5_len\tmpt_cds_len\tmpt_utr3_len\tpc_avg_score\tpc_std_score\tpp_avg_score\tpp_std_score\n")

        for gid in gid2gio_dic:
            # Get gene info.
            gio = gid2gio_dic[gid]
            gene_name = gio.gene_name
            gene_biotype = gio.gene_biotype
            chr_id = gio.chr_id
            gene_len = gio.gene_e - gio.gene_s + 1
            # Get transcript IDs.
            tr_ids = gio.tr_ids
            c_tr_ids = len(tr_ids)
            # Get transcript biotypes.
            tr_biotypes = gio.tr_biotypes
            # Get transcript lengths (combined exon lengths).
            tr_lengths = gio.tr_lengths
            # MPT.
            mpt_id = gid2mpt_dic[gid]
            mpt_exon_coords = tid2tio_dic[mpt_id].exon_coords
            mpt_exon_c = tid2tio_dic[mpt_id].exon_c
            mpt_intron_c = mpt_exon_c - 1
            mpt_intron_coords = tid2tio_dic[mpt_id].intron_coords
            mpt_tr_length = tid2tio_dic[mpt_id].tr_length

            mpt_exon_lengths = []
            mpt_intron_lengths = []
            for exon_se in mpt_exon_coords:
                exon_len = exon_se[1] - exon_se[0] + 1
                mpt_exon_lengths.append(exon_len)
            for intron_se in mpt_intron_coords:
                intron_len = intron_se[1] - intron_se[0] + 1
                mpt_intron_lengths.append(intron_len)

            mpt_avg_exon_len = 0.0
            mpt_std_exon_len = 0.0
            mpt_median_exon_len = 0.0
            mpt_total_exon_len = 0
            mpt_max_exon_len = 0
            mpt_min_exon_len = 0
            mpt_avg_intron_len = 0.0
            mpt_std_intron_len = 0.0
            mpt_median_intron_len = 0.0
            mpt_total_intron_len = 0
            mpt_max_intron_len = 0
            mpt_min_intron_len = 0
            if mpt_exon_lengths:
                mpt_avg_exon_len = statistics.mean(mpt_exon_lengths)
                mpt_std_exon_len = statistics.stdev(mpt_exon_lengths) if len(mpt_exon_lengths) > 1 else 0.0
                mpt_median_exon_len = statistics.median(mpt_exon_lengths)
                mpt_total_exon_len = sum(mpt_exon_lengths)
                mpt_max_exon_len = max(mpt_exon_lengths)
                mpt_min_exon_len = min(mpt_exon_lengths)
            if mpt_intron_lengths:
                mpt_avg_intron_len = statistics.mean(mpt_intron_lengths)
                mpt_std_intron_len = statistics.stdev(mpt_intron_lengths) if len(mpt_intron_lengths) > 1 else 0.0
                mpt_median_intron_len = statistics.median(mpt_intron_lengths)
                mpt_total_intron_len = sum(mpt_intron_lengths)
                mpt_max_intron_len = max(mpt_intron_lengths)
                mpt_min_intron_len = min(mpt_intron_lengths)
            
            # Get mRNA region lengths.
            mpt_utr5_len = 0
            mpt_cds_len = 0
            mpt_utr3_len = 0
            if mpt_id in tid2regl_dic:
                mpt_utr5_len, mpt_cds_len, mpt_utr3_len = tid2regl_dic[mpt_id]

            pc_exon_sc = []
            pp_exon_sc = []

            for exon_se in mpt_exon_coords:

                # Get exon start and end coordinates.
                start = exon_se[0] - 1
                end = exon_se[1]

                # phastCons.
                if pc_bw:
                    try:
                        # Get conservation scores for the region.
                        scores = pc_bw.values(chr_id, start, end, numpy=False)
                        # Convert NaN values to 0.0.
                        scores = [0.0 if np.isnan(s) else s for s in scores]
                        for sc in scores:
                            pc_exon_sc.append(sc)
                    except RuntimeError:
                        print(f"Skipping exon region {chr_id}:{start}-{end} (coordinates not in phastCons bigWig)")

                # phyloP.
                if pp_bw:
                    try:
                        # Get conservation scores for the region.
                        scores = pp_bw.values(chr_id, start, end, numpy=False)
                        # Convert NaN values to 0.0.
                        scores = [0.0 if np.isnan(s) else s for s in scores]
                        for sc in scores:
                            pp_exon_sc.append(sc)
                    except RuntimeError:
                        print(f"Skipping exon region {chr_id}:{start}-{end} (coordinates not in phyloP bigWig)")

            pc_avg_score = statistics.mean(pc_exon_sc) if pc_exon_sc else 0.0
            pc_std_score = statistics.stdev(pc_exon_sc) if len(pc_exon_sc) > 1 else 0.0
            pp_avg_score = statistics.mean(pp_exon_sc) if pp_exon_sc else 0.0
            pp_std_score = statistics.stdev(pp_exon_sc) if len(pp_exon_sc) > 1 else 0.0

            # Write gene features to table.
            f.write("%s\t%s\t%s\t%i\t%i\t%s\t%i\t%i\t%.2f\t%.2f\t%.2f\t%i\t%i\t%i\t%i\t%.2f\t%.2f\t%.2f\t%i\t%i\t" 
            "%i\t%i\t%i\t%.2f\t%.2f\t%.2f\t%.2f\n" % (gid, gene_name, gene_biotype, gene_len, c_tr_ids,
            mpt_id, mpt_tr_length, mpt_exon_c, mpt_avg_exon_len, mpt_median_exon_len, mpt_std_exon_len, mpt_max_exon_len, mpt_min_exon_len, 
            mpt_intron_c, mpt_total_intron_len, mpt_avg_intron_len, mpt_median_intron_len, mpt_std_intron_len, mpt_max_intron_len, mpt_min_intron_len,
            mpt_utr5_len, mpt_cds_len, mpt_utr3_len,
            pc_avg_score, pc_std_score, pp_avg_score, pp_std_score)) #24

    pc_bw.close()
    pp_bw.close()

    print("Output written to: %s" % (args.out_table))
