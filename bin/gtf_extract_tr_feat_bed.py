#!/usr/bin/env python3

import argparse
import os
from rbpbench import benchlib


###############################################################################

def setup_argument_parser():
    """Setup argparse parser."""
    help_description = """
    Extract --feat transcript feature (e.g. --feat stop_codon) from --gtf GTF file 
    and store in --out BED file. If start_codon or stop_codon are used, length 
    check of 3 nt is performed and filtered. Note that only genes on standard
    chromosomes (1,2,..,X,Y,MT) are currently used.

    """
    # Define argument parser.
    p = argparse.ArgumentParser(add_help=False,
                                prog="gtf_extract_tr_feat_bed.py",
                                description=help_description,
                                formatter_class=argparse.MetavarTypeHelpFormatter)

    # Required arguments.
    p.add_argument("-h", "--help",
                   action="help",
                   help="Print help message")
    p.add_argument("--feat",
                   dest="feat_in",
                   type=str,
                   metavar='str',
                   required = True,
                   help = "Supply transcript related feature to be extracted from GTF file and stored as BED file (e.g. transcript, exon, CDS, start_codon, stop_codon)")
    p.add_argument("--gtf",
                   dest="in_gtf",
                   type=str,
                   metavar='str',
                   required = True,
                   help = "Input GTF file with genomic annotations to extract exon + intron regions from. Note that only genes on standard chromosomes (1,2,..,X,Y,MT) are currently used")
    p.add_argument("--out",
                   dest="out_bed",
                   type=str,
                   metavar='str',
                   required=True,
                   help="Output BED file to store transcript feature regions in")
    p.add_argument("--chr-id-style",
                   dest="chr_id_style",
                   type=int,
                   default=1,
                   choices=[1, 2, 3],
                   help="Define to which chromosome ID style to convert chromosome IDs to. 1: do not change chromosome IDs. 2: convert to chr1,chr2,...,chrM style. 3: convert to 1,2,...,MT style (default: 1)")
    p.add_argument("--tr-list",
                   dest="tr_list",
                   type=str,
                   metavar='str',
                   default = False,
                   help = "Supply file with transcript IDs (one ID per row) to define which transcripts to extract feature regions from --gtf. By default, all transcripts are used")
    p.add_argument("--uniq-reg",
                   dest="unique_regions",
                   default = False,
                   action = "store_true",
                   help = "Only output unique feature regions. I.e., if genomic region already encountered, do not output again. E.g. stop_codon can be the same region for different transcripts (default: False)")

    return p


################################################################################

if __name__ == '__main__':

    parser = setup_argument_parser()
    args = parser.parse_args()

    assert os.path.exists(args.in_gtf), "--gtf GTF file \"%s\" not found" % (args.in_gtf)

    tr_ids_dic = {}
    if args.tr_list:
        assert os.path.exists(args.tr_list), "given --tr-list file \"%s\" not found" % (args.tr_list)
        tr_ids_dic = benchlib.read_ids_into_dic(args.tr_list,
                                                check_dic=False)
        assert tr_ids_dic, "no IDs read in from provided --tr-list file. Please provide a valid IDs file (one ID per row)"
        print("# of transcript IDs (read in from --tr-list): ", len(tr_ids_dic))

    tr_feat = args.feat_in

    benchlib.gtf_tr_feat_to_bed(args.in_gtf, args.out_bed, tr_feat,
                                chr_style=args.chr_id_style,
                                uniq_reg=args.unique_regions,
                                tr_ids_dic=tr_ids_dic)


################################################################################
