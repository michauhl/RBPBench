#!/usr/bin/env python3

import argparse
import pyBigWig
import numpy as np
from rbpbench import benchlib
import statistics
import os

"""
Dependencies:
conda install bioconda::pybigwig

"""

###############################################################################

def setup_argument_parser():
    """Setup argparse parser."""
    help_description = """
    Given a BED file with transcript regions, add gene infos (gene ID, gene name)
    as additional columns.

    """
    # Define argument parser.
    p = argparse.ArgumentParser(add_help=False,
                                prog="bed_tr_context_add_gene_infos.py",
                                description=help_description,
                                formatter_class=argparse.MetavarTypeHelpFormatter)

    # Required arguments.
    p.add_argument("-h", "--help",
                   action="help",
                   help="Print help message")
    p.add_argument("--bed",
                   dest="in_bed",
                   type=str,
                   metavar='str',
                   required = True,
                   help = "Input transcript regions BED file to add gene info columns to")
    p.add_argument("--gtf",
                   dest="in_gtf",
                   type=str,
                   metavar='str',
                   required = True,
                   help = "Input GTF file with genomic annotations to extract gene infos from. Note that only genes on standard chromosomes (1,2,..,X,Y,MT) are currently used")
    p.add_argument("--out",
                   dest="out_bed",
                   type=str,
                   metavar='str',
                   required = True,
                   help = "Output BED file with added gene info columns")
    p.add_argument("--chr-id-style",
                   dest="chr_id_style",
                   type=int,
                   default=1,
                   choices=[1, 2, 3],
                   help="Define to which chromosome ID style to convert chromosome IDs to. 1: do not change chromosome IDs. 2: convert to chr1,chr2,...,chrM style. 3: convert to 1,2,...,MT style (default: 1)")
    return p


################################################################################

if __name__ == '__main__':

    parser = setup_argument_parser()
    args = parser.parse_args()

    assert os.path.exists(args.in_bed), "given --bed BED file \"%s\" not found" % (args.in_bed)
    assert os.path.exists(args.in_gtf), "given --gtf GTF file \"%s\" not found" % (args.in_gtf)

    print("Read in --bed transcript IDs (BED column 1) ... ")
    tr_ids_dic = benchlib.bed_read_chr_ids_dic(args.in_bed)
    assert tr_ids_dic, "no transcript IDs read in from provided --bed file (column 1). Please provide a valid BED file with transcript regions"

    print("Read in gene features from --gtf ... ")
    tr2gid_dic = {}
    gid2gio_dic = benchlib.gtf_read_in_gene_infos(args.in_gtf,
                                                tr2gid_dic=tr2gid_dic,
                                                chr_style=args.chr_id_style,
                                                empty_check=False)

    for tid in tr_ids_dic:
        assert tid in tr2gid_dic, "transcript ID \"%s\" from --bed file not found in --gtf file! Please provide compatible --bed and --gtf files" %(tid)

    print("Add gene infos to BED file ... ")

    BEDOUT = open(args.out_bed, "w")
    
    with open(args.in_bed) as f:
        for line in f:
            cols = line.strip().split("\t")
            assert len(cols) >= 6, "invalid --in BED format. Please provide valid BED file (i.e., >= 6 column format)"
            tr_id = cols[0]

            gene_id = tr2gid_dic[tr_id]
            
            gene_name = gid2gio_dic[gene_id].gene_name

            cols.append(gene_id)
            cols.append(gene_name)

            BEDOUT.write("\t".join(cols) + "\n")

    f.closed
    BEDOUT.close()

