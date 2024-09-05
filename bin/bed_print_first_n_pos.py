#!/usr/bin/env python3

import argparse
import os
from rbpbench import benchlib


###############################################################################

def setup_argument_parser():
    """Setup argparse parser."""
    help_description = """
    Print first n (set via --ext) positions of each region from --in BED file. 
    For minus strand regions, the last n positions are printed.

    """
    # Define argument parser.
    p = argparse.ArgumentParser(add_help=False,
                                prog="bed_print_first_n_pos.py",
                                description=help_description,
                                formatter_class=argparse.MetavarTypeHelpFormatter)

    # Required arguments.
    p.add_argument("-h", "--help",
                   action="help",
                   help="Print help message")
    p.add_argument("--in",
                   dest="in_bed",
                   type=str,
                   metavar='str',
                   required = True,
                   help = "Input BED file with regions to extract first n positions from")
    p.add_argument("--ext",
                   dest="ext",
                   type=int,
                   metavar='int',
                   required=True,
                   help="Print first --ext positions of --in BED regions")
    return p


################################################################################

if __name__ == '__main__':

    parser = setup_argument_parser()
    args = parser.parse_args()

    assert os.path.exists(args.in_bed), "--in BED file \"%s\" not found" % (args.in_bed)
    assert benchlib.boundary_check(args.ext, 1, 1000000), "set --ext expected to be >= 1 and <= 1000000"

    with open(args.in_bed) as f:
        for line in f:
            cols = line.strip().split("\t")
            chr_id = cols[0]
            reg_s = int(cols[1])
            reg_e = int(cols[2])
            reg_id = cols[3]
            sc = cols[4]
            strand = cols[5]

            new_reg_s = reg_s
            new_reg_e = reg_s + args.ext

            if strand == "-":
                new_reg_s = reg_e - args.ext
                new_reg_e = reg_e
                if new_reg_s < 0:
                    new_reg_s = 0

            print("%s\t%i\t%i\t%s\t%s\t%s" % (chr_id, new_reg_s, new_reg_e, reg_id, sc, strand))

    f.closed


################################################################################
