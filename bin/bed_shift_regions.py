#!/usr/bin/env python3

import argparse
import os
import re

###############################################################################

def setup_argument_parser():
    """Setup argparse parser."""
    help_description = """
    Given a BED file, shift all regions by a given number of nucleotides.

    """
    # Define argument parser.
    p = argparse.ArgumentParser(add_help=False,
                                prog="bed_shift_regions.py",
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
                   help = "Input BED file with regions to shift")
    p.add_argument("--num",
                   dest="shift_num",
                   type=int,
                   metavar='int',
                   required=True,
                   help="Number to shift regions by. E.g. --num 10 will shift all regions by 10 nucleotides downstream, --num -10 will shift all regions by 10 nucleotides upstream")
    return p


################################################################################

def boundary_check(value, min, max):
    """
    Return if value is within integer/float boundaries.

    >>> value = 0.5
    >>> min = 1E-9
    >>> max = 1.0
    >>> boundary_check(value, min, max)
    True
    >>> value = 0.0
    >>> min = 1E-9
    >>> max = 1.0
    >>> boundary_check(value, min, max)
    False

    """
    if value >= min and value <= max:
        return True
    else:
        return False


################################################################################

if __name__ == '__main__':

    parser = setup_argument_parser()
    args = parser.parse_args()

    assert os.path.exists(args.in_bed), "--in BED file \"%s\" not found" % (args.in_bed)
    assert boundary_check(args.shift_num, -100000, 100000), "set --num expected to be >= -100000 and <= 100000"

    with open(args.in_bed) as f:
        for line in f:
            cols = line.strip().split("\t")
            chr_id = cols[0]
            reg_s = cols[1]
            reg_e = cols[2]
            reg_id = cols[3]
            sc = cols[4]
            strand = cols[5]

            new_reg_s = reg_s
            new_reg_e = reg_e

            if strand == "+":
                new_reg_s = int(reg_s) + args.shift_num
                new_reg_e = int(reg_e) + args.shift_num
            else:
                new_reg_s = int(reg_s) - args.shift_num
                new_reg_e = int(reg_e) - args.shift_num

            if new_reg_s < 0:
                new_reg_s = 0

            print("%s\t%i\t%i\t%s\t%s\t%s" % (chr_id, new_reg_s, new_reg_e, reg_id, sc, strand))

    f.closed


################################################################################
