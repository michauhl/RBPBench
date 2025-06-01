#!/usr/bin/env python3

import argparse
import os
import re

###############################################################################

def setup_argument_parser():
    """Setup argparse parser."""
    help_description = """
    Given a BED file, extend all regions via --ext, e.g. --ext 20 (different up- 
    and downstream extension is possible too (--ext 30,10).

    """
    # Define argument parser.
    p = argparse.ArgumentParser(add_help=False,
                                prog="bed_extend_regions.py",
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
                   help = "Input BED file with regions to extend")
    p.add_argument("--ext",
                   dest="ext_up_down",
                   type=str,
                   metavar='str',
                   default="0",
                   help="Up- and downstream extension of --in sites in nucleotides (nt). Set e.g. --ext 30 for 30 nt on both sides, or --ext 20,10 for different up- and downstream extension (default: 0)")
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

    # Process extension info.
    ext_parts = args.ext_up_down.split(",")
    c_ext_parts = len(ext_parts)
    ext_up = 0
    ext_down = 0
    if c_ext_parts == 1:
        ext_up = int(ext_parts[0])
        ext_down = int(ext_parts[0])
    elif c_ext_parts == 2:
        ext_up = int(ext_parts[0])
        ext_down = int(ext_parts[1]) 
    else:
        assert False, "invalid --ext argument provided (correct format: --ext 10 OR --ext 20,10)"

    assert boundary_check(ext_up, 0, 100000), "set --ext upstream extension expected to be >= 0 and <= 100000"
    assert boundary_check(ext_down, 0, 100000), "set --ext downstream extension expected to be >= 0 and <= 100000"

    with open(args.in_bed) as f:
        for line in f:
            cols = line.strip().split("\t")
            chr_id = cols[0]
            reg_s = int(cols[1])
            reg_e = int(cols[2])
            reg_id = cols[3]
            sc = cols[4]
            strand = cols[5]

            # Extend.
            new_reg_s = reg_s - ext_up
            new_reg_e = reg_e + ext_down
            if strand == "-":
                new_reg_s = reg_s - ext_down
                new_reg_e = reg_e + ext_up

            # Simple boundary check.
            if new_reg_s < 0:
                new_reg_s = 0

            print("%s\t%i\t%i\t%s\t%s\t%s" % (chr_id, new_reg_s, new_reg_e, reg_id, sc, strand))

    f.closed


################################################################################

