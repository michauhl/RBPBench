#!/usr/bin/env python3

import argparse
import pyBigWig
import numpy as np
import matplotlib.pyplot as plt
import os

"""
Dependencies:
conda install bioconda::pybigwig

"""

###############################################################################

def setup_argument_parser():
    """Setup argparse parser."""
    help_description = """
    Given a set of genomic regions of interest (--in BED) and a file with phylogenomic 
    conservation scores (phastCons, phyloP, --con bigWig), extract conservation scores 
    for these regions.

    """
    # Define argument parser.
    p = argparse.ArgumentParser(add_help=False,
                                prog="get_genomic_conservation_scores.py",
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
                   help = "Input genomic regions of interest in BED format to extract conservation scores for")
    p.add_argument("--con",
                   dest="con_sc_file",
                   type=str,
                   metavar='str',
                   required = True,
                   help = "Genomic .bigWig file (phastCons or phyloP) with conservation scores")
    p.add_argument("--out",
                   dest="out_file",
                   type=str,
                   metavar='str',
                   required = True,
                   help = "Output table file with conservation scores")
    p.add_argument("--use-regions",
                   dest="use_regions",
                   default = False,
                   action = "store_true",
                   help="Use genomic regions as --bed / --control-bed input region IDs instead of BED col4 IDs (default: False)")
    p.add_argument("--no-id-check",
                   dest="no_id_check",
                   default = False,
                   action = "store_true",
                   help="Do not check region IDs, instead overwriting existing regions if they have identical IDs (default: False)")
    p.add_argument("--keep-nan",
                   dest="keep_nan",
                   default = False,
                   action = "store_true",
                   help="Do not convert nan conservation values to zero (0.0) (default: False)")
    return p

################################################################################

def bed_read_in_regions(in_bed,
                        no_id_check=False,
                        use_regions_as_ids=False):
    """
    Read in BED regions into a dictionary.

    """

    bed_regions_dic = {}

    with open(in_bed) as f:
        for line in f:
            cols = line.strip().split("\t")
            chr_id = cols[0]
            start = int(cols[1])
            end = int(cols[2])
            reg_id = cols[3]
            pol = cols[5]
            
            if use_regions_as_ids:
                reg_id = chr_id + ":" + str(start) + "-" + str(end) + "," + pol

            if no_id_check:
                if reg_id in bed_regions_dic:
                    # Overwrite existing region.
                    print("WARNING: region ID \"%s\" already found in \"%s\". Overwriting existing region .." % (reg_id, in_bed))
            else:
                assert reg_id not in bed_regions_dic, "region ID \"%s\" already found in \"%s\". Please provide BED file with unique col4 IDs or set --use-regions" % (reg_id, in_bed)
            bed_regions_dic[reg_id] = [chr_id, start, end, pol]

    f.closed
    return bed_regions_dic


################################################################################

if __name__ == '__main__':

    parser = setup_argument_parser()
    args = parser.parse_args()

    assert os.path.exists(args.in_bed), "given --in BED file \"%s\" not found" % (args.in_bed)
    assert os.path.exists(args.con_sc_file), "given --con bigWig file \"%s\" not found" % (args.con_sc_file)

    print("Read in --in BED regions  ... ")

    # Read in regions of interest.
    in_regions_dic = bed_read_in_regions(args.in_bed,
                                         no_id_check=args.no_id_check,
                                         use_regions_as_ids=args.use_regions)
    
    print("# --in BED regions: %i" % (len(in_regions_dic)))

    """
    Get conservation scores for --in BED region positions.
    
    """

    print("Read in conservation scores ... ")

    con_sc_data = pyBigWig.open(args.con_sc_file)

    RESOUT = open(args.out_file, "w")

    for reg_id, reg_info in in_regions_dic.items():

        chr_id = reg_info[0]
        start = reg_info[1]
        end = reg_info[2]

        try:
            # Get conservation scores for the region.
            scores = con_sc_data.values(chr_id, start, end, numpy=False)
            # Convert NaN values to 0.0.
            if not args.keep_nan:
                scores = [0.0 if np.isnan(s) else s for s in scores]

            scores_str = ",".join([str(s) for s in scores])
            # Write out conservation scores.
            RESOUT.write("%s\t%s\t%s\t%s\t%s\n" % (reg_id, chr_id, start, end, scores_str))

        except RuntimeError:
            print(f"Skipping --bed region {chr_id}:{start}-{end} (coordinates not in bigWig)")


    RESOUT.close()
    
    print("Conservation scores written to:\n%s" %(args.out_file))
    
