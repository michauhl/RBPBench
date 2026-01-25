#!/usr/bin/env python3

import argparse
import pyBigWig
import numpy as np
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
    Given a BED file with genomic regions, get the average conservation score 
    for each of these regions (given phyloP or phastCons bigWig files via --con).

    """
    # Define argument parser.
    p = argparse.ArgumentParser(add_help=False,
                                prog="bed_gen_context_add_avg_con_scores.py",
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
                   help = "Input genomic regions BED file to extract average conservation scores for")
    p.add_argument("--con",
                   dest="con_sc_file",
                   type=str,
                   metavar='str',
                   required = True,
                   help = "Genomic .bigWig file (phastCons or phyloP) with conservation scores (has to be compatible with --bed)")
    p.add_argument("--out",
                   dest="out_bed",
                   type=str,
                   metavar='str',
                   required = True,
                   help = "Output BED file with added average genomic region conservation scores column")
    return p


################################################################################

if __name__ == '__main__':

    parser = setup_argument_parser()
    args = parser.parse_args()

    assert os.path.exists(args.in_bed), "given --bed BED file \"%s\" not found" % (args.in_bed)
    assert os.path.exists(args.con_sc_file), "given --con bigWig file \"%s\" not found" % (args.con_sc_file)

    """
    Extract conservation scores for given genomic regions.
    
    """

    print("Extract conservation scores for --bed genomic regions ... ")

    # Open conservation scores.
    con_sc_data = pyBigWig.open(args.con_sc_file)
    
    BEDOUT = open(args.out_bed, "w")
    
    with open(args.in_bed) as f:
        for line in f:
            cols = line.strip().split("\t")
            # assert len(cols) >= 6, "invalid --in BED format. Please provide valid BED file (i.e., >= 6 column format)"
            chr_id = cols[0]
            start = int(cols[1])  # 0-based.
            end = int(cols[2])
            reg_len = end - start
            
            avg_score = 0.0
            
            try:
                # Get conservation scores for the region.
                scores = con_sc_data.values(chr_id, start, end, numpy=False)
                c_scores = len(scores)
                assert c_scores == reg_len, "# of extracted scores != region length (%i != %i) for line:\n%s" %(c_scores, reg_len, line)
                # Convert NaN values to 0.0.
                scores = [0.0 if np.isnan(s) else s for s in scores]
                avg_score = statistics.mean(scores) if scores else 0.0

            except RuntimeError:
                print(f"Skipping chunk region {chr_id}:{start}-{end} (coordinates not in bigWig)")

            cols.append(str(avg_score))

            BEDOUT.write("\t".join(cols) + "\n")

    f.closed
    BEDOUT.close()








