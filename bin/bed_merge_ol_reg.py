#!/usr/bin/env python3

from rbpbench import benchlib
import argparse
import os


###############################################################################

def setup_argument_parser():
    """Setup argparse parser."""
    help_description = """
    Given --in BED file, merge bookend or overlapping regions and output merged 
    regions to --out BED file. Optionally, extend --in regions by --ext
    nucleotides up- and/or downstream before merging. Note that new unique BED
    column 4 region IDs will be generated, unless disabled.

    """
    # Define argument parser.
    p = argparse.ArgumentParser(add_help=False,
                                prog="bed_merge_ol_reg.py",
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
                   required=True,
                   help="Input BED file to merge overlapping or bookend regions in")
    p.add_argument("--out",
                   dest="out_bed",
                   type=str,
                   metavar='str',
                   required=True,
                   help="Output BED file to store merged regions in")
    p.add_argument("--ext",
                   dest="ext_up_down",
                   type=str,
                   metavar='str',
                   default="0",
                   help="Up- and downstream extension of --in sites in nucleotides (nt). Set e.g. --ext 30 for 30 nt on both sides, or --ext 20,10 for different up- and downstream extension (default: 0)")
    p.add_argument("--core-reg-id",
                   dest="core_reg_id",
                   type=str,
                   metavar='str',
                   default = "s",
                   help="New core BED column 4 region ID, so by default new region IDs will have format s1, s2, s3, ... (default: s)")
    p.add_argument("--disable-new-ids",
                   dest="disable_new_ids",
                   default = False,
                   action = "store_true",
                   help = "Disable generation of new BED column 4 region IDs. NOTE that if --in IDs are same, --disable-new-ids messes up results (default: False)")

    return p


################################################################################

if __name__ == '__main__':

    parser = setup_argument_parser()
    args = parser.parse_args()

    assert os.path.exists(args.in_bed), "--in BED file \"%s\" not found" % (args.in_bed)


    tmp_bed = args.out_bed + ".tmp"

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

    # Preprocess --in BED.
    print("Preprocess --in regions ... ")

    reg2pol_dic = {}

    new_reg_ids = True
    if args.disable_new_ids:
        new_reg_ids = False

    reg_stats_dic = benchlib.bed_filter_extend_bed(args.in_bed, tmp_bed,
                                          ext_up=ext_up,
                                          ext_down=ext_down,
                                          remove_dupl=False,
                                          reg2pol_dic=reg2pol_dic,
                                          core_reg_id=args.core_reg_id,
                                          new_reg_ids=new_reg_ids)

    print("# --in regions pre-filtering:  ", reg_stats_dic["c_in"])
    print("# --in regions post-filtering: ", reg_stats_dic["c_out"])
    print("# regions with invalid chr_id: ", reg_stats_dic["c_chr_filter"])
    print("# duplicated regions removed:  ", reg_stats_dic["c_dupl_filter"])
    print("# regions filtered by score:   ", reg_stats_dic["c_sc_thr"])

    assert reg_stats_dic["c_out"], "no --in BED sites remain after chromosome ID (or optionally score) filtering. If caused by invalid chr_id filtering, make sure chromosome IDs in --genome FASTA and --in BED files are compatible (i.e., \"chr1\" vs. \"1\" notation). If --in regions are on transcripts, use rbpbench searchrna"

    print("Merge overlapping and bookend regions ... ")
    benchlib.bed_get_effective_reg_bed(tmp_bed, args.out_bed, reg2pol_dic)

    c_merged_reg = benchlib.count_lines_in_file(args.out_bed)

    # Delete temporary BED file.
    os.remove(tmp_bed)

    print("# merged regions:              ", c_merged_reg)
    print("Merged regions written to --out BED file:\n", args.out_bed)
               