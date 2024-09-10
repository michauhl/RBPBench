#!/usr/bin/env python3

from rbpbench import benchlib
import argparse
import os


###############################################################################

def setup_argument_parser():
    """Setup argparse parser."""
    help_description = """
    Extract exon, intron or both types of region from GTF file to BED file.
    Note that most prominent transcript (MPT) regions are used (unless 
    --tr-list is provided).

    """
    # Define argument parser.
    p = argparse.ArgumentParser(add_help=False,
                                prog="gtf_extract_exon_intron_region_bed.py",
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
                   help = "Input GTF file with genomic annotations to extract exon + intron regions from. Note that only genes on standard chromosomes (1,2,..,X,Y,MT) are currently used")
    p.add_argument("--out",
                   dest="out_bed",
                   type=str,
                   metavar='str',
                   required=True,
                   help="Output BED file to store exon and/or intron regions in")
    p.add_argument("--type-mode",
                   dest="type_mode",
                   type=int,
                   default=1,
                   choices=[1, 2, 3],
                   help="Define what regions to output to --out BED 1: exon + intron regions. 2: exon regions only. 3: intron regions only (default: 1)")
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
                   help = "Supply file with transcript IDs (one ID per row) to define which transcripts to extract exon + intron regions from in --gtf. NOTE that this overwrites most prominent transcript selection")

    return p


################################################################################

if __name__ == '__main__':

    parser = setup_argument_parser()
    args = parser.parse_args()

    assert os.path.exists(args.in_gtf), "--gtf file \"%s\" not found" % (args.in_gtf)

    print("Read in gene features from --gtf ... ")
    tr2gid_dic = {}
    tr_types_dic = {}  # Store transcript biotypes in GTF file.
    gid2gio_dic = benchlib.gtf_read_in_gene_infos(args.in_gtf,
                                                    tr2gid_dic=tr2gid_dic,
                                                    chr_style=args.chr_id_style,
                                                    empty_check=False)

    assert gid2gio_dic, "no gene infos read in from --gtf. Please provide a valid/compatible GTF file (e.g. from Ensembl or ENCODE)"

    c_gene_infos = len(gid2gio_dic)
    print("# gene features read in from --gtf:", c_gene_infos)

    tr_ids_dic = {}
    if args.tr_list:
        tr_ids_dic = benchlib.read_ids_into_dic(args.tr_list,
                                                check_dic=False)
        assert tr_ids_dic, "no IDs read in from provided --tr-list file. Please provide a valid IDs file (one ID per row)"
        for tr_id in tr_ids_dic:
            assert tr_id in tr2gid_dic, "transcript ID \"%s\" from provided --tr-list file does not appear in --gtf file. Please provide compatible files" %(tr_id)
            tr_ids_dic[tr_id] = tr2gid_dic[tr_id]
        print("# of transcript IDs (read in from --tr-list):", len(tr_ids_dic))
    else:
        # Get most prominent transcripts from gene infos.
        tr_ids_dic = benchlib.select_mpts_from_gene_infos(gid2gio_dic,
                                basic_tag=False,  # do not be strict (only_tsl=False too).
                                ensembl_canonical_tag=False,
                                prior_basic_tag=True,  # Prioritize basic tag transcript.
                                only_tsl=False)
        assert tr_ids_dic, "most prominent transcript selection from gene infos failed. Please contact developers"
        print("# of transcript IDs (most prominent transcripts):", len(tr_ids_dic))

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
                                                        tr_ids_dic=tr_ids_dic,
                                                        correct_min_ex_order=correct_min_ex_order,
                                                        chr_style=args.chr_id_style,
                                                        empty_check=False)

    assert tid2tio_dic, "no transcript infos read in from --gtf. Please provide a valid/compatible GTF file (e.g. from Ensembl or ENCODE)"

    # (in)sanity checks.
    for tr_id in tr_ids_dic:
        assert tr_id in tid2tio_dic, "transcript ID %s not in tid2tio_dic"
    for tr_id in tid2tio_dic:
        assert tr_id in tr_ids_dic, "transcript ID %s not in tr_ids_dic"

    c_tr_infos = len(tid2tio_dic)
    print("# transcript features read in from --gtf:", c_tr_infos)

    print("Output exon and/or intron annotations to --out BED ... ")

    benchlib.output_transcript_info_intron_exon_to_bed(tid2tio_dic, args.out_bed,
                                        output_mode=args.type_mode,
                                        report_counts=True,
                                        add_tr_id=True,
                                        empty_check=False)

    print("Done.")


################################################################################
