#!/usr/bin/env python3

import argparse
import os
from rbpbench import benchlib

"""
Example call:
python gtf_extract_transcript_data.py --gtf /path/to/Homo_sapiens.GRCh38.112.gtf.gz --genome /path/to/hg38.fa --out wanted_genes_out --gene-list wanted_genes.txt

"""


###############################################################################

def setup_argument_parser():
    """Setup argparse parser."""
    help_description = """
    Given a list of gene IDs (--in) and a --gtf file, output (--out) list of
    transcript IDs belonging to input gene IDs from GTF file.

    """
    # Define argument parser.
    p = argparse.ArgumentParser(add_help=False,
                                prog="gtf_get_gene_transcripts.py",
                                description=help_description,
                                formatter_class=argparse.MetavarTypeHelpFormatter)

    # Required arguments.
    p.add_argument("-h", "--help",
                   action="help",
                   help="Print help message")
    p.add_argument("--in",
                   dest="gene_list",
                   type=str,
                   metavar='str',
                   required = True,
                   help = "Supply file with gene IDs (one ID per row) to define for which genes to extract transcript IDs from --gtf")
    p.add_argument("--gtf",
                   dest="in_gtf",
                   type=str,
                   metavar='str',
                   required = True,
                   help = "Input GTF file with genomic annotations to extract transcript data from. Note that only genes on standard chromosomes (1,2,..,X,Y,MT) are currently used")
    p.add_argument("--out",
                   dest="out_tr_list",
                   type=str,
                   metavar='str',
                   required=True,
                   help="Output transcript IDs table file (including gene IDs, gene names and gene biotypes)")
    p.add_argument("--ignore-version-numbers",
                   dest="ignore_version_numbers",
                   default = False,
                   action = "store_true",
                   help = "Set to ignore ID version numbers in --gtf file, i.e., read in gene and transcript IDs without version numbers. This has to be set if input IDs have no version number but GTF file has (default: False)")
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

    assert os.path.exists(args.in_gtf), "--gtf file \"%s\" not found" % (args.in_gtf)
    assert os.path.exists(args.gene_list), "--gene-list file \"%s\" not found" % (args.gene_list)

    print("Read in gene IDs from --in file ...")
    gene_ids_dic = benchlib.read_ids_into_dic(args.gene_list,
                                              check_dic=False)
    
    assert gene_ids_dic, "no IDs read in from provided --gene-list file. Please provide a valid IDs file (one ID per row)"

    tr2gid_dic = {}
    print("Read in gene infos from --gtf file ... ")
    gid2gio_dic = benchlib.gtf_read_in_gene_infos(args.in_gtf,
                                                  tr2gid_dic=tr2gid_dic,
                                                  chr_style=args.chr_id_style,
                                                  gene_ids_dic=gene_ids_dic,
                                                  remove_version_numbers=args.ignore_version_numbers,
                                                  empty_check=False)

    assert gid2gio_dic, "no gene infos read in from --gtf. Please provide a valid/compatible GTF file (e.g. from Ensembl or ENCODE)"

    for gid in gene_ids_dic:
        if gid not in gid2gio_dic:
            print("WARNING: gene ID \"%s\" not found in GTF file!" % (gid))

    # Collect transcript IDs.
    tr_ids_dic = {}
    for gid in gid2gio_dic:
        tr_ids = gid2gio_dic[gid].tr_ids
        for idx, tid in enumerate(tr_ids):
            tr_ids_dic[tid] = 1
    
    print("Found %d transcript IDs in GTF file belonging to %d genes" % (len(tr_ids_dic), len(gid2gio_dic)))

    # Get most prominent transcripts from gene infos.
    print("Select most prominent transcript (MPT) for each gene ... ")
    mpt_ids_dic = benchlib.select_mpts_from_gene_infos(gid2gio_dic,
                            basic_tag=False,  # do not be strict (only_tsl=False too).
                            ensembl_canonical_tag=False,
                            prior_basic_tag=True,  # Prioritize basic tag transcript.
                            prior_mane_select=True,  # mane select if set trumps all.
                            prior_lncrna_primary_tag=True,  # for lncRNA genes prioritize gencode primary tagged transcripts (mane select still better but should not occur together for lncRNAs).
                            only_tsl=False,
                            gene_ids_dic=gene_ids_dic)

    assert mpt_ids_dic, "most prominent transcript selection from gene infos failed. Please contact developers"
    print("# of transcript IDs (most prominent transcripts): ", len(mpt_ids_dic))


    # Check exon order (return True if minus strand exon 1 is most downstream, not most upstream, which is the correct way).
    print("Check minus-strand exon order in --gtf ... ")
    correct_min_ex_order = benchlib.gtf_check_exon_order(args.in_gtf)
    if correct_min_ex_order:
        print("Correct order encountered ... ")
    else:
        print("Reverse order encountered ... ")

    # Get transcript infos.
    print("Read in transcript infos from --gtf to get spliced lengths ... ")
    tid2tio_dic = benchlib.gtf_read_in_transcript_infos(args.in_gtf, 
                                                        tr_ids_dic=tr_ids_dic,
                                                        correct_min_ex_order=correct_min_ex_order,
                                                        chr_style=args.chr_id_style,
                                                        remove_version_numbers=args.ignore_version_numbers,
                                                        empty_check=False)

    OUT = open(args.out_tr_list, "w")
    OUT.write("transcript_id\ttranscript_biotype\ttranscript_length\texon_count\trepresentative_transcript\tgene_name\tgene_id\tgene_biotype\n")

    for gid in gid2gio_dic:
        gene_name = gid2gio_dic[gid].gene_name
        gene_biotype = gid2gio_dic[gid].gene_biotype
        tr_ids = gid2gio_dic[gid].tr_ids
        tr_biotypes = gid2gio_dic[gid].tr_biotypes
        for idx, tid in enumerate(tr_ids):
            tr_biotype = tr_biotypes[idx]
            tr_length = tid2tio_dic[tid].tr_length
            c_exons = tid2tio_dic[tid].exon_c
            is_mpt = "No"
            if tid in mpt_ids_dic:
                is_mpt = "Yes"

            OUT.write("%s\t%s\t%i\t%i\t%s\t%s\t%s\t%s\n" % (tid, tr_biotype, tr_length, c_exons, is_mpt, gene_name, gid, gene_biotype))
    OUT.close()

    print("Output transcript IDs written to:\n%s" % (args.out_tr_list))
    print("")

