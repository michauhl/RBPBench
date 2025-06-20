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
                                                  empty_check=False)

    assert gid2gio_dic, "no gene infos read in from --gtf. Please provide a valid/compatible GTF file (e.g. from Ensembl or ENCODE)"

    for gid in gene_ids_dic:
        if gid not in gid2gio_dic:
            print("WARNING: gene ID \"%s\" not found in GTF file!" % (gid))

    OUT = open(args.out_tr_list, "w")
    OUT.write("transcript_id\ttranscript_biotype\tgene_name\tgene_id\tgene_biotype\n")

    for gid in gid2gio_dic:
        gene_name = gid2gio_dic[gid].gene_name
        gene_biotype = gid2gio_dic[gid].gene_biotype
        tr_ids = gid2gio_dic[gid].tr_ids
        tr_biotypes = gid2gio_dic[gid].tr_biotypes
        for idx, tid in enumerate(tr_ids):
            tr_biotype = tr_biotypes[idx]
            OUT.write("%s\t%s\t%s\t%s\t%s\n" % (tid, tr_biotype, gene_name, gid, gene_biotype))
    OUT.close()

    print("Output transcript IDs written to:\n%s" % (args.out_tr_list))
    print("")