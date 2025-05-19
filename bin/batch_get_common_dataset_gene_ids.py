#!/usr/bin/env python3

from rbpbench import benchlib
import pandas as pd

import argparse
import os

"""
Example --in file format (gene_region_occupancies.tsv):
gene_id	gene_name	transcript_id	IGF2BP1,clipper_idr,hepg2_eclip	IGF2BP1,clipper_idr,k562_eclip
ENSG00000173153	ESRRA	ENST00000000442	1	0
ENSG00000004478	FKBP4	ENST00000001008	0	1
ENSG00000003249	DBNDD1	ENST00000002501	1	1
...

"""

###############################################################################

def setup_argument_parser():
    """Setup argparse parser."""
    help_description = """
    Extract gene IDs which occur in all datasets, given the gene_region_occupancies.tsv
    table from RBPBench batch output folder.
    
    """
    p = argparse.ArgumentParser(add_help=False,
                                prog="batch_get_common_dataset_gene_ids.py",
                                description=help_description,
                                formatter_class=argparse.MetavarTypeHelpFormatter)
    p.add_argument("-h", "--help",
                   action="help",
                   help="Print help message")
    p.add_argument("--in",
                   dest="in_file",
                   type=str,
                   metavar='str',
                   required = True,
                   help = "Gene region occupancies .tsv file from RBPBench batch output folder")
    p.add_argument("--out",
                   dest="out_file",
                   type=str,
                   metavar='str',
                   required=True,
                   help="Output text file to store common gene IDs (optionally further filtered by --gene-list) in")
    p.add_argument("--gene-list",
                   dest="gene_ids_list",
                   type=str,
                   metavar='str',
                   default = False,
                   help = "Supply file with gene IDs (one ID per row) to further filter obtained --in gene IDs, i.e., keep only gene IDs from --in that occur in --gene-list")
    return p


###############################################################################

if __name__ == '__main__':

    parser = setup_argument_parser()
    args = parser.parse_args()

    assert os.path.exists(args.in_file), "--in file \"%s\" not found" % (args.in_file)

    # Optionally read in gene IDs to filter the output by.
    gene_ids_dic = {}
    if args.gene_ids_list:
        assert os.path.exists(args.gene_ids_list), "--gene-list file \"%s\" not found" % (args.gene_ids_list)
        gene_ids_dic = benchlib.read_ids_into_dic(args.gene_ids_list,
                                                  check_dic=False)
        print("# of gene IDs (read in from --gene-list):", len(gene_ids_dic))

    # Read --in table.
    df = pd.read_csv(args.in_file, sep='\t')

    # Number of gene IDs in the input file.
    num_gene_ids = df.shape[0]
    print("# of gene IDs in --in file):", num_gene_ids)

    # Get rows that have all 1s after the 'transcript_id' column.
    data_columns = df.columns.tolist()[3:]

    mask = (df[data_columns] == 1).all(axis=1)

    gene_ids_with_all_ones = df.loc[mask, 'gene_id'].tolist()
    print("# of gene IDs (with all 1s):", len(gene_ids_with_all_ones))

    if gene_ids_dic:
        # Filter the gene IDs by the ones in the gene_ids_dic.
        gene_ids_with_all_ones = [gene_id for gene_id in gene_ids_with_all_ones if gene_id in gene_ids_dic]
        print("# of gene IDs (after filtering by --gene-list):", len(gene_ids_with_all_ones))

    # Output the gene IDs to the output file.
    with open(args.out_file, 'w') as f:
        for gene_id in gene_ids_with_all_ones:
            f.write(f"{gene_id}\n")
    print("Output written to:", args.out_file)
    print("Done.")
