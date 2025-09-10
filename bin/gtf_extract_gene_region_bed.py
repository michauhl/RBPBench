#!/usr/bin/env python3

import argparse
import os
import re
import gzip
from rbpbench import benchlib


###############################################################################

def setup_argument_parser():
    """Setup argparse parser."""
    help_description = """
    Extract gene regions from GTF file and create BED file with gene regions.

    """
    # Define argument parser.
    p = argparse.ArgumentParser(add_help=False,
                                prog="gtf_extract_gene_region_bed.py",
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
                   help = "Input GTF file with genomic annotations to extract gene regions from. Note that only genes on standard chromosomes (1,2,..,X,Y,MT) are currently used")
    p.add_argument("--out",
                   dest="out_bed",
                   type=str,
                   metavar='str',
                   required=True,
                   help="Output BED file to store gene regions in")
    p.add_argument("--gene-list",
                   dest="gene_list",
                   type=str,
                   metavar='str',
                   default = False,
                   help = "Supply file with gene IDs (one ID per row) to define which genes to extract from --gtf")
    p.add_argument("--bed-col4-infos",
                   dest="bed_col4_infos",
                   type=int,
                   default=1,
                   choices=[1, 2, 3],
                   help="Define what to store in BED column 4. 1: store gene_id. 2: store gene_id;gene_name. 3: store gene_id;gene_name;gene_biotype (default: 1)")
    p.add_argument("--chr-id-style",
                   dest="chr_id_style",
                   type=int,
                   default=1,
                   choices=[1, 2, 3],
                   help="Define to which chromosome ID style to convert chromosome IDs to. 1: do not change chromosome IDs. 2: convert to chr1,chr2,...,chrM style. 3: convert to 1,2,...,MT style (default: 1)")
    return p


################################################################################

def gtf_output_gene_regions_to_bed(in_gtf, out_bed,
                                   bed_col6_infos=1,
                                   gene_ids_dic=False,
                                   chr_id_style=1):
    """
    Read in gene infos into GeneInfo objects, including information on 
    transcript isoforms for the gene. Note that only features on standard 
    chromosomes (1,2,...,X Y MT) are currently used.

    Assuming gtf file with order: gene,transcript(s),exon(s) ...

    chr_style:
        0: do not change
        1: change to chr1, chr2 ...
        2: change to 1, 2, 3, ...

    """

    OUTBED = open(out_bed, "w")
    c_gene_regions = 0

    if re.search(r".+\.gz$", in_gtf):
        f = gzip.open(in_gtf, 'rt')
    else: 
        f = open(in_gtf, "r")
    for line in f:
        # Skip header.
        if re.search("^#", line):
            continue

        cols = line.strip().split("\t")
        chr_id = cols[0]
        feature = cols[2]
        feat_s = int(cols[3])  # 1-based index (see e.g. start_codon feature for proof).
        feat_e = int(cols[4])
        feat_pol = cols[6]
        infos = cols[8]

        chr_id = benchlib.check_convert_chr_id(chr_id, id_style=chr_id_style)
        # If not one of standard chromosomes, continue.
        if not chr_id:
            continue

        assert feat_e >= feat_s, "feature end < feature start in GTF file \"%s\", line \"%s\". Since both coordinates are expected to have 1-based index, this should not happen" %(in_gtf, line)

        if feature == "gene":

            m = re.search(r'gene_id "(.+?)"', infos)
            assert m, "gene_id entry missing in GTF file \"%s\", line \"%s\"" %(in_gtf, line)
            gene_id = m.group(1)

            if gene_ids_dic:
                if gene_id not in gene_ids_dic:
                    continue

            m = re.search(r'gene_name "(.+?)"', infos)
            gene_name = "-"  # optional.
            if m:
                gene_name = m.group(1)
            gene_biotype = "-"  # # optional.
            m = re.search(r'gene_biotype "(.+?)"', infos)
            if not m:
                m = re.search('gene_type "(.+?)"', infos)
            if m:
                gene_biotype = m.group(1)

            if bed_col6_infos == 1:
                bed_col6_str = gene_id
            elif bed_col6_infos == 2:
                bed_col6_str = gene_id + ";" + gene_name
            elif bed_col6_infos == 3:
                bed_col6_str = gene_id + ";" + gene_name + ";" + gene_biotype

            OUTBED.write("%s\t%i\t%i\t%s\t0\t%s\n" %(chr_id, feat_s-1, feat_e, bed_col6_str, feat_pol))

            c_gene_regions += 1

    f.close()
    OUTBED.close()

    assert c_gene_regions, "no gene regions extracted from GTF file \"%s\"" %(in_gtf)

    print("# gene regions extracted: %i" %(c_gene_regions))


################################################################################

if __name__ == '__main__':

    parser = setup_argument_parser()
    args = parser.parse_args()

    assert os.path.exists(args.in_gtf), "--gtf file \"%s\" not found" % (args.in_gtf)

    gene_ids_dic = False
    if args.gene_list:
        print("Using gene IDs list from --gene-list ... ")
        gene_ids_dic = benchlib.read_ids_into_dic(args.gene_list,
                                                  check_dic=False)
        assert gene_ids_dic, "no IDs read in from provided --gene-list file. Please provide a valid IDs file (one ID per row)"
        print("# of gene IDs (read in from --gene-list): ", len(gene_ids_dic))

    print("Read in gene features from --gtf ... ")


    print("Read in GTF and output gene regions to BED file ... ")
    gtf_output_gene_regions_to_bed(args.in_gtf, args.out_bed,
                                   bed_col6_infos=args.bed_col4_infos,
                                   gene_ids_dic=gene_ids_dic,
                                   chr_id_style=args.chr_id_style)
    print("Done.")


################################################################################