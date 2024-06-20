#!/usr/bin/env python3

import argparse
import os
import re
import subprocess


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
                   help="Define to which chromosome ID style to convert chromosomes to. 1: do not change chromosome IDs. 2: convert to chr1,chr2,...,chrM style. 3: convert to 1,2,...,MT style (default: 1)")
    return p


################################################################################

"""
Define additional chromosome names to be supported.

From Drosophila melanogaster:
chr2L
chr2R
chr3L
chr3R
2L
2R
3L
3R

"""

add_chr_names_dic = {
    "chr2L" : 1,
    "chr2R" : 1,
    "chr3L" : 1,
    "chr3R" : 1,
    "2L" : 1,
    "2R" : 1,
    "3L" : 1,
    "3R" : 1
}

################################################################################

def check_convert_chr_id(chr_id,
                         chr_id_style=1):
    """
    Check and convert chromosome IDs to format:
    chr1, chr2, chrX, ...
    If chromosome IDs like 1,2,X, .. given, convert to chr1, chr2, chrX ..
    Return False if given chr_id not standard and not convertable.

    Filter out scaffold IDs like:
    GL000009.2, KI270442.1, chr14_GL000009v2_random
    chrUn_KI270442v1 ...

    chr_id_style:
        Defines to which style chromosome IDs should be converted to.
        1: Do not convert at all, just return chr_id.
        2: to chr1,chr2,...,chrM style.
        3: to 1,2,...,MT style.

    """
    assert chr_id, "given chr_id empty"

    if chr_id_style == 1: # If id_style == 1.
        return chr_id

    elif chr_id_style == 2:
        if re.search("^chr", chr_id):
            if chr_id in add_chr_names_dic or re.search("^chr[\dMXY]+$", chr_id):
                return chr_id
            else:
                return False
        else:
            # Convert to "chr" IDs.
            if chr_id == "MT": # special case MT -> chrM.
                return "chrM"
            if chr_id in add_chr_names_dic or re.search("^[\dXY]+$", chr_id):
                return "chr" + chr_id
            else:
                return False

    elif chr_id_style == 3:

        if re.search("^chr", chr_id):
            if chr_id == "chrM": # special case chrM -> MT.
                return "MT"
            if chr_id in add_chr_names_dic or re.search("^chr[\dXY]+$", chr_id):
                # Cut out chr suffix.
                m = re.search("chr(.+)", chr_id)
                assert m, "no match for regex search"
                chr_suffix = m.group(1)
                return chr_suffix
            else:
                return False

        else:
            if chr_id == "MT": # special case MT.
                return chr_id
            if chr_id in add_chr_names_dic or re.search("^[\dXY]+$", chr_id):
                return chr_id
            else:
                return False
    else:
        assert False, "invalid chr_id_style set"


################################################################################

def gtf_output_gene_regions_to_bed(in_gtf, out_bed,
                                   bed_col6_infos=1,
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

    if re.search(".+\.gz$", in_gtf):
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

        chr_id = check_convert_chr_id(chr_id, chr_id_style=chr_id_style)
        # If not one of standard chromosomes, continue.
        if not chr_id:
            continue

        assert feat_e >= feat_s, "feature end < feature start in GTF file \"%s\", line \"%s\". Since both coordinates are expected to have 1-based index, this should not happen" %(in_gtf, line)

        if feature == "gene":

            m = re.search('gene_id "(.+?)"', infos)
            assert m, "gene_id entry missing in GTF file \"%s\", line \"%s\"" %(in_gtf, line)
            gene_id = m.group(1)
            m = re.search('gene_name "(.+?)"', infos)
            gene_name = "-"  # optional.
            if m:
                gene_name = m.group(1)
            gene_biotype = "-"  # # optional.
            m = re.search('gene_biotype "(.+?)"', infos)
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

    print("Read in GTF and output gene regions to BED file ... ")
    gtf_output_gene_regions_to_bed(args.in_gtf, args.out_bed,
                                   bed_col6_infos=args.bed_col4_infos,
                                   chr_id_style=args.chr_id_style)
    print("Done.")


################################################################################