from rbpbench import benchlib
import sys

"""
Check GTF file for uniqueness of MANE select and primary tags in lncRNA genes.

"""

if len(sys.argv) != 2:
    print("Usage: python gtf_check_mane_primary_tags.py gtf_file_path", file=sys.stderr)
    sys.exit(1)

in_gtf = sys.argv[1]

gid2gio_dic = benchlib.gtf_read_in_gene_infos(in_gtf,
                                              empty_check=False)

for gid in gid2gio_dic:
    gene_biotype = gid2gio_dic[gid].gene_biotype
    c_mane_select_tags = 0
    for tag in gid2gio_dic[gid].tr_mane_select_tags:
        c_mane_select_tags += tag
    
    if c_mane_select_tags > 1:
        print("Gene %s has %d MANE select tags!" % (gid, c_mane_select_tags))

    # if not c_mane_select_tags:
    #     print("Gene %s has no MANE select tags!" % (gid))

    c_primary_tags = 0
    for tag in gid2gio_dic[gid].tr_primary_tags:
        c_primary_tags += tag
    if c_primary_tags > 1 and gene_biotype == "lncRNA":
        print("lncRNA gene %s has %d primary tags!" % (gid, c_primary_tags))

    if gene_biotype == "lncRNA" and c_primary_tags == 0:
        print("lncRNA gene %s has no primary tags!" % (gid))