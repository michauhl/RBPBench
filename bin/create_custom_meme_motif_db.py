#!/usr/bin/env python3

from rbpbench import benchlib
import argparse
import pandas as pd
import os


###############################################################################

def setup_argument_parser():
    """Setup argparse parser."""
    help_description = """
    Create MEME motif database from input table with regex infos, 
    which can be used as custom motif database in RBPBench (--custom-db option).

    """
    # Define argument parser.
    p = argparse.ArgumentParser(add_help=False,
                                prog="create_custom_meme_motif_db.py",
                                description=help_description,
                                formatter_class=argparse.MetavarTypeHelpFormatter)

    # Required arguments.
    p.add_argument("-h", "--help",
                   action="help",
                   help="Print help message")
    p.add_argument("--in",
                   dest="in_table",
                   type=str,
                   metavar='str',
                   required = True,
                   help = "Input tab-separated table file containing regex and motif infos to create MEME motif database (mandatory columns: rbp_id, motif_id, regex)")
    p.add_argument("--out",
                   dest="out_folder",
                   type=str,
                   metavar='str',
                   required=True,
                   help="Output folder to store MEME motif database genrated from --in regex infos in. Contains MEME XML file and table file used as input files (--custom-db-meme-xml xml_file --custom-db-info table_file) or folder (--custom-db db_folder) to RBPBench)")
    return p


################################################################################

if __name__ == '__main__':

    parser = setup_argument_parser()
    args = parser.parse_args()

    assert os.path.exists(args.in_table), "--in table file \"%s\" not found" % (args.in_table)

    if not os.path.exists(args.out_folder):
        os.makedirs(args.out_folder)

    out_xml = args.out_folder + "/seq_motifs.meme"
    out_info_tsv = args.out_folder + "/info.txt"

    mandatory_columns = ["rbp_id", "motif_id", "regex"]
    optional_columns = ["organism", "gene_id", "function_ids", "reference", "experiments", "comments"]

    # # Get database paths + function ID -> function descriptions mapping.
    # benchlib_path = os.path.dirname(benchlib.__file__)
    # db_path = benchlib_path + "/content"

    # fid2desc_file = db_path + "/rbp_function_list.32728246.tsv"
    # fid2desc_dic, desc2fid_dic = benchlib.get_fid2desc_mapping(fid2desc_file)

    df = pd.read_csv(args.in_table, sep="\t")

    # Check if mandatory columns are present.
    for col in mandatory_columns:
        assert col in df.columns, f"column \"{col}\" not found in --in table {args.in_table}"

    mid2rid_dic = {}
    query_motif_blocks_dic = {}
    OUTINFO = open(out_info_tsv, "w")
    OUTINFO.write("RBP_motif_ID\tRBP_name\tMotif_type\tOrganism\tGene_ID\tFunction_IDs\tReference\tExperiment\tComments\n")

    for index, row in df.iterrows():
        rbp_id = row["rbp_id"]
        motif_id = row["motif_id"]
        if motif_id in mid2rid_dic:
            assert False, "non-unique motif ID \"%s\" found in --in table %s. Please provide unique motif IDs" % (motif_id, args.in_table)
        mid2rid_dic[motif_id] = rbp_id

        print("Read in RBP ID %s, motif ID %s ..." % (rbp_id, motif_id))

        motif_type = "meme_xml"
        motif_stats = benchlib.MotifInfos(rbp_id, motif_id, motif_type)

        function_ids = "-"
        gene_id = "-"
        organism = "-"
        reference = "-"
        experiments = "-"
        comments = "-"

        if "organism" in df.columns:
            if pd.isna(row["organism"]):
                motif_stats.organism = "-"
            else:
                motif_stats.organism = str(row["organism"])
            organism = motif_stats.organism
        if "gene_id" in df.columns:
            if pd.isna(row["gene_id"]):
                motif_stats.gene_id = "-"
            else:
                motif_stats.gene_id = str(row["gene_id"])
            gene_id = motif_stats.gene_id
        if "function_ids" in df.columns:
            if pd.isna(row["function_ids"]):
                motif_stats.function_ids = "-"
                function_ids = "-"
            else:
                function_ids_list = str(row["function_ids"]).split(";")
                motif_stats.function_ids = function_ids_list
                function_ids = ";".join(function_ids_list)
        if "reference" in df.columns:
            if pd.isna(row["reference"]):
                motif_stats.reference = "-"
            else:
                motif_stats.reference = str(row["reference"])
            reference = motif_stats.reference
        if "experiments" in df.columns:
            if pd.isna(row["experiments"]):
                motif_stats.experiments = "-"
            else:
                motif_stats.experiments = str(row["experiments"])
            experiments = motif_stats.experiments
        if "comments" in df.columns:
            if pd.isna(row["comments"]):
                motif_stats.comments = "-"
            else:
                motif_stats.comments = str(row["comments"])
            comments = motif_stats.comments

        function_ids_str = ";".join(motif_stats.function_ids)
        if not function_ids_str:
            function_ids_str = "-"

        OUTINFO.write(f"{motif_id}\t{rbp_id}\t{motif_type}\t{organism}\t{gene_id}\t{function_ids}\t{reference}\t{experiments}\t{comments}\n")

        regex = row["regex"]

        assert benchlib.is_valid_regex(regex), "given --in table regex \"%s\" is not a valid regular expression. Please provide valid expression (only square brackets allowed as special characters, e.g. AC[AC]GTA)" % (regex)

        special_chars = r"[,\d.^$*+?{}()|\\]"

        regex = benchlib.remove_special_chars_from_str(regex,
                                                    reg_ex=special_chars)

        assert regex, "empty string after removing special chars \"%s\" from --regex. Please provide a valid regex with DNA letters and optionally square brackets" %(special_chars)

        print("Regex after removing special chars:", regex)

        # Get all sequence parts from regex.
        seq_parts_list = benchlib.get_seq_parts_from_regex(regex)

        for seq in seq_parts_list:
            assert benchlib.seq_check_alphabet(seq, alphabet=["A", "C", "G", "T"]), "sequence \"%s\" derived from --in string has non-DNA letters in it. Please provide DNA sequences / regular expressions" %(seq)

            seq_motif_block = benchlib.seq_parts_to_motif_block(seq_parts_list)

            query_motif_blocks_dic[motif_id] = seq_motif_block

    OUTINFO.close()

    out_str, c_added_motifs = benchlib.blocks_to_xml_string(query_motif_blocks_dic, query_motif_blocks_dic,
                                                            mid2rid_dic=mid2rid_dic)
    benchlib.output_string_to_file(out_str, out_xml)
    print("# of added motifs to XML from --in:", c_added_motifs)

    assert os.path.exists(out_xml), "no motifs file written to output folder. Please contact developers"

    print("")
    print("Motif database written to folder:\n%s" % (args.out_folder))
    print("Motif database MEME XML file:\n%s" % (out_xml))
    print("Motif database info table file:\n%s" % (out_info_tsv))
    print("")

################################################################################