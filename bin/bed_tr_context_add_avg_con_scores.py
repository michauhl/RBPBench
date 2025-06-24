#!/usr/bin/env python3

import argparse
import pyBigWig
import numpy as np
from rbpbench import benchlib
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
    Given a BED file with transcript context regions (regions on transcripts / 
    with transcript coordinates, not genomic coordinates), map these regions 
    back to the genome to get an average conservation score for each of these 
    regions (given phyloP or phastCons bigWig files, plus a GTF file to map 
    transcript to genomic regions). 

    """
    # Define argument parser.
    p = argparse.ArgumentParser(add_help=False,
                                prog="bed_tr_context_add_avg_con_scores.py",
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
                   help = "Input transcript regions BED file to extract genomic conservation scores for")
    p.add_argument("--con",
                   dest="con_sc_file",
                   type=str,
                   metavar='str',
                   required = True,
                   help = "Genomic .bigWig file (phastCons or phyloP) with conservation scores (has to be compatible with --bed + --gtf)")
    p.add_argument("--gtf",
                   dest="in_gtf",
                   type=str,
                   metavar='str',
                   required = True,
                   help = "Input GTF file with genomic annotations to extract transcript data from. Note that only genes on standard chromosomes (1,2,..,X,Y,MT) are currently used")
    p.add_argument("--out",
                   dest="out_bed",
                   type=str,
                   metavar='str',
                   required = True,
                   help = "Output BED file with added transcript region conservation scores column")
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

    assert os.path.exists(args.in_bed), "given --bed BED file \"%s\" not found" % (args.in_bed)
    assert os.path.exists(args.con_sc_file), "given --con bigWig file \"%s\" not found" % (args.con_sc_file)
    assert os.path.exists(args.in_gtf), "given --gtf GTF file \"%s\" not found" % (args.in_gtf)

    print("Read in transcript IDs (col1) from --bed file ...") 

    tr_ids_dic = benchlib.bed_read_chr_ids_dic(args.in_bed)
    
    assert tr_ids_dic, "no transcript IDs read in from provided --bed file (column 1). Please provide a valid BED file with transcript regions"

    print("# of transcript IDs from --bed file: ", len(tr_ids_dic))

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

    assert tid2tio_dic, "no transcript infos read in from --gtf. Please provide compatible --bed and --gtf files"

    print("# of transcript infos read in: ", len(tid2tio_dic))

    for tid in tr_ids_dic:
        assert tid in tid2tio_dic, "transcript ID \"%s\" from --bed file not found in --gtf file! Please provide compatible --bed and --gtf files" %(tid)

    """
    Extract genomic conservation scores for transcript regions.
    
    """

    print("Extract genomic conservation scores for transcript regions ... ")

    # Open conservation scores.
    con_sc_data = pyBigWig.open(args.con_sc_file)
    
    BEDOUT = open(args.out_bed, "w")
    
    with open(args.in_bed) as f:
        for line in f:
            cols = line.strip().split("\t")
            assert len(cols) >= 6, "invalid --in BED format. Please provide valid BED file (i.e., >= 6 column format)"
            tr_id = cols[0]
            tr_s = int(cols[1]) + 1  # make 1-based.
            tr_e = int(cols[2])
            reg_id = cols[3]
            reg_sc = cols[4]
            reg_pol = cols[5]
            reg_len = tr_e - tr_s + 1
            
            # Transcript features.
            exon_coords = tid2tio_dic[tr_id].exon_coords
            chr_id = tid2tio_dic[tr_id].chr_id
            gen_s = tid2tio_dic[tr_id].tr_s  # 1-based genomic transcript start.
            gen_e = tid2tio_dic[tr_id].tr_e  # 1-based genomic transcript end.
            gen_pol = tid2tio_dic[tr_id].tr_pol

            # gen_coords format: [[1901, 1903], [1501, 1600], [1198, 1200]]
            gen_coords = benchlib.map_transcript_to_genomic(gen_s, gen_e, gen_pol, exon_coords, tr_s, tr_e)

            assert gen_coords, "no genomic coordinates found for transcript ID \"%s\" in --bed file. Please provide compatible --bed and --gtf files" %(tr_id)

            all_chunk_scores = []

            for gen_chunk in gen_coords:
                gen_chunk_s = gen_chunk[0] - 1  # make 0-based for conservation score extraction.
                gen_chunk_e = gen_chunk[1]
                                
                try:
                    # Get conservation scores for the region.
                    scores = con_sc_data.values(chr_id, gen_chunk_s, gen_chunk_e, numpy=False)
                    c_scores = len(scores)
                    assert c_scores == reg_len, "# of extracted scores != region length (%i != %i) for line:\n%s" %(c_scores, reg_len, line)
                    # Convert NaN values to 0.0.
                    scores = [0.0 if np.isnan(s) else s for s in scores]
                    all_chunk_scores.extend(scores)

                except RuntimeError:
                    print(f"Skipping chunk region {chr_id}:{gen_chunk_s}-{gen_chunk_e} (coordinates not in bigWig)")

            # Average score over all positions.
            avg_score = statistics.mean(all_chunk_scores) if all_chunk_scores else 0.0

            cols.append(str(avg_score))

            BEDOUT.write("\t".join(cols) + "\n")

    f.closed
    BEDOUT.close()







