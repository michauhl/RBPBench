#!/usr/bin/env python3

import argparse
import os
from rbpbench import benchlib
import statistics


###############################################################################

def setup_argument_parser():
    """Setup argparse parser."""
    help_description = """
    Get nucleotide frequencies from most prominent transcript (MPT) sequences 
    (introns excluded) extracted from --gtf and --genome.

    """
    # Define argument parser.
    p = argparse.ArgumentParser(add_help=False,
                                prog="gtf_get_mpt_nt_freqs.py",
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
                   help = "Input GTF file with genomic annotations to extract transcript regions from. Note that only genes on standard chromosomes (1,2,..,X,Y,MT) are currently used")
    p.add_argument("--genome",
                   dest="in_genome",
                   type=str,
                   metavar='str',
                   required = True,
                   help = "Genomic sequences file (currently supported formats: FASTA)")
    p.add_argument("--out",
                   dest="out_folder",
                   type=str,
                   metavar='str',
                   required=True,
                   help="Output folder to store sequences in")
    return p



################################################################################

if __name__ == '__main__':

    parser = setup_argument_parser()
    args = parser.parse_args()

    assert os.path.exists(args.in_gtf), "--gtf file \"%s\" not found" % (args.in_gtf)
    assert os.path.exists(args.in_genome), "--genome file \"%s\" not found" % (args.in_genome)

    if not os.path.exists(args.out_folder):
        os.makedirs(args.out_folder)

    tr_seqs_fa = args.out_folder + "/transcript_sequences.fa"

    """
    Get chromosome IDs from --genome.
    """
    print("Get --genome FASTA headers ... ")
    chr_ids_dic = benchlib.get_fasta_headers(args.in_genome)

    """
    Guess chromosome ID style.

    chr_style:
        1: chr1, chr2, ..., chrX, chrM
        2: 1, 2, ... , X, MT

    """
    chr_style = 0  # no changes to chromosome IDs in GTF files.
    
    if "1" in chr_ids_dic:
        assert "chr1" not in chr_ids_dic, "inconsistent chromosome IDs in --genome FASTA file (both chr1 and 1)"
    if "chr1" in chr_ids_dic:
        assert "1" not in chr_ids_dic, "inconsistent chromosome IDs in --genome FASTA file (both chr1 and 1)"
    if "chr1" in chr_ids_dic:
        chr_style = 1
    if "1" in chr_ids_dic:
        chr_style = 2

    print("Read in gene features from --gtf ... ")
    tr2gid_dic = {}
    gid2gio_dic = benchlib.gtf_read_in_gene_infos(args.in_gtf,
                                                  tr2gid_dic=tr2gid_dic,
                                                  check_chr_ids_dic=chr_ids_dic,
                                                  chr_style=chr_style,
                                                  empty_check=False)

    assert gid2gio_dic, "no gene infos read in from --gtf. Please provide a valid/compatible GTF file (e.g. from Ensembl or ENCODE)"
    c_gene_infos = len(gid2gio_dic)
    print("# gene features read in from --gtf:", c_gene_infos)

    # Get most prominent transcripts from gene infos.
    print("Select most prominent transcript (MPT) for each gene ... ")
    tr_ids_dic = benchlib.select_mpts_from_gene_infos(gid2gio_dic,
                            basic_tag=False,  # do not be strict (only_tsl=False too).
                            ensembl_canonical_tag=False,
                            prior_basic_tag=True,  # Prioritize basic tag transcript.
                            only_tsl=False)
    assert tr_ids_dic, "most prominent transcript selection from gene infos failed. Please contact developers"
    print("# of transcript IDs (most prominent transcripts): ", len(tr_ids_dic))

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
                                                        chr_style=chr_style,
                                                        empty_check=False)

    # Get transcript sequences.
    print("Extract transcript sequences ... ")
    tr_seqs_dic = benchlib.get_transcript_sequences_from_gtf(tid2tio_dic, args.in_genome,
                                                             tr_ids_dic=tr_ids_dic,
                                                             dna=True,
                                                             all_uc=True,
                                                             tmp_out_folder=args.out_folder)

    # Give some length statistics on the extracted sequences.
    tr_lens = [len(seq) for seq in tr_seqs_dic.values()]
    print("Transcript sequence length statistics:")
    print("Min:   ", min(tr_lens))
    print("Max:   ", max(tr_lens))
    print("Mean:  ", statistics.mean(tr_lens))
    print("Median:", statistics.median(tr_lens))

    # Output sequences to FASTA.
    print("Output transcript sequences to FASTA ... ")
    benchlib.fasta_output_dic(tr_seqs_dic, tr_seqs_fa,
                              split=True)

    print("Determine nucleotide frequencies of MPT sequences ... ")

    nt_counts = {"A": 0, "C": 0, "G": 0, "T": 0}

    for seq_id, seq in tr_seqs_dic.items():
        for nt in seq:
            if nt in nt_counts:
                nt_counts[nt] += 1

    # Total number of nucleotides.
    total_nt = sum(nt_counts.values())

    # Calculate the frequency of each nucleotide.
    if total_nt > 0:
        nt_freqs = {nuc: count / total_nt for nuc, count in nt_counts.items()}
    else:
        nt_freqs = {nuc: 0 for nuc in nt_counts}

    # Print the nucleotide frequencies
    print("Nucleotide frequencies of MPT sequences:")
    for nt, freq in nt_freqs.items():
        print(f"{nt}: {freq:.4f}")
    print("Done.")


################################################################################







