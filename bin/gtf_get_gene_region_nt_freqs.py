#!/usr/bin/env python3

import argparse
import os
from rbpbench import benchlib
import statistics


###############################################################################

def setup_argument_parser():
    """Setup argparse parser."""
    help_description = """
    Get nucleotide frequencies from all gene regions extracted from 
    --gtf and --genome.

    """
    # Define argument parser.
    p = argparse.ArgumentParser(add_help=False,
                                prog="gtf_get_gene_region_nt_freqs.py",
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

    gene_reg_bed = args.out_folder + "/gene_regions.bed"
    gene_reg_fa = args.out_folder + "/gene_regions.fa"

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

    print("Output gene regions to BED ... ")

    OUTBED = open(gene_reg_bed, "w")
    # Output transcript regions.
    for gid in gid2gio_dic:
        chr_id = gid2gio_dic[gid].chr_id
        gene_s = gid2gio_dic[gid].gene_s - 1
        gene_e = gid2gio_dic[gid].gene_e
        gene_pol = gid2gio_dic[gid].gene_pol
        OUTBED.write("%s\t%i\t%i\t%s\t0\t%s\n" %(chr_id, gene_s, gene_e, gid, gene_pol))
    OUTBED.close()

    print("Extract gene region sequences from --genome ... ")
    benchlib.bed_extract_sequences_from_fasta(gene_reg_bed, 
                                              args.in_genome, gene_reg_fa,
                                              add_param="-name",
                                              print_warnings=True)

    print("Read in extracted sequences ... ")
    gene_seqs_dic = benchlib.read_fasta_into_dic(gene_reg_fa,
                                                 dna=True,
                                                 all_uc=True,
                                                 id_check=True,
                                                 name_bed=True,
                                                 empty_check=False,
                                                 skip_n_seqs=False)

    # Give some length statistics on the extracted sequences.
    gen_lens = [len(seq) for seq in gene_seqs_dic.values()]
    print("Gene sequence length statistics:")
    print("Min:   ", min(gen_lens))
    print("Max:   ", max(gen_lens))
    print("Mean:  ", statistics.mean(gen_lens))
    print("Median:", statistics.median(gen_lens))

    print("Determine nucleotide frequencies of gene sequences ... ")

    nt_counts = {"A": 0, "C": 0, "G": 0, "T": 0}

    for seq_id, seq in gene_seqs_dic.items():
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
    print("Nucleotide frequencies of gene sequences:")
    for nt, freq in nt_freqs.items():
        print(f"{nt}: {freq:.4f}")
    print("Done.")


################################################################################







