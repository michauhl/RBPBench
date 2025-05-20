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
    Extract transcript data (regions, sequences) from GTF file. By default 
    extracts most prominent transcript (MPT) regions for each gene. Change 
    behavior with --tr-list, --mrna-only or --gene-list.

    """
    # Define argument parser.
    p = argparse.ArgumentParser(add_help=False,
                                prog="gtf_extract_transcript_data.py",
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
                   help = "Input GTF file with genomic annotations to extract transcript data from. Note that only genes on standard chromosomes (1,2,..,X,Y,MT) are currently used")
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
                   help="Output folder to store transcript data")
    p.add_argument("--tr-list",
                   dest="tr_list",
                   type=str,
                   metavar='str',
                   default = False,
                   help = "Supply file with transcript IDs (one ID per row) to define which transcripts to extract from --gtf (overrides MPT selection)")
    p.add_argument("--gene-list",
                   dest="gene_list",
                   type=str,
                   metavar='str',
                   default = False,
                   help = "Supply file with gene IDs (one ID per row) to define which genes to extract from --gtf for subsequent MPT selection")
    p.add_argument("--mrna-only",
                   dest="only_mrna",
                   default = False,
                   action = "store_true",
                   help = "Set if only mRNAs should be extracted from --gtf file. Removes all non-mRNA transcripts (default: False)")
    # p.add_argument("--bed-col4-infos",
    #                dest="bed_col4_infos",
    #                type=int,
    #                default=1,
    #                choices=[1, 2],
    #                help="Define what to store in BED column 4. 1: store transcript_id. 2: store transcript_id;gene_id (default: 1)")
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
    assert os.path.exists(args.in_genome), "--genome file \"%s\" not found" % (args.in_genome)

    if not os.path.exists(args.out_folder):
        os.makedirs(args.out_folder)

    mrna_regions_bed = os.path.join(args.out_folder, "mrna_regions.bed")  # UTR CDS regions on mRNAs (transcript context).
    out_tr_bed = os.path.join(args.out_folder, "transcript_regions.bed")  # transcript regions on genome.
    out_exon_intron_bed = os.path.join(args.out_folder, "exon_intron_regions.bed")  # Exon regions of transcripts on genome.
    tr_seqs_fa = os.path.join(args.out_folder, "transcript_seqs.fa")  # Transcript sequences (spliced) in FASTA format.
    tr_seqs_len_out = os.path.join(args.out_folder, "transcript_seqs_len.txt")  # Transcript lengths.

    gene_ids_dic = False
    if args.gene_list:
        print("Using gene IDs list from --gene-ids-list ... ")
        gene_ids_dic = benchlib.read_ids_into_dic(args.gene_list,
                                                  check_dic=False)
        assert gene_ids_dic, "no IDs read in from provided --gene-list file. Please provide a valid IDs file (one ID per row)"
        print("# of gene IDs (read in from --gene-list): ", len(gene_ids_dic))

    print("Read in gene features from --gtf ... ")
    tr2gid_dic = {}
    tr_types_dic = {}  # Store transcript biotypes in GTF file.
    gid2gio_dic = benchlib.gtf_read_in_gene_infos(args.in_gtf,
                                                tr2gid_dic=tr2gid_dic,
                                                tr_types_dic=tr_types_dic,
                                                chr_style=args.chr_id_style,
                                                gene_ids_dic=gene_ids_dic,
                                                empty_check=False)

    assert gid2gio_dic, "no gene infos read in from --gtf. Please provide a valid/compatible GTF file (e.g. from Ensembl or ENCODE)"
    c_gene_infos = len(gid2gio_dic)
    print("# gene features read in from --gtf:", c_gene_infos)

    # Get transcript ID -> gene name mapping.
    tr2gn_dic = {}
    for tr_id in tr2gid_dic:
        gene_id = tr2gid_dic[tr_id]
        gene_name = gid2gio_dic[gene_id].gene_name
        tr2gn_dic[tr_id] = gene_name

    # Get most prominent transcripts or if --tr-list is set, read in transcript IDs.
    tr_ids_dic = {}
    if args.tr_list:
        tr_ids_dic = benchlib.read_ids_into_dic(args.tr_list,
                                                check_dic=False)
        assert tr_ids_dic, "no IDs read in from provided --tr-list file. Please provide a valid IDs file (one ID per row)"
        for tr_id in tr_ids_dic:
            assert tr_id in tr2gid_dic, "transcript ID \"%s\" from provided --tr-list file does not appear in --gtf file (or if --gene-ids-list supplied not in resulting subset). Please provide compatible settings" %(tr_id)
            tr_ids_dic[tr_id] = tr2gid_dic[tr_id]
        print("# of transcript IDs (read in from --tr-list): ", len(tr_ids_dic))
    else:
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
                                                        chr_style=args.chr_id_style,
                                                        empty_check=False)

    assert tid2tio_dic, "no transcript infos read in from --gtf. Please provide a valid/compatible GTF file (e.g. from Ensembl or ENCODE)"

    # (in)sanity checks.
    for tr_id in tr_ids_dic:
        assert tr_id in tid2tio_dic, "transcript ID %s not in tid2tio_dic"
    for tr_id in tid2tio_dic:
        assert tr_id in tr_ids_dic, "transcript ID %s not in tr_ids_dic"

    c_tr_infos = len(tid2tio_dic)
    print("# transcript features read in from --gtf:", c_tr_infos)

    # If --only-mrna, only select mRNAs, which also triggers mRNA region occupancy plots generation.
    tid2regl_dic = {}
    if args.only_mrna:
        # Get mRNA region lengths (from transcript IDs with CDS feature).
        print("Get mRNA region lengths ... ")
        # tid2tio_dic contains infos for transcripts in tr_ids_dic (i.e. either MPT or --tr-list defined ones).
        tid2regl_dic = benchlib.get_mrna_region_lengths(tid2tio_dic)
        assert tid2regl_dic, "tid2regl_dic empty. If --tr-list was set, this means none of supplied transcript IDs contain a CDS. If --tr-list was not set, this means that none of the MPTs contain a CDS. In this case please provide a valid/compatible GTF file (e.g. from Ensembl or ENCODE) or contact developers"
        c_mrna_tids = len(tid2regl_dic)
        if args.tr_list:
            print("# mRNA transcripts (containing CDS) from --tr-list:", c_mrna_tids)
        else:
            print("# mRNA transcripts (containing CDS, out of MPT selected set):", c_mrna_tids)

        assert c_mrna_tids, "no mRNA transcripts (containing CDS) found in --gtf. Please provide a valid/compatible GTF file (e.g. from Ensembl or ENCODE). Alternatively, if --tr-list was given, make sure that the list contains mRNA transcripts, or do not set --only-mrna"

        # Output mRNA regions (5'UTR CDS 3'UTR) to BED.
        print("Output mRNA regions to BED ... ")
        benchlib.output_mrna_regions_to_bed(tid2regl_dic, mrna_regions_bed)
    
        # Assign mRNA transcripts as transcripts to extract sequences for.
        tr_ids_dic = tid2regl_dic
    else:
        # Also extract mRNA regions to BED file, but do not restrict transcript IDs to mRNAs.
        print("Get mRNA region lengths ... ")
        tid2regl_dic = benchlib.get_mrna_region_lengths(tid2tio_dic)
        print("Output mRNA regions to BED ... ")
        benchlib.output_mrna_regions_to_bed(tid2regl_dic, mrna_regions_bed,
                                            tr2gid_dic=tr2gid_dic,
                                            tr2gn_dic=tr2gn_dic,
                                            empty_check=False)
    

    """
    Output transcript sequences to FASTA file.

    """

    # Get transcript sequences.
    print("Extract transcript sequences ... ")
    tr_seqs_dic = benchlib.get_transcript_sequences_from_gtf(tid2tio_dic, args.in_genome,
                                                             tr_ids_dic=tr_ids_dic,
                                                             tmp_out_folder=args.out_folder)

    # Output sequences to FASTA.
    print("Output transcript sequences to FASTA ... ")
    benchlib.fasta_output_dic(tr_seqs_dic, tr_seqs_fa,
                              tr2gid_dic=tr2gid_dic,  # add gene ID to header.
                              tr2gn_dic=tr2gn_dic,  # add gene name to header.
                              split=True)

    # Transcript sequence lengths.
    tr_seq_len_dic = {}
    for tr_id in tr_seqs_dic:
        tr_seq_len_dic[tr_id] = len(tr_seqs_dic[tr_id])

    # Output transcript lengths to tr_seqs_len_out.
    print("Output transcript lengths to file ... ")
    benchlib.output_tr_lengths(tr_seq_len_dic, tr_seqs_len_out,
                               tr2gid_dic=tr2gid_dic,  # add gene ID.
                               tr2gn_dic=tr2gn_dic)  # add gene name.


    """
    Output exon and intron regions to BED file.
    
    """

    print("Output exon and intron regions to BED ... ")

    OUTBED = open(out_exon_intron_bed, "w")

    for tr_id in tr_ids_dic:
        tio = tid2tio_dic[tr_id]
        chr_id = tio.chr_id
        gene_id = tio.gene_id
        gene_name = tr2gn_dic[tr_id]

        # Loop over intron regions.
        for intron in tio.intron_coords:
            intron_s = intron[0] - 1
            intron_e = intron[1]
            OUTBED.write("%s\t%i\t%i\t%s;%s;%s;intron\t0\t%s\n" %(chr_id, intron_s, intron_e, tr_id, gene_id, gene_name, tio.tr_pol))
        # Loop over exon regions.
        for exon in tio.exon_coords:
            exon_s = exon[0] - 1
            exon_e = exon[1]
            exon_len = exon_e - exon_s
            OUTBED.write("%s\t%i\t%i\t%s;%s;%s;exon\t0\t%s\n" %(chr_id, exon_s, exon_e, tr_id, gene_id, gene_name, tio.tr_pol))

    OUTBED.close()


    """
    Output transcript regions to BED file.

    """

    print("Output transcript regions to BED ... ")

    c_out = 0
    OUTBED = open(out_tr_bed, "w")

    for tid in tid2tio_dic:
        tio = tid2tio_dic[tid]
        chr_id = tio.chr_id
        tr_s = tio.tr_s - 1
        tr_e = tio.tr_e
        tr_pol = tio.tr_pol
        gene_id = tio.gene_id
        gene_name = tr2gn_dic[tid]

        if tid not in tr_ids_dic:
            continue

        reg_id = tid + ";" + gene_id + ";" + gene_name

        c_out += 1

        OUTBED.write("%s\t%i\t%i\t%s\t0\t%s\n" %(chr_id, tr_s, tr_e, reg_id, tr_pol))

    OUTBED.close()

    print("Output FASTA file with transcript sequences:\n%s" %(tr_seqs_fa))

    print("Output BED file with exon and intron regions:\n%s" %(out_exon_intron_bed))

    print("# transcript regions output to BED file:", c_out)

    print("Output BED file with transcript regions:\n%s" %(out_tr_bed))

    print("Output BED file with mRNA regions:\n%s" %(mrna_regions_bed))

    print("Done.")


################################################################################

