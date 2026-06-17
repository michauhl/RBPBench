#!/usr/bin/env python3

import argparse
import os
from rbpbench import benchlib


###############################################################################

def setup_argument_parser():
    """Setup argparse parser."""
    help_description = """
    Extract transcript data (regions, sequences) from GTF file. By default 
    extracts most prominent transcript (MPT) regions for each gene. Change 
    behavior with --tr-select-mode, --tr-list, --mrna-only, --gene-list ...

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
                   nargs='+',
                   default = False,
                   help = "Supply list of transcript IDs, either via --tr-list ENST1 ENST2 or via --tr-list txt_file (one transcript ID per row) to define which transcripts to extract from --gtf (overrides MPT selection)")
    p.add_argument("--gene-list",
                   dest="gene_list",
                   type=str,
                   metavar='str',
                   nargs='+',
                   default = False,
                   help = "Supply file with gene IDs, either via --gene-list ENSG1 ENSG2 or via --gene-list txt_file (one ID per row) to define which genes to extract from --gtf for subsequent MPT selection")
    p.add_argument("--tr-select-mode",
                   dest="tr_select_mode",
                   type=int,
                   default=1,
                   choices=[1, 2, 3],
                   help="Define how to select a representative transcript from each gene 1: MPT (most prominent transcript) selection 2: Longest transcript selection. 3: Keep ALL transcripts (no selection!). Note that if --gene-list is given, only transcripts from genes in the list are considered for selection. If --tr-list is given, only transcripts from the list are used (overrides --tr-select-mode) (default: 1)")
    p.add_argument("--mrna-only",
                   dest="only_mrna",
                   default = False,
                   action = "store_true",
                   help = "Set if only mRNAs should be extracted from --gtf file. Removes all non-mRNA transcripts (default: False)")
    p.add_argument("--tr-types",
                   dest="tr_types_list",
                   type=str,
                   metavar='str',
                   nargs='+',
                   help="List of transcript biotypes to extract transcripts for. Useful for keeping only specific transcript types (e.g. --tr-types protein_coding lncRNA) from --gtf file. Note that filtering is applied on selected transcript IDs (after optional --tr-list or --gene-list specifications)")
    p.add_argument("--tr-ids-only",
                   dest="tr_ids_only",
                   default = False,
                   action = "store_true",
                   help = "Only store transcript IDs in FASTA header. By default, also add gene IDs and gene names (default: False)")
    p.add_argument("--add-ei-numbers",
                   dest="add_ei_numbers",
                   default = False,
                   action = "store_true",
                   help = "Add exon and intron numbers to BED file (default: False)")
    p.add_argument("--ignore-version-numbers",
                   dest="ignore_version_numbers",
                   default = False,
                   action = "store_true",
                   help = "Set to ignore ID version numbers in --gtf file, i.e., read in gene and transcript IDs without version numbers. This has to be set if input IDs have no version number but GTF file has (default: False)")
    p.add_argument("--skip-tec",
                   dest="skip_tec",
                   default = False,
                   action = "store_true",
                   help = "Skip genes with TEC (To be Experimentally Confirmed) gene biotype (default: False)")
    p.add_argument("--ignore-miss-tr",
                   dest="ignore_miss_tr",
                   default = False,
                   action = "store_true",
                   help = "Skip transcript IDs provided via --tr-list that are not in --gtf. By default throws an error (default: False)")
    # Promoter region extraction settings.
    p.add_argument("--prom-min-tr-len",
                   dest="prom_min_tr_len",
                   type=int,
                   metavar='int',
                   default=False,
                   help="Minimum transcript length for promoter region extraction. By default consider all selectedtranscript regions")
    p.add_argument("--prom-mrna-only",
                   dest="prom_mrna_only",
                   default = False,
                   action = "store_true",
                   help="Consider only mRNA transcript regions for promoter region extraction")
    p.add_argument("--prom-ext",
                   dest="prom_ext_up_down",
                   type=str,
                   metavar='str',
                   default="1000,100",
                   help="Up- and downstream extension of transcript start site (TSS) to define putative promoter regions, e.g. --prom-ext 500,0 for 500 upstream and 0 downstream extension (default: 1000,100)")
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
    utr5_seqs_fa = os.path.join(args.out_folder, "transcript_seqs.utr5.fa")  # 5'UTR sequences.
    utr3_seqs_fa = os.path.join(args.out_folder, "transcript_seqs.utr3.fa")  # 3'UTR sequences.
    cds_seqs_fa = os.path.join(args.out_folder, "transcript_seqs.cds.fa")  # CDS sequences.
    out_tr_list = os.path.join(args.out_folder, "transcript_infos.tsv")  # Transcript infos table.
    promoter_regions_bed = os.path.join(args.out_folder, "promoter_regions.bed")

    # If transcript list is given.
    tr_ids_dic = {}
    if args.tr_list:
        print("Using transcript IDs from --tr-list ... ")

        if os.path.isfile(args.tr_list[0]):
            tr_ids_dic = benchlib.read_ids_into_dic(args.tr_list[0], check_dic=False)
            assert tr_ids_dic, "no IDs read in from provided --tr-list file. Please provide a valid IDs file (one ID per row)"
        else:
            for tr_id in args.tr_list:
                tr_ids_dic[tr_id] = 1
            assert tr_ids_dic, "no IDs read in from provided --tr-list. Please provide transcript IDs either via --tr-list ENST1 ENST2 ... or via --tr-list txt_file (one transcript ID per row)"

        print("# of transcript IDs read in from --tr-list:", len(tr_ids_dic))

    keep_tr_types_dic = {}
    if args.tr_types_list:
        print("Using transcript biotypes from --tr-types ... ")
        for tr_type in args.tr_types_list:
            keep_tr_types_dic[tr_type] = 1
        assert keep_tr_types_dic, "no transcript biotypes read in from provided --tr-types. Please provide a valid list of transcript biotypes (e.g. protein_coding, lncRNA, etc.) to keep from --gtf file"
        print("# of transcript biotypes to keep from --gtf:", len(keep_tr_types_dic))

    gene_ids_dic = {}
    if args.gene_list:
        print("Using gene IDs list from --gene-list ... ")

        if os.path.isfile(args.gene_list[0]):
            gene_ids_dic = benchlib.read_ids_into_dic(args.gene_list,
                                                    check_dic=False)
            assert gene_ids_dic, "no IDs read in from provided --gene-list file. Please provide a valid IDs file (one ID per row)"
        else:
            for gene_id in args.gene_list:
                gene_ids_dic[gene_id] = 1
            assert gene_ids_dic, "no IDs read in from provided --gene-list. Please provide gene IDs either via --gene-list ENSG1 ENSG2 ... or via --gene-list txt_file (one gene ID per row)"

        print(f"# of gene IDs read in from --gene-list: {len(gene_ids_dic)}")

    print("Read in gene features from --gtf ... ")

    skip_gene_biotype_dic = {}
    if args.skip_tec:
        skip_gene_biotype_dic = {"TEC" : 1}

    tr2gid_dic = {}
    tr_types_dic = {}  # Store transcript biotypes in GTF file.
    gid2gio_dic = benchlib.gtf_read_in_gene_infos(args.in_gtf,
                                                  tr2gid_dic=tr2gid_dic,
                                                  tr_types_dic=tr_types_dic,
                                                  chr_style=args.chr_id_style-1,
                                                  gene_ids_dic=gene_ids_dic,
                                                  skip_gene_biotype_dic=skip_gene_biotype_dic,
                                                  remove_version_numbers=args.ignore_version_numbers,
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

    # If --tr-list, only keep transcript IDs from list (overrides --tr-select-mode setting).
    if tr_ids_dic:
        c_found = 0
        for tr_id in tr_ids_dic:
            if tr_id not in tr2gid_dic:
                if args.ignore_miss_tr:
                    print("Warning: transcript ID \"%s\" from provided --tr-list file does not appear in --gtf file (or if --gene-ids-list supplied not in resulting subset). Skipping this transcript ID ..." %(tr_id))
                else:
                    assert tr_id in tr2gid_dic, "transcript ID \"%s\" from provided --tr-list file does not appear in --gtf file (or if --gene-ids-list supplied not in resulting subset). Please provide compatible settings" %(tr_id)
            else:
                tr_ids_dic[tr_id] = tr2gid_dic[tr_id]
                c_found += 1
        assert c_found, "none of the transcript IDs from provided --tr-list file appear in --gtf file (or if --gene-ids-list supplied not in resulting subset). Please provide compatible inputs!"

        survival_pct = c_found / len(tr_ids_dic) * 100
        print(f"% of surviving transcript IDs from --tr-list: {survival_pct:.2f}% ({c_found} out of {len(tr_ids_dic)})")

    else:
        """
        Transcript selection from GTF file

        --tr-select-mode defines how to select a representative transcript for each gene from the GTF file. 
        1: Select most prominent transcript (MPT) for each gene based on various criteria such as the presence of specific tags in the GTF file.
        2: Select longest transcript for each gene.
        3: Keep all transcripts (no selection, i.e. all transcripts from GTF file are used for subsequent steps). 
        
        Note that if --gene-list is given, only transcripts from genes in the list are considered for selection. 
        
        If --tr-list is given, this step is skipped and only transcripts from the list are used (overrides --tr-select-mode).

        """

        if args.tr_select_mode == 1:

            print("Select most prominent transcript (MPT) for each gene ... ")

            tr_ids_dic = benchlib.select_mpts_from_gene_infos(gid2gio_dic,
                                    basic_tag=False,  # do not be strict (only_tsl=False too).
                                    ensembl_canonical_tag=False,
                                    prior_basic_tag=True,  # Prioritize basic tag transcript.
                                    prior_mane_select=True,  # mane select if set trumps all.
                                    prior_lncrna_primary_tag=True,  # for lncRNA genes prioritize gencode primary tagged transcripts (mane select still better but should not occur together for lncRNAs).
                                    only_tsl=False)

            assert tr_ids_dic, "most prominent transcript selection from gene infos failed. Please contact developers"
            print("# of transcript IDs (most prominent transcripts): ", len(tr_ids_dic))

        elif args.tr_select_mode == 2:

            print("Select longest transcript for each gene ... ")

            tr_ids_dic = benchlib.select_longest_tr_from_gene_infos(gid2gio_dic)

            assert tr_ids_dic, "longest transcript selection from gene infos failed. Please contact developers"
            print("# of transcript IDs (longest transcripts): ", len(tr_ids_dic))


        else:

            print("Keep all transcripts (no selection) ... ")

            tr_ids_dic = benchlib.select_all_trs_from_gene_infos(gid2gio_dic)

            assert tr_ids_dic, "selecting all transcripts from gene infos failed. Please contact developers"
            print("# of transcript IDs (all transcripts): ", len(tr_ids_dic))


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
                                                        chr_style=args.chr_id_style-1,
                                                        remove_version_numbers=args.ignore_version_numbers,
                                                        empty_check=False)

    assert tid2tio_dic, "no transcript infos read in from --gtf. Please provide a valid/compatible GTF file (e.g. from Ensembl or ENCODE) and if --tr-list is provided make sure given IDs are in GTF file"

    # (in)sanity checks.
    if not args.ignore_miss_tr:
        for tr_id in tr_ids_dic:
            assert tr_id in tid2tio_dic, "transcript ID %s not in tid2tio_dic" %(tr_id)
    else:
        # Make sure tr_ids_dic only contains transcript IDs present in tid2tio_dic.
        tr_ids_dic = {}
        for tr_id in tid2tio_dic:
            tr_ids_dic[tr_id] = 1

    for tr_id in tid2tio_dic:
        assert tr_id in tr_ids_dic, "transcript ID %s not in tr_ids_dic" %(tr_id)

    c_tr_infos = len(tid2tio_dic)
    print("# transcript features read in from --gtf:", c_tr_infos)

    """
    Optional filtering by transcript biotypes.
    I.e. update tr_ids_dic and tid2tio_dic to only contain transcripts with biotypes in keep_tr_types_dic (if provided).

    At this tr_ids_dic and tid2tio_dic should have same keys (tr_ids).
    
    """
    if keep_tr_types_dic:
        print("Filtering transcripts by biotypes from --tr-types ... ")
        tr_ids_dic, tid2tio_dic = benchlib.filter_transcripts_by_biotype(tr_ids_dic, tid2tio_dic, keep_tr_types_dic)
        c_tr_infos = len(tid2tio_dic)
        assert c_tr_infos, "no transcript features left after filtering by --tr-types. Please provide compatible settings (e.g. valid/compatible GTF file and/or --tr-types list)"
        print("# transcript features after filtering by --tr-types:", c_tr_infos)

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

        assert c_mrna_tids, "no mRNA transcripts (containing CDS) found in --gtf. Please provide a valid/compatible GTF file (e.g. from Ensembl or ENCODE). Alternatively, if --tr-list was given, make sure that the list contains mRNA transcripts, or do not set --only-mrna. Also if --mrna-only and --tr-types are specified, make sure to include mRNA transcript biotypes (e.g. protein_coding) in the list!"

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
        
        c_mrna_tids = len(tid2regl_dic)
        print("# mRNA transcripts:", c_mrna_tids)


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
                              seq_ids_only=args.tr_ids_only,  # only store transcript IDs in FASTA header.
                              to_upper=True,  # convert sequences to upper case.
                              split_size=60,  # split sequences into lines of 60 characters.
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

    # Output promoter regions to BED file.

    prom_ext_parts = args.prom_ext_up_down.split(",")
    c_prom_ext_parts = len(prom_ext_parts)
    assert c_prom_ext_parts == 2, "invalid --prom-ext argument provided (correct format: --prom-ext 1000,100, i.e., please provide two integers separated by a comma)"

    prom_ext_up = int(prom_ext_parts[0])
    prom_ext_down = int(prom_ext_parts[1])

    assert benchlib.boundary_check(prom_ext_up, 1, 100000), "set promoter upstream extension expected to be >= 1 and <= 100000"
    assert benchlib.boundary_check(prom_ext_down, 0, 100000), "set promoter downstream extension expected to be >= 0 and <= 100000"

    print("Output putative promoter regions to BED ... ")
    benchlib.output_promoter_regions_to_bed(tid2tio_dic, promoter_regions_bed,
                                    prom_min_tr_len=args.prom_min_tr_len,
                                    prom_mrna_only=args.prom_mrna_only,
                                    tr_ids_dic=tr_ids_dic,
                                    mrna_biotype_label="protein_coding",
                                    prom_ext_up=prom_ext_up,
                                    prom_ext_down=prom_ext_down)

    # Output mRNA regions.
    print("Output 5'UTR, CDS and 3'UTR sequences to FASTA ... ")

    benchlib.fasta_output_mrna_regions(tid2regl_dic, "utr5", tr_seqs_dic, utr5_seqs_fa,
                                       split=True,
                                       split_size=60,
                                       to_upper=True,
                                       id_sep=";",
                                       reg_sep=",",
                                       tr2gid_dic=tr2gid_dic,
                                       tr2gn_dic=tr2gn_dic)

    benchlib.fasta_output_mrna_regions(tid2regl_dic, "cds", tr_seqs_dic, cds_seqs_fa,
                                       split=True,
                                       split_size=60,
                                       to_upper=True,
                                       id_sep=";",
                                       reg_sep=",",
                                       tr2gid_dic=tr2gid_dic,
                                       tr2gn_dic=tr2gn_dic)

    benchlib.fasta_output_mrna_regions(tid2regl_dic, "utr3", tr_seqs_dic, utr3_seqs_fa,
                                       split=True,
                                       split_size=60,
                                       to_upper=True,
                                       id_sep=";",
                                       reg_sep=",",
                                       tr2gid_dic=tr2gid_dic,
                                       tr2gn_dic=tr2gn_dic)


    """
    Output exon and intron regions to BED file.
    
    """

    print("Output exon and intron regions to BED ... ")

    OUTBED = open(out_exon_intron_bed, "w")

    intron_label = "intron"
    exon_label = "exon" 

    for tr_id in tr_ids_dic:
        tio = tid2tio_dic[tr_id]
        chr_id = tio.chr_id
        gene_id = tio.gene_id
        gene_name = tr2gn_dic[tr_id]

        # Loop over intron regions.
        for idx, intron in enumerate(tio.intron_coords):
            intron_s = intron[0] - 1
            intron_e = intron[1]

            intron_id = intron_label
            if args.add_ei_numbers:
                intron_num = idx + 1
                intron_id = f"intron{intron_num}"

            OUTBED.write("%s\t%i\t%i\t%s;%s;%s,%s\t0\t%s\n" %(chr_id, intron_s, intron_e, tr_id, gene_id, gene_name, intron_id, tio.tr_pol))
        # Loop over exon regions.
        for idx, exon in enumerate(tio.exon_coords):
            exon_s = exon[0] - 1
            exon_e = exon[1]
            exon_len = exon_e - exon_s

            exon_id = exon_label
            if args.add_ei_numbers:
                exon_num = idx + 1
                exon_id = f"exon{exon_num}"

            OUTBED.write("%s\t%i\t%i\t%s;%s;%s,%s\t0\t%s\n" %(chr_id, exon_s, exon_e, tr_id, gene_id, gene_name, exon_id, tio.tr_pol))

    OUTBED.close()


    """
    Output transcript regions to BED file.

    """

    print("Output transcript regions to BED (BED12) ... ")

    c_out = 0
    OUTBED = open(out_tr_bed, "w")

    for tid in tr_ids_dic:
        tio = tid2tio_dic[tid]
        chr_id = tio.chr_id
        tr_s = tio.tr_s  # 1-based genomic start.
        tr_e = tio.tr_e
        tr_pol = tio.tr_pol
        gene_id = tio.gene_id
        gene_name = tr2gn_dic[tid]
        exon_coords = tio.exon_coords

        if tid not in tr_ids_dic:
            continue

        reg_id = tid + ";" + gene_id + ";" + gene_name

        c_out += 1

        # BED12 format.
        bed12_line = benchlib.transcript_to_bed12(chr_id, tr_s, tr_e, tr_pol, exon_coords, 
                                                  cds_start=tio.cds_s, 
                                                  cds_end=tio.cds_e, 
                                                  transcript_id=reg_id)

        OUTBED.write(bed12_line + "\n")
        # OUTBED.write("%s\t%i\t%i\t%s\t0\t%s\n" %(chr_id, tr_s, tr_e, reg_id, tr_pol))

    OUTBED.close()

    print("Output transcript infos to TSV file ... ")

    OUT = open(out_tr_list, "w")
    OUT.write("transcript_id\ttranscript_biotype\ttranscript_length\texon_count\tgene_name\tgene_id\tgene_biotype\tbasic_tag\tensembl_canonical\tmane_select\tprimary_tag\ttsl_id\n")

    for tid in tr_ids_dic:

        gene_name = tr2gn_dic[tid]
        gene_id = tr2gid_dic[tid]
        gene_biotype = gid2gio_dic[gene_id].gene_biotype
        tr_length = tid2tio_dic[tid].tr_length
        tr_biotype = tid2tio_dic[tid].tr_biotype
        exon_c = tid2tio_dic[tid].exon_c
        basic_tag = tid2tio_dic[tid].basic_tag
        ensembl_canonical = tid2tio_dic[tid].ensembl_canonical
        mane_select = tid2tio_dic[tid].mane_select
        primary_tag = tid2tio_dic[tid].primary_tag
        tsl_id = tid2tio_dic[tid].tsl_id

        OUT.write("%s\t%s\t%i\t%i\t%s\t%s\t%s\t%i\t%i\t%i\t%i\t%s\n" % (tid, tr_biotype, tr_length, exon_c, gene_name, gene_id, gene_biotype, basic_tag, ensembl_canonical, mane_select, primary_tag, tsl_id))

    OUT.close()

    print("Output FASTA file with transcript sequences:\n%s" %(tr_seqs_fa))
    print("Output FASTA file with 5'UTR sequences:\n%s" %(utr5_seqs_fa))
    print("Output FASTA file with CDS sequences:\n%s" %(cds_seqs_fa))
    print("Output FASTA file with 3'UTR sequences:\n%s" %(utr3_seqs_fa))
    print("Output BED file with exon and intron regions:\n%s" %(out_exon_intron_bed))
    print("Output BED file with promoter regions:\n%s" %(promoter_regions_bed))
    print("# transcript regions output to BED file:", c_out)
    print("Output BED file with transcript regions:\n%s" %(out_tr_bed))
    print("Output BED file with mRNA regions:\n%s" %(mrna_regions_bed))
    print("Output transcript infos table:\n%s" %(out_tr_list))

    print("Done.")


################################################################################

