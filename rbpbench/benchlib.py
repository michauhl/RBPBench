from distutils.spawn import find_executable
from typing import Optional
import os
import re
import subprocess
import gzip
import matplotlib.pyplot as plt
from matplotlib_venn import venn2, venn3
from venn import venn
from logomaker import Logo
from markdown import markdown
import pandas as pd
import plotly.express as px
from math import log10
import textdistance
from upsetplot import plot


"""

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
~~~~~~ OPEN FOR BUSINESS ~~~~~~
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


~~~~~~~~~~~~~
Run doctests
~~~~~~~~~~~~~

python3 -m doctest benchlib.py



"""


################################################################################

def calc_edit_dist_query_list(query_str, lst):
    """
    Calculate edit distance between query string and list of strings.

    >>> query_str = 'AAAC'
    >>> lst = ['AAAC', 'AAAA', 'AACA']
    >>> calc_edit_dist_query_list(query_str, lst)
    {'AAAA': 1, 'AACA': 2}

    """

    assert lst, "given lst empty"
    assert query_str, "given query_str empty"
    comp_pairs = []
    for str in lst:
        if str != query_str:
            comp_pairs.append([query_str, str])
    pair_dist_dic = {}
    for pair in comp_pairs:
        s1 = pair[0]
        s2 = pair[1]
        d = textdistance.levenshtein(s1, s2)
        pair_dist_dic[s2] = d

    return pair_dist_dic


################################################################################

def calc_edit_dist(s1, s2):
    """
    Calculate edit distance (Levenshtein distance) between two strings s1, s2.

    >>> s1 = 'HALLO'
    >>> s2 = 'HELLO'
    >>> calc_edit_dist(s1, s1)
    0
    >>> calc_edit_dist(s1, s2)
    1
    
    """
    assert s1, "given s1 empty"
    assert s2, "given s2 empty"
    return textdistance.levenshtein(s1, s2)


################################################################################

def output_con_table_results(out_tsv, pval_ll, rbp_list):
    """
    Output contingency table test p-values to file.
    
    """
    TABOUT = open(out_tsv, "w")

    heads = ["#rbp_id"]
    for rbp_id in rbp_list:
        heads.append(rbp_id)
    heads_str = "\t".join(heads)
    TABOUT.write("%s\n" %(heads_str))
    for i,r in enumerate(pval_ll):
        # print("i:", i)
        # print("rbp_list[i]:", rbp_list[i])
        # sys.exit()
        row = [rbp_list[i]]
        for e in pval_ll[i]:
            row.append(str(e))
        row_str = "\t".join(row)
        TABOUT.write("%s\n" %(row_str))
    TABOUT.close()


################################################################################

def dir_get_files(file_dir,
                  file_ending=False,
                  check=True):
    """
    Return list of files from given file_dir.
    E.g. file_ending="bed" to filter for .bed files.

    >>> test_dir = "test_data"
    >>> dir_get_files(test_dir, file_ending="fasta")
    ['test.fasta']

    """

    from os import listdir
    from os.path import isfile, join
    dir_files = [f for f in listdir(file_dir) if isfile(join(file_dir, f))]
    if check:
        assert dir_files, "given directory \"%s\" contains no files" %(file_dir)
    # If filter for file ending true.
    if file_ending:
        new_files = []
        for df in dir_files:
            if re.search(".+\.%s" %(file_ending), df):
                new_files.append(df)
        if check:
            assert new_files, "no files left after filtering by file ending \"%s\"" %(file_ending)
        return sorted(new_files)
    else:
        return sorted(dir_files)


################################################################################

def make_contingency_table_2x2(region_labels_dic, idx1, idx2):
    """
    Make a contingency table 2x2, using region_labels_dic.
    region_labels_dic format:
    region_id -> [False, True, False ... ] 
    with list number of RBP IDs (len_rbp_list), alphabetically sorted.
    True: region covered by RBP at idx (i.e. motif hits)
    False: region not covered by RBP at idx (i.e. no motif hits)
    
    Return table format:
    table = [[A, B], [C, D]]

    https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.fisher_exact.html
                   List 1              Not in List 1
    List 2         A                   B
    Not in List 2  C                   D

    table = [[A, B], [C, D]]

    >>> region_labels_dic = {'r1': [True, False], 'r2': [False, True], 'r3': [True, True], 'r4': [False, False], 'r5': [True, True]}
    >>> make_contingency_table_2x2(region_labels_dic, 0, 1)
    [[2, 1], [1, 1]]

    """

    cont_a = 0
    cont_b = 0
    cont_c = 0
    cont_d = 0

    for reg_id in region_labels_dic:
        in1 = region_labels_dic[reg_id][idx1]
        in2 = region_labels_dic[reg_id][idx2]
        if in1 and in2:
            cont_a += 1
        elif not in1 and in2:
            cont_b += 1
        elif in1 and not in2:
            cont_c += 1
        else: # not in1 and not in2
            cont_d += 1

    table = [[cont_a, cont_b], [cont_c, cont_d]]
    return table


################################################################################

def read_in_cm_blocks(cm_file):
    """
    Read in .cm file, containing covariance model(s). Return single model 
    text blocks in dictionary, with accession IDs as keys.

    """

    idx = -1
    blocks_list = []
    acc_list = []

    with open(cm_file) as f:
        for line in f:
            if re.search("^INFERNAL", line):
                blocks_list.append(line)
                idx += 1
            elif re.search("^ACC\s+\w+", line):
                m = re.search("ACC\s+(\w+)", line)
                acc_id = m.group(1)
                acc_list.append(acc_id)
                blocks_list[idx] += line
            else:
                blocks_list[idx] += line
    f.closed

    assert blocks_list, "no covariance models read in from %s (blocks_list empty)" %(cm_file)
    assert 2*len(blocks_list) == len(acc_list), "invalid cm_file %s format (2*len(blocks_list) != len(acc_list))" %(cm_file)


    blocks_dic = {}
    for i, acc_id in enumerate(acc_list[::2]):
        if not acc_id in blocks_dic:
            blocks_dic[acc_id] = blocks_list[i]
        else:
            assert False, "structure motif accession (ACC) ID %s appears twice in cm_file %s" %(acc_id, cm_file)
    return blocks_dic


################################################################################

def output_cm_blocks_to_file(cm_blocks_dic, motif_ids_dic, out_file):
    """
    Output covariance models to file.

    """

    OUTSTR = open(out_file, "w")
    for motif_id in motif_ids_dic:
        if motif_id in cm_blocks_dic:
            OUTSTR.write(cm_blocks_dic[motif_id])
    OUTSTR.close()


################################################################################

def check_cm_file(cm_file, cmstat_out,
                  check=True):
    """
    Check if .cm file is valid, and return accession IDs, using cmstat.

    """
    assert is_tool("cmstat"), "cmstat not in PATH"
    assert os.path.exists(cm_file), "cm_file %s does not exist" %(cm_file)

    check_cmd = "cmstat " + cm_file + " > " + cmstat_out
    output = subprocess.getoutput(check_cmd)

    error = False
    if output:
        error = True
    assert error == False, "cmstat is complaining (supplied --user-cm file probably has invalid format):\n%s\n%s\nPlease provide valid covariance model file" %(check_cmd, output)

    """
    # idx   name                  accession      nseq  eff_nseq   clen      W   bps  bifs  model     cm    hmm
    # ----  --------------------  ---------  --------  --------  -----  -----  ----  ----  -----  -----  -----
         1  Histone3              RF00032          46     46.00     46     57     6     0     cm  0.699  0.604
    #
    """

    # Check for models / get model accession IDs from cmstat.
    acc_ids_dic = {}
    with open(cmstat_out) as f:
        for line in f:
            if re.search("^#", line):
                continue
            cols_pre = line.strip().split(" ")
            # Remove empty column values.
            cols = []
            for c in cols_pre:
                if c:
                    cols.append(c)
            acc_id = cols[2]
            acc_ids_dic[acc_id] = 1
    f.closed

    if check:
        assert acc_ids_dic, "no valid covariance models found in supplied --user-cm file (checked by cmstat). Please provide valid covariance model file"

    return acc_ids_dic


################################################################################

def read_cm_acc(in_cm):
    """
    Read in accession IDs from covariance model file.
    Accession ID expected to be appearing in one line, with format:
    ACC      accession_id
    Return IDs dictionary, with values == appearance counts (expected to be 
    2 for each individual ID)

    >>> in_cm = "test_data/SLBP_USER.cm"
    >>> read_cm_acc(in_cm)
    {'RF00032': 2}

    """
    assert os.path.exists(in_cm), "in_cm %s does not exist" %(in_cm)
    
    acc_dic = {}

    with open(in_cm) as f:
        for line in f:
            if re.search("^ACC\s+\w+", line):
                m = re.search("ACC\s+(\w+)", line)
                acc_id = m.group(1)
                if acc_id in acc_dic:
                    acc_dic[acc_id] += 1
                else:
                    acc_dic[acc_id] = 1
    f.closed

    return acc_dic


################################################################################

def get_fasta_headers(in_fa):
    """
    Get FASTA header IDs.
    This grep version appears to be much faster than reading in file via 
    Python line by line.

    >>> in_fa = "test_data/test.fa"
    >>> get_fasta_headers(in_fa)
    {'seq1': 1, 'seq2': 1}
    
    """
    assert os.path.exists(in_fa), "in_fa %s does not exist" %(in_fa)
    check_cmd = "grep '>' " + in_fa 
    output = subprocess.getoutput(check_cmd)

    seq_ids_dic = {}
    for line in output.split('\n'):
        if re.search("^>", line):
            m = re.search(">(.+)", line)
            seq_id = m.group(1)
            seq_ids_dic[seq_id] = 1

    assert seq_ids_dic, "no FASTA header IDs read in from in_fa %s" %(in_fa)
    return seq_ids_dic


################################################################################

def run_cmsearch(in_fa, in_cm, out_tab,
                 error_check=False,
                 call_dic=None,
                 params="-g --tformat fasta --toponly --incT 1 -T 1 --default"):
    """
    Run cmsearch on a given covariance model (.cm) file and a given FASTA 
    file. Note that .cm file can contain several models.

    cmsearch parameters:
    -g global alignment

    """
    assert is_tool("cmsearch"), "cmsearch not in PATH"
    assert os.path.exists(in_fa), "in_fa %s does not exist" %(in_fa)
    assert os.path.exists(in_cm), "in_cm %s does not exist" %(in_cm)

    params += " --tblout %s " %(out_tab)

    check_cmd = "cmsearch " + params + " " + in_cm + " " + in_fa
    output = subprocess.getoutput(check_cmd)

    if call_dic is not None:
        call_dic["cmsearch_call"] = check_cmd

    if error_check:
        error = False
        if output:
            error = True
        assert error == False, "cmsearch is complaining:\n%s\n%s" %(check_cmd, output)


################################################################################

def run_fast_fimo(in_fa, in_meme_xml, out_tsv,
                  pval_thr=0.001,
                  params="--norc --verbosity 1 --skip-matched-sequence --text",
                  nt_freqs_file=False,
                  call_dic=None,
                  error_check=True):
    """
    Run FIMO with --skip-matched-sequence --text, to increase speed (disables 
    q-value and matched_sequence calculations), also no other output files.

    Example:
    fimo --norc --verbosity 1 --skip-matched-sequence --text test_pum2.meme test_pum2.fa > check.tsv
    motif_id	motif_alt_id	sequence_name	start	stop	strand	score	p-value	q-value	matched_sequence
    PUM2_2		chr19:100000-100085(+)	21	28	+	15.7959	0.0000171		
    PUM2_2		chr18:100000-100096(+)	17	24	+	15.7959	0.0000171		
    PUM2_2		chr18:100000-100096(+)	75	82	+	15.7959	0.0000171
    
    """
    assert is_tool("fimo"), "fimo not in PATH"
    assert os.path.exists(in_fa), "in_fa %s does not exist" %(in_fa)
    assert os.path.exists(in_meme_xml), "in_meme_xml %s does not exist" %(in_meme_xml)
    
    if nt_freqs_file:
        params += " --bfile %s " %(nt_freqs_file)

    check_cmd = "fimo " + params + " --thresh " + str(pval_thr) + " " + in_meme_xml + " " + in_fa + " > " + out_tsv
    output = subprocess.getoutput(check_cmd)

    if call_dic is not None:
        call_dic["fimo_call"] = check_cmd

    if error_check:
        error = False
        if output:
            error = True
        assert error == False, "fimo is complaining:\n%s\n%s" %(check_cmd, output)


################################################################################

def run_fimo(in_fa, in_meme_xml, out_folder,
             pval_thr=0.001,
             nt_freqs_file=False,
             params="--norc --verbosity 1"):
    """
    
    Example call:
    fimo -oc test2_out --norc --verbosity 1 test_pum2.meme test_pum2.fa
    fimo -oc test3_out --norc --verbosity 1 --skip-matched-sequence --text test_pum2.meme test_pum2.fa

    fimo -oc test3_out --norc --verbosity 1 test_pum2.meme test_pum2.fa

    fimo --norc --verbosity 1 --skip-matched-sequence --text test_pum2.meme test_pum2.fa > check.tsv
    motif_id	motif_alt_id	sequence_name	start	stop	strand	score	p-value	q-value	matched_sequence
    PUM2_2		chr19:100000-100085(+)	21	28	+	15.7959	0.0000171		
    PUM2_2		chr18:100000-100096(+)	17	24	+	15.7959	0.0000171		
    PUM2_2		chr18:100000-100096(+)	75	82	+	15.7959	0.0000171		

    

    Some FIMO options:    
     --no-qvalue:
        Do not compute a q-value for each p-value. The q-value calculation 
        is that of Benjamini and Hochberg (1995).
    To speed up, add:
    --skip-matched-sequence
    --text
    
    """
    assert is_tool("fimo"), "fimo not in PATH"
    assert os.path.exists(in_fa), "in_fa %s does not exist" %(in_fa)
    assert os.path.exists(in_meme_xml), "in_meme_xml %s does not exist" %(in_meme_xml)

    if nt_freqs_file:
        params += " --bfile %s " %(nt_freqs_file)

    check_cmd = "fimo " + params + " --thresh " + str(pval_thr) + " -oc " + out_folder + " " + in_meme_xml + " " + in_fa
    output = subprocess.getoutput(check_cmd)
    error = False
    if output:
        error = True
    assert error == False, "fimo is complaining:\n%s\n%s" %(check_cmd, output)


################################################################################

def bed_extract_sequences_from_fasta(in_bed, in_fa, out_fa,
                                     print_warnings=False,
                                     ignore_errors=False):
    """
    Extract sequences from genome (provide genome .fasta file).
    Gzipped FASTA not supported by bedtools getfasta so far ...

    print_warnings:
        Instead of raising an error when output is encountered,
        print output (i.e. bedtools getfasta warnings)

    # Example call.
    bedtools getfasta -fi hg38.fa -bed PUM1_K562_IDR_peaks.uniq_ids.bed -s -fo test_out.fa

    # Content example:
    >chr8:9772198-9772297(+)
    gtttcactaaggagaagaatatttgcgtggttttcttattacagacttactttcaaagggggctggggggaggagtaattatatattggagaaagggga
    >chr11:65962035-65962092(+)
    GCAGAGCCCTCCGAGCGGCGCGTGAAGCGGGAGAAGCGCGATGACGGCTACGAGGCC
    ...
    There is no line break in sequence output ...
    bedtools getfasta -fi hg38.fa.gz -bed PUM1_K562_IDR_peaks.uniq_ids.bed -s

    >>> in_bed = "test_data/test2.bed"
    >>> in_fa = "test_data/test2.fa"
    >>> out_fa = "test_data/test2.tmp.fa"
    >>> bed_extract_sequences_from_fasta(in_bed, in_fa, out_fa, ignore_errors=True)
    >>> read_fasta_into_dic(out_fa, dna=True, all_uc=True)
    {'chr1:5-10(+)': 'GGGTT', 'chr2:5-10(-)': 'AACCC'}
    
    """

    assert is_tool("bedtools"), "bedtools not in PATH"
    assert os.path.exists(in_bed), "in_bed %s does not exist" %(in_bed)
    assert os.path.exists(in_fa), "in_fa %s does not exist" %(in_fa)

    check_cmd = "bedtools getfasta -fi " + in_fa +  " -bed " + in_bed + " -s  -fo " +  out_fa
    output = subprocess.getoutput(check_cmd)
    if print_warnings:
        if output:
            print("bedtools getfasta warnings:\n", output)
    else:
        error = False
        if output:
            error = True
        if not ignore_errors:
            assert error == False, "bedtools getfasta is complaining:\n%s\n%s" %(check_cmd, output)


################################################################################

def read_fasta_into_dic(fasta_file,
                        seqs_dic=False,
                        ids_dic=False,
                        dna=False,
                        report=1,
                        all_uc=False,
                        empty_check=True,
                        id_check=True,
                        skip_data_id="set",
                        skip_n_seqs=True):
    """
    Read in FASTA sequences, store in dictionary and return dictionary.
    FASTA file can be plain text or gzipped (watch out for .gz ending).

    >>> test_fasta = "test_data/test.fa"
    >>> read_fasta_into_dic(test_fasta)
    {'seq1': 'acguACGUacgu', 'seq2': 'ugcaUGCAugcaACGUacgu'}

    """
    if not seqs_dic:
        seqs_dic = {}
    seq_id = ""

    # Open FASTA either as .gz or as text file.
    if re.search(".+\.gz$", fasta_file):
        f = gzip.open(fasta_file, 'rt')
    else:
        f = open(fasta_file, "r")
    for line in f:
        if re.search(">.+", line):
            m = re.search(">(.+)", line)
            seq_id = m.group(1)
            if id_check:
                assert seq_id not in seqs_dic, "non-unique FASTA header \"%s\" in \"%s\"" % (seq_id, fasta_file)
            if ids_dic:
                if seq_id in ids_dic:
                    seqs_dic[seq_id] = ""
            else:
                seqs_dic[seq_id] = ""
        elif re.search("[ACGTUN]+", line, re.I):
            m = re.search("([ACGTUN]+)", line, re.I)
            seq = m.group(1)
            if seq_id in seqs_dic:
                if dna:
                    # Convert to DNA, concatenate sequence.
                    seq = seq.replace("U","T").replace("u","t")
                else:
                    # Convert to RNA, concatenate sequence.
                    seq = seq.replace("T","U").replace("t","u")
                if all_uc:
                    seq = seq.upper()
                seqs_dic[seq_id] += seq
    f.close()

    # Check if sequences read in.
    if empty_check:
        assert seqs_dic, "no sequences read in (input FASTA file \"%s\" empty or mal-formatted?)" %(fasta_file)
    # If sequences with N nucleotides should be skipped.
    c_skipped_n_ids = 0
    if skip_n_seqs:
        del_ids = []
        for seq_id in seqs_dic:
            seq = seqs_dic[seq_id]
            if re.search("N", seq, re.I):
                if report == 1:
                    print ("WARNING: sequence with seq_id \"%s\" in file \"%s\" contains N nucleotides. Discarding sequence ... " % (seq_id, fasta_file))
                c_skipped_n_ids += 1
                del_ids.append(seq_id)
        for seq_id in del_ids:
            del seqs_dic[seq_id]
        assert seqs_dic, "no sequences remaining after deleting N containing sequences (input FASTA file \"%s\")" %(fasta_file)
        if c_skipped_n_ids:
            if report == 2:
                print("# of N-containing %s regions discarded:  %i" %(skip_data_id, c_skipped_n_ids))
    return seqs_dic


################################################################################

def fasta_output_dic(fasta_dic, fasta_out,
                     split=False,
                     out_ids_dic=False,
                     header_add_sc_dic=False,
                     to_upper=False,
                     split_size=60):
    """
    Output FASTA sequences dictionary (sequence_id -> sequence) to fasta_out.

    split:
        Split FASTA sequence for output to file
    split_size:
        Split size (row width)
    to_upper:
        Convert sequences to uppercase.
    out_ids_dic:
        IDs to output dictionary.
    header_add_sc_dic:
        ID to scoring mapping.
        Add a score to the header, so header format changes from "id1"
        to "id1,10"

    """
    # Check.
    assert fasta_dic, "given dictionary fasta_dic empty"
    # Write sequences to FASTA file.
    OUTFA = open(fasta_out,"w")
    for seq_id in fasta_dic:
        seq = fasta_dic[seq_id]
        if out_ids_dic:
            if seq_id not in out_ids_dic:
                continue
        if to_upper:
            seq = seq.upper()
        out_id = seq_id
        if header_add_sc_dic:
            out_id = out_id + "," + str(header_add_sc_dic[seq_id])
        if split:
            OUTFA.write(">%s\n" %(out_id))
            for i in range(0, len(seq), split_size):
                OUTFA.write("%s\n" %((seq[i:i+split_size])))
        else:
            OUTFA.write(">%s\n%s\n" %(out_id, seq))
    OUTFA.close()


################################################################################

def read_in_xml_motifs(meme_xml_file, 
                       check=True):
    """
    Read in XML motifs, store in blocks dictionary.

    """

    assert os.path.exists(meme_xml_file), "meme_xml_file %s not found" % (meme_xml_file)

    raw_text = ""
    with open(meme_xml_file, "r") as f:
        raw_text = f.read()
    assert raw_text, "nothing read in from MEME XML file %s" %(meme_xml_file)

    # Get motif blocks.
    motif_blocks_dic = extract_motif_blocks(raw_text)

    if check:
        assert motif_blocks_dic, "motif_blocks_dic empty (malformatted MEME XML file provided?)"

    return motif_blocks_dic


################################################################################

def extract_motif_blocks(raw_text):
    """
    Extract MEME XML motif blocks, store in dictionary:
    motif_id -> motif text block, as list of lines

    """
    motif_blocks_dic = {}
    motif_id = ""
    lines = raw_text.strip().split('\n')
    for l in lines:
        if re.search("MOTIF\s\w+", l):
            m = re.search("MOTIF (\w+)", l)
            motif_id = m.group(1)
        else:
            if motif_id and l:
                if motif_id in motif_blocks_dic:
                    motif_blocks_dic[motif_id].append(l)
                else:
                    motif_blocks_dic[motif_id] = [l]

    return motif_blocks_dic


################################################################################

def blocks_to_xml_string(motif_blocks_dic, motif_ids_dic,
                         out_xml=False):
    """
    Return MEME XML string based on given motif IDs dictionary and available
    motif blocks dictionary.
    """

    header = """MEME version 5

ALPHABET= ACGT

strands: +

Background letter frequencies
A 0.250000 C 0.250000 G 0.250000 T 0.250000

"""
    if out_xml:
        OUTXML = open(out_xml, "w")

    xml_str = header
    c_added_motifs = 0
    for motif_id in motif_ids_dic:
        if motif_id in motif_blocks_dic: # Only for sequence motifs.
            block = motif_blocks_dic[motif_id]
            block_str = "\n".join(block)
            block_str_final = "MOTIF " + motif_id + " \n" + block_str + "\n\n"
            xml_str += block_str_final
            c_added_motifs += 1
            if out_xml:
                OUTXML.write(block_str_final)
    if out_xml:
        OUTXML.close()
    return xml_str, c_added_motifs


################################################################################

def output_string_to_file(out_str, out_file):
    """
    Output string to file.
    """

    OUTSTR = open(out_file, "w")
    OUTSTR.write(out_str)
    OUTSTR.close()


################################################################################

def get_rbp_id_mappings(rbp2ids_file):
    """
    Read in file mapping RBP names to motif IDs and motif types.
    Return dictionaries with:
    rbp_name -> motif_ids_list
    motif_id -> motif_type

    FILE FORMAT:

    RBP_motif_id	RBP_name	Motif_type	Organism
    AGGF1_1	AGGF1	meme_xml	human
    AGGF1_2	AGGF1	meme_xml	human
    AKAP1_1	AKAP1	meme_xml	human
    BCCIP_1	BCCIP	meme_xml	human
    BUD13_1	BUD13	meme_xml	human
    ...
    RF00032	SLBP	cm	human

    """
    name2ids_dic = {}
    id2type_dic = {}
    id2org_dic = {}

    with open(rbp2ids_file) as f:
        for line in f:
            if re.search("^RBP_motif_id", line):
                continue
            cols = line.strip().split("\t")
            motif_id = cols[0]
            rbp_name = cols[1]
            motif_type = cols[2]
            organism = cols[3]
            
            if organism != "human":
                rbp_name = rbp_name + "_" + organism

            id2org_dic[motif_id] = organism
            id2type_dic[motif_id] = motif_type
            if rbp_name in name2ids_dic:
                name2ids_dic[rbp_name].append(motif_id)
            else:
                name2ids_dic[rbp_name] = [motif_id]
    f.closed

    assert name2ids_dic, "no RBP IDs read in from %s" %(rbp2ids_file)

    return name2ids_dic, id2type_dic, id2org_dic


################################################################################

def get_uniq_gen_size(gen_sites_bed):
    """
    Get unique genomic space size, which the genomic sites inside
    gen_sites_bed cover.

    >>> gen_sites_bed = "test_data/test_gen_size.bed"
    >>> get_uniq_gen_size(gen_sites_bed)
    2500

    sort -k1,1 -k2,2n in_sites.filtered.bed | mergeBed -i stdin -s -c 4 -o distinct -delim ";"
    """

    params_str = '-s -c 4 -o distinct -delim ";"'
    check_cmd = "sort -k1,1 -k2,2n " + gen_sites_bed + " | mergeBed -i stdin " + params_str
    output = subprocess.getoutput(check_cmd)

    reg_len_sum = 0
    for line in output.split('\n'):
        cols = line.strip().split("\t")
        reg_s = int(cols[1])
        reg_e = int(cols[2])
        reg_len = reg_e - reg_s
        reg_len_sum += reg_len

    assert reg_len_sum, "no merged regions obtained from \"%s\"" %(check_cmd)
    return reg_len_sum


################################################################################

def is_tool(name):
    """Check whether tool "name" is in PATH."""
    return find_executable(name) is not None


################################################################################

def bed_check_format(bed_file, gimme=False):
    """
    Check whether given BED file is not empty + has >= 6 columns.

    >>> bed_file = "test_data/test2.bed"
    >>> bed_check_format(bed_file, gimme=True)
    True
    
    """
    pols = ["+", "-"]
    okay = False
    with open(bed_file) as f:
        for line in f:
            cols = line.strip().split("\t")
            assert len(cols) >= 6, "invalid --in BED format. Please provide valid BED file (i.e., >= 6 column format)"
            reg_pol = cols[5]
            assert reg_pol in pols, "invalid polarity in --in BED column 6 (found: \"%s\", expected \"+\" or \"-\"). Please provide valid BED file (i.e., >= 6 column format), or use --unstranded option" %(reg_pol)
            okay = True
            break
    f.closed
    assert okay, "invalid --in BED format (file empty?)"
    if gimme:
        return okay


################################################################################

def check_report_in_file(in_file):
    """
    Check if in_file is RBP or motif stats file.
    
    >>> in_file = 'test_data/rbp_hit_stats.tsv'
    >>> check_report_in_file(in_file)
    'rbp_stats'
    >>> in_file = 'test_data/motif_hit_stats.tsv'
    >>> check_report_in_file(in_file)
    'motif_stats'
    >>> in_file = 'test_data/test2.bed'
    >>> check_report_in_file(in_file)
    '?'

    """
    type = "?"

    with open(in_file) as f:
        for line in f:
            cols = line.strip().split("\t")
            if len(cols) == 27: # RBP stats.
                if cols[0] == "data_id" and cols[26] == "internal_id":
                    type = "rbp_stats"
            elif len(cols) == 20:
                if cols[0] == "data_id" and cols[19] == "internal_id":
                    type = "motif_stats"
            break
    f.closed
    return type


################################################################################

def reg_get_core_id(reg_id):
    """
    From region ID reg_id, get core ID without strand.

    So from chr1:100-101:+ to chr1:100-101

    >>> reg_id = "chr1:100-101(+)"
    >>> reg_get_core_id(reg_id)
    'chr1:100-101'
    
    """

    if re.search("\w+:\d+-\d+\([+|-]\)", reg_id):
        m = re.search("(\w+):(\d+)-(\d+)\(", reg_id)
        core_id = "%s:%s-%s" %(m.group(1), m.group(2), m.group(3))
        return core_id
    else:
        assert False, "region ID has invalid format (%s). Please contact developers" %(reg_id)


################################################################################

def get_hit_id_elements(hit_id):
    """
    From  hit ID to ID elements list.

    Hit ID format:
    chr1:100-110(+),motif_id

    >>> hit_id = "chr1:100-110(+),motif_id"
    >>> get_hit_id_elements(hit_id)
    ['chr1', '100', '110', '+', 'motif_id']
    
    """

    if re.search("\w+:\d+-\d+\([+|-]\),\w+", hit_id):
        m = re.search("(\w+):(\d+)-(\d+)\(([+|-])\),(.+)", hit_id)
        id_elements = [m.group(1), m.group(2), m.group(3), m.group(4), m.group(5)]
        return id_elements
    else:
        assert False, "hit ID has invalid format (%s). Please contact developers" %(hit_id)


################################################################################

def get_center_position(start, end):
    """
    Get center position (1-based), given a (genomic) start (0-based) and
    end coordinate (1-based).

    >>> get_center_position(10, 11)
    11
    >>> get_center_position(1000,2000)
    1501
    >>> get_center_position(11, 20)
    17

    """
    # If region has length of 1, return end position.
    center_pos = end
    # Otherwise calculate new center position.
    if not end - start == 1:
        center_pos = round( ( (end - start) / 2 ) + start ) + 1
    return center_pos


################################################################################

def diff_two_files_identical(file1, file2):
    """
    Check whether two files are identical. Return true if diff reports no
    differences.

    >>> file1 = "test_data/file1"
    >>> file2 = "test_data/file2"
    >>> diff_two_files_identical(file1, file2)
    True
    >>> file1 = "test_data/test1.bed"
    >>> diff_two_files_identical(file1, file2)
    False

    """
    same = True
    check_cmd = "diff " + file1 + " " + file2
    output = subprocess.getoutput(check_cmd)
    if output:
        same = False
    return same


################################################################################

def bed_extend_bed(in_bed, out_bed, 
                   ext_lr=10,
                   cp_mode=1,
                   remove_dupl=False,
                   chr_ids_dic=None):
    """
    Center + extend BED file regions (rbpbench dist), optionally filter by 
    chromosome IDs.

    >>> in_bed = 'test_data/test.bed'
    >>> out_bed = 'test_data/test.tmp.bed'
    >>> exp_bed = 'test_data/test.exp.bed'
    >>> reg_stats_dic = bed_extend_bed(in_bed, out_bed, ext_lr=2, cp_mode=1)
    >>> diff_two_files_identical(out_bed, exp_bed)
    True
    >>> in_fa = "test_data/test2.fa"
    >>> out_fa = "test_data/test.tmp.fa"
    >>> bed_extract_sequences_from_fasta(exp_bed, in_fa, out_fa, ignore_errors=True)
    >>> read_fasta_into_dic(out_fa, dna=True, all_uc=True)
    {'chr1:3-8(+)': 'AAGGG', 'chr2:7-12(-)': 'GGAAC'}

    """

    OUTBED = open(out_bed, "w")

    c_in = 0
    c_out = 0
    c_chr_filter = 0
    c_dupl_filter = 0
    reg_len_sum = 0
    seen_reg_dic = {}

    with open(in_bed) as f:
        for line in f:
            cols = line.strip().split("\t")
            chr_id = cols[0]
            reg_s = int(cols[1])
            reg_e = int(cols[2])
            reg_id = cols[3]
            reg_sc = cols[4]
            reg_pol = cols[5]

            c_in += 1

            # Check pol.
            assert reg_pol == "+" or reg_pol == "-", "invalid region polarity (strand) given in --in BED file %s column 6 (found \"%s\", expecting \"+\" or \"-\"). For unstranded search see --unstranded option" %(in_bed, reg_pol)
            # Check chr_id
            if chr_ids_dic is not None:
                if chr_id not in chr_ids_dic:
                    c_chr_filter += 1
                    continue
            # Remove duplicates.
            if remove_dupl:
                reg_str = "%s:%i-%i(%s)" % (chr_id, reg_s, reg_e, reg_pol)
                if reg_str not in seen_reg_dic:
                    seen_reg_dic[reg_str] = 1
                else:
                    c_dupl_filter += 1
                    continue

            # Center position.
            center_s = reg_s
            center_e = center_s + 1
            if cp_mode == 1:
                if reg_pol == "-":
                    center_s = reg_e - 1
                    center_e = center_s + 1
            elif cp_mode == 2:
                center_e = get_center_position(reg_s, reg_e)
                center_s = center_e - 1
            elif cp_mode == 3:
                center_s = reg_e - 1
                center_e = center_s + 1
                if reg_pol == "-":
                    center_s = reg_s
                    center_e = center_s + 1
            else:
                assert False, "invalid cp_mode given"

            assert center_e - center_s == 1, "invalid center positioning encountered"

            # Extend.
            new_s = center_s - ext_lr
            new_e = center_e + ext_lr
            if new_s < 0:
                new_s = 0
            reg_len = new_e - new_s
            reg_len_sum += reg_len

            c_out += 1
            
            OUTBED.write("%s\t%s\t%s\t%s\t%s\t%s\n" % (chr_id, new_s, new_e, reg_id, reg_sc, reg_pol))

    f.closed
    OUTBED.close()

    stats_dic = {}
    stats_dic["c_in"] = c_in
    stats_dic["c_out"] = c_out
    stats_dic["c_chr_filter"] = c_chr_filter
    stats_dic["c_dupl_filter"] = c_dupl_filter
    stats_dic["reg_len_sum"] = reg_len_sum

    return stats_dic


################################################################################

def bed_read_rows_into_dic(in_bed, to_list=False):
    """
    Read in .bed file rows into dictionary.
    Mapping is region ID -> bed row.

    """

    id2row_dic = {}

    with open(in_bed) as f:
        for line in f:
            cols = line.strip().split("\t")
            chr_id = cols[0]
            site_s = cols[1]
            site_e = cols[2]
            site_id = cols[3]
            site_sc = cols[4]
            site_pol = cols[5]

            row = "%s\t%s\t%s\t%s\t%s\t%s" %(chr_id, site_s, site_e, site_id, site_sc, site_pol)
            if to_list:
                id2row_dic[site_id] = cols
            else:
                id2row_dic[site_id] = row

    f.closed

    return id2row_dic


################################################################################

def bed_filter_extend_bed(in_bed, out_bed,
                          ext_up=0,
                          ext_down=0,
                          remove_dupl=True,
                          reg2sc_dic=None,
                          score_col=5,
                          chr_ids_dic=None,
                          unstranded=False):
    """
    Filter / extend in_bed BED file, output to out_bed.

    remove_dupl:
        Remove duplicated regions, i.e., regions with same 
        chr_id, reg_s, reg_e, reg_pol
    unstranded:
        If True, output both strands of each region (+ and -), so given one 
        input region, two output regions are generated.
        
    """

    OUTBED = open(out_bed, "w")

    assert score_col > 0, "invalid score column (< 1)"

    c_in = 0
    c_out = 0
    c_chr_filter = 0
    c_dupl_filter = 0
    reg_len_sum = 0
    if reg2sc_dic is None:
        reg2sc_dic = {}

    with open(in_bed) as f:
        for line in f:
            cols = line.strip().split("\t")
            chr_id = cols[0]
            reg_s = int(cols[1])
            reg_e = int(cols[2])
            reg_id = cols[3]
            reg_sc = float(cols[score_col-1])
            reg_pol = cols[5]

            assert len(cols) >= score_col, "given score column for --in BED file out of range (# columns = %i, set score column = %i)" %(len(cols), score_col)

            c_in += 1

            # Check chr_id.
            if chr_ids_dic is not None:
                if chr_id not in chr_ids_dic:
                    c_chr_filter += 1
                    continue

            """
            If unstranded = True, add both strands for each region.

            """
            if unstranded:

                # Extend (use + region style extension for both).
                new_s = reg_s - ext_up
                if new_s < 0:
                    new_s = 0
                new_e = reg_e + ext_down

                if remove_dupl:
                    reg_str = "%s:%i-%i(+)" % (chr_id, new_s, new_e)
                    if reg_str not in reg2sc_dic:
                        reg2sc_dic[reg_str] = reg_sc
                    else:
                        c_dupl_filter += 1
                        continue
                    reg_str = "%s:%i-%i(-)" % (chr_id, new_s, new_e)
                    if reg_str not in reg2sc_dic:
                        reg2sc_dic[reg_str] = reg_sc
                    else:
                        c_dupl_filter += 1
                        continue
                
                """
                Add stats separately for strands.
                This could be adapted, counting only once per two-strand region.

                """
                reg_len = new_e - new_s
                reg_len_sum += reg_len*2
                c_out += 2

                OUTBED.write("%s\t%s\t%s\t%s\t%s\t+\n" % (chr_id, new_s, new_e, reg_id, reg_sc))
                OUTBED.write("%s\t%s\t%s\t%s\t%s\t-\n" % (chr_id, new_s, new_e, reg_id, reg_sc))

            else:

                assert reg_pol == "+" or reg_pol == "-", "invalid region polarity (strand) given in --in BED file %s column 6 (found \"%s\", expecting \"+\" or \"-\"). For unstranded search see --unstranded option" %(in_bed, reg_pol)

                # Extend.
                new_s = reg_s - ext_up
                new_e = reg_e + ext_down
                if reg_pol == "-":
                    new_s = reg_s - ext_down
                    new_e = reg_e + ext_up
                # Lower bound check.
                if new_s < 0:
                    new_s = 0

                if remove_dupl:
                    reg_str = "%s:%i-%i(%s)" % (chr_id, new_s, new_e, reg_pol)
                    if reg_str not in reg2sc_dic:
                        reg2sc_dic[reg_str] = reg_sc
                    else:
                        c_dupl_filter += 1
                        continue

                reg_len = new_e - new_s
                reg_len_sum += reg_len

                c_out += 1
                
                OUTBED.write("%s\t%s\t%s\t%s\t%s\t%s\n" % (chr_id, new_s, new_e, reg_id, reg_sc, reg_pol))

    f.closed
    OUTBED.close()

    stats_dic = {}
    stats_dic["c_in"] = c_in
    stats_dic["c_out"] = c_out
    stats_dic["c_chr_filter"] = c_chr_filter
    stats_dic["c_dupl_filter"] = c_dupl_filter
    stats_dic["reg_len_sum"] = reg_len_sum 

    return stats_dic


################################################################################

def get_genomic_coords_from_seq_name(seq_name, motif_s, motif_e,
                                     one_based_start=False):
    """
    Get genomic coordinates for FIMO motif hit, with seq_name format:
    chr6:35575787-35575923(-)
    as well as motif start (0-based), end (1-based) and strand info.
    seq_name coordinates BED format, i.e. start 0-based, end 1-based.

    one_based_start:
        Set to get one-based start (BED format: 0-based)

    >>> seq_name = "chr1:100-200(+)"
    >>> motif_s = 10
    >>> motif_e = 20
    >>> get_genomic_coords_from_seq_name(seq_name, motif_s, motif_e, one_based_start=False)
    ['chr1', 110, 120, '+']
    >>> get_genomic_coords_from_seq_name(seq_name, motif_s, motif_e, one_based_start=True)
    ['chr1', 111, 120, '+']
    >>> seq_name = "chr1:100-200(-)"
    >>> get_genomic_coords_from_seq_name(seq_name, motif_s, motif_e, one_based_start=False)
    ['chr1', 180, 190, '-']
    >>> get_genomic_coords_from_seq_name(seq_name, motif_s, motif_e, one_based_start=True)
    ['chr1', 181, 190, '-']
    
    """

    if re.search("\w+:\d+-\d+\([+|-]\)", seq_name):
        m = re.search("(\w+):(\d+)-(\d+)\(([+|-])\)", seq_name)
        chr_id = m.group(1)
        reg_s = int(m.group(2))
        reg_e = int(m.group(3))
        strand = m.group(4)

        gen_motif_s = reg_s + motif_s
        gen_motif_e = reg_s + motif_e
        if strand == "-":
            gen_motif_s = reg_e - motif_e
            gen_motif_e = reg_e - motif_s

        if one_based_start:
            gen_motif_s += 1

        return [chr_id, gen_motif_s, gen_motif_e, strand]
    else:
        assert False, "invalid seq_name format given (%s)" %(seq_name)


################################################################################

def get_length_from_seq_name(seq_name):
    """
    Get length of genomic region, encoded in seq_name.
    seq_name format:
    chr6:35575787-35575923(-)
    Assumes genomic start is 0-based.

    >>> seq_name = "chr1:100-200(-)"
    >>> get_length_from_seq_name(seq_name)
    100

    """
    if re.search("\w+:\d+-\d+\(", seq_name):
        m = re.search("\w+:(\d+)-(\d+)\(", seq_name)
        reg_s = int(m.group(1))
        reg_e = int(m.group(2))
        return reg_e - reg_s
    else:
        assert False, "invalid seq_name format given (%s)" %(seq_name)


################################################################################

class MotifStats:
    """
    Stores motif hit stats.

    data_id, method_id, run_id, motif_db stored in RBPStats object, linked 
    by dictionary (common internal_id).
    
    hit_id : chr:s-e(+),motif_id

    """

    def __init__(self,
                 hit_id: str,
                 internal_id: str,
                 region_id = "-",
                 rbp_id = "-",
                 motif_id = "-",
                 chr_id = "chr1",
                 gen_s = 0, # 1-based genomic motif start.
                 gen_e = 0,
                 strand = "+",
                 region_s = 0, # 1-based region motif start.
                 region_e = 0,
                 region_len = 0,
                 uniq_count = 0,
                 fimo_score: Optional[float] = None,
                 fimo_pval: Optional[float] = None,
                 cms_score: Optional[float] = None,
                 cms_eval: Optional[float] = None) -> None:
        self.hit_id = hit_id
        self.internal_id = internal_id
        self.region_id = region_id
        self.rbp_id = rbp_id
        self.motif_id = motif_id
        self.chr_id = chr_id
        self.gen_s = gen_s
        self.gen_e = gen_e
        self.strand = strand
        self.region_s = region_s
        self.region_e = region_e
        self.region_len = region_len
        self.uniq_count = uniq_count
        self.fimo_score = fimo_score
        self.fimo_pval = fimo_pval
        self.cms_score = cms_score
        self.cms_eval = cms_eval


################################################################################

def read_in_motif_stats(in_file,
                        motif_stats_dic=None,
                        store_uniq_only=True):
    """
    Read in motif stats file into dictionary of MotifStats objects.
    Each object is one row of motif stats.

    store_uniq_only:
        Store unique genomic motif hits only (i.e. do not store same genomic hit 
        twice if it appears in input file).

    """

    if motif_stats_dic is None:
        motif_stats_dic = {}
    seen_int_hit_ids_dic = {}

    with open(in_file) as f:
        for line in f:
            cols = line.strip().split("\t")
            internal_id = cols[19]
            if internal_id == "internal_id":
                continue
            hit_id = "%s:%s-%s(%s),%s" %(cols[7], cols[8], cols[9], cols[10], cols[6])
            
            int_hit_id = internal_id + "," + hit_id
            if store_uniq_only:
                if int_hit_id in seen_int_hit_ids_dic:
                    continue
                else:
                    seen_int_hit_ids_dic[int_hit_id] = 1
                
            motif_stats = MotifStats(hit_id, internal_id)
            motif_stats.region_id = cols[4]
            motif_stats.rbp_id = cols[5]
            motif_stats.motif_id = cols[6]
            motif_stats.chr_id = cols[7]
            motif_stats.gen_s = int(cols[8])
            motif_stats.gen_e = int(cols[9])
            motif_stats.strand = cols[10]
            motif_stats.region_s = int(cols[11])
            motif_stats.region_e = int(cols[12])
            motif_stats.region_len = int(cols[13])
            motif_stats.uniq_count = int(cols[14])
            fimo_score = cols[15]
            fimo_pval = cols[16]
            cms_score = cols[17]
            cms_eval = cols[18]
            if fimo_score != "-":
                motif_stats.fimo_score = float(fimo_score)
                motif_stats.fimo_pval = float(fimo_pval)
            if cms_score != "-":
                motif_stats.cms_score = float(cms_score)
                motif_stats.cms_eval = float(cms_eval)           

            if internal_id in motif_stats_dic:
                motif_stats_dic[internal_id].append(motif_stats)
            else:
                motif_stats_dic[internal_id] = [motif_stats]
    f.closed


################################################################################

class RBPStats:
    """
    Stores RBP hit stats.

    data_id, method_id, run_id, motif_db stored in RBPStats object, linked 
    by dictionary (common internal_id) to MotifStats objects.
    
    """

    def __init__(self,
                 internal_id: str,
                 data_id: str,
                 method_id: str,
                 run_id: str,
                 motif_db: str,
                 rbp_id = "-",
                 c_regions = 0,
                 mean_reg_len = 0.0,
                 median_reg_len = 0,
                 min_reg_len = 0,
                 max_reg_len = 0,
                 called_reg_size = 0,
                 effective_reg_size = 0,
                 c_reg_with_hits = 0,
                 perc_reg_with_hits = 0.0,
                 c_motif_hits = 0,
                 c_uniq_motif_hits = 0,
                 c_uniq_motif_nts = 0,
                 perc_uniq_motif_nts_cal_reg = 0.0,
                 perc_uniq_motif_nts_eff_reg = 0.0,
                 uniq_motif_hits_cal_1000nt = 0.0,
                 uniq_motif_hits_eff_1000nt = 0.0,
                 wc_pval = 1.0,
                 seq_motif_ids = None,
                 str_motif_ids = None,
                 seq_motif_hits = None,
                 str_motif_hits = None) -> None:
        self.internal_id = internal_id
        self.data_id = data_id
        self.method_id = method_id
        self.run_id = run_id
        self.motif_db = motif_db
        self.rbp_id = rbp_id
        self.c_regions = c_regions
        self.mean_reg_len = mean_reg_len
        self.median_reg_len = median_reg_len
        self.min_reg_len = min_reg_len
        self.max_reg_len = max_reg_len
        self.called_reg_size = called_reg_size
        self.effective_reg_size = effective_reg_size
        self.c_reg_with_hits = c_reg_with_hits
        self.perc_reg_with_hits = perc_reg_with_hits
        self.c_motif_hits = c_motif_hits
        self.c_uniq_motif_hits = c_uniq_motif_hits
        self.c_uniq_motif_nts = c_uniq_motif_nts
        self.perc_uniq_motif_nts_cal_reg = perc_uniq_motif_nts_cal_reg
        self.perc_uniq_motif_nts_eff_reg = perc_uniq_motif_nts_eff_reg
        self.uniq_motif_hits_cal_1000nt = uniq_motif_hits_cal_1000nt
        self.uniq_motif_hits_eff_1000nt = uniq_motif_hits_eff_1000nt
        self.wc_pval = wc_pval
        if seq_motif_ids is None:
            self.seq_motif_ids = []
        else:
            self.seq_motif_ids = seq_motif_ids
        if str_motif_ids is None:
            self.str_motif_ids = []
        else:
            self.str_motif_ids = str_motif_ids
        if seq_motif_hits is None:
            self.seq_motif_hits = []
        else:
            self.seq_motif_hits = seq_motif_hits
        if str_motif_hits is None:
            self.str_motif_hits = []
        else:
            self.str_motif_hits = str_motif_hits


################################################################################

def read_in_rbp_stats(in_file,
                      rbp_stats_dic=None):
    """
    Read in RBP stats file into dictionary of RBPStats objects.
    Each object is one row of RBP stats.

    """

    if rbp_stats_dic is None:
        rbp_stats_dic = {}
    id_check_dic = {}
    with open(in_file) as f:
        for line in f:
            cols = line.strip().split("\t")
            internal_id = cols[26]
            if internal_id == "internal_id":
                continue
            assert internal_id not in id_check_dic, "internal_id %s (supposed to be unique) appears > 1 in %s. Please contact developers!" %(internal_id, in_file)
            id_check_dic[internal_id] = 1
            rbp_stats = RBPStats(internal_id, cols[0], cols[1], cols[2], cols[3])
            rbp_stats.rbp_id = cols[4]
            rbp_stats.c_regions = int(cols[5])
            rbp_stats.mean_reg_len = float(cols[6])
            rbp_stats.median_reg_len = int(cols[7])
            rbp_stats.min_reg_len = int(cols[8])
            rbp_stats.max_reg_len = int(cols[9])
            rbp_stats.called_reg_size = int(cols[10])
            rbp_stats.effective_reg_size = int(cols[11])
            rbp_stats.c_reg_with_hits = int(cols[12])
            rbp_stats.perc_reg_with_hits = float(cols[13])
            rbp_stats.c_motif_hits = int(cols[14])
            rbp_stats.c_uniq_motif_hits = int(cols[15])
            rbp_stats.c_uniq_motif_nts = int(cols[16])
            rbp_stats.perc_uniq_motif_nts_cal_reg = float(cols[17])
            rbp_stats.perc_uniq_motif_nts_eff_reg = float(cols[18])
            rbp_stats.uniq_motif_hits_cal_1000nt = float(cols[19])
            rbp_stats.uniq_motif_hits_eff_1000nt = float(cols[20])
            rbp_stats.wc_pval = float(cols[21])
            seq_motif_ids = cols[22]
            seq_motif_hits = cols[23]
            str_motif_ids = cols[24]
            str_motif_hits = cols[25]
            if seq_motif_ids != "-":
                for motif_id in seq_motif_ids.split(","):
                    rbp_stats.seq_motif_ids.append(motif_id)
                for hit_c in seq_motif_hits.split(","):
                    rbp_stats.seq_motif_hits.append(int(hit_c))
            if str_motif_ids != "-":
                for motif_id in str_motif_ids.split(","):
                    rbp_stats.str_motif_ids.append(motif_id)
                for hit_c in str_motif_hits.split(","):
                    rbp_stats.str_motif_hits.append(int(hit_c))
            assert internal_id not in rbp_stats_dic, "internal_id %s appears > 1 in provided input files. Please provide input files with unique internal IDs"
            rbp_stats_dic[internal_id] = rbp_stats
    f.closed


################################################################################

class RBP:
    """
    RBP class.

    Note that all motifs of an RBP are counting (not just single ones).

    """
    def __init__(self, name: str,
                 internal_id: str,
                 seq_motif_ids = None,
                 str_motif_ids = None,
                 seq_motif_hits = None,
                 str_motif_hits = None,
                 c_hit_reg = 0, # # regions with motif hits.
                 c_no_hit_reg = 0, # # region with no motif hits.
                 perc_hit_reg = 0.0, # % hit regions over all regions (i.e. how many input regions contain >= 1 RBP motif).
                 c_motif_hits = 0, # # motif hits.
                 c_uniq_motif_hits = 0, # # unique motif hits.
                 c_uniq_motif_nts = 0, # # unique motif nucleotides.
                 perc_uniq_motif_nts_eff_reg = 0.0, # % unique motif nts over effective region length.
                 perc_uniq_motif_nts_cal_reg = 0.0, # % unique motif nts over called region length.
                 uniq_motif_hits_eff_1000nt = 0.0, # unique motif hits per effective 1000 nt.
                 uniq_motif_hits_cal_1000nt = 0.0, # unique motif hits per called 1000 nt.
                 # ks_pval = 1.0, # Kolmogorov-Smirnov (KS) statistic p-value (are higher scoring sites enriched with motifs).
                 wc_pval = 1.0, # Wilcoxon rank-sum test (Mann-Whitney U test) statistic p-value (are higher scoring sites enriched with motifs).
                 wc_pval_less = 1.0, # Wilcoxon rank-sum test (Mann-Whitney U test) statistic p-value (are lower scoring sites enriched with motifs).
                 organism: Optional[str] = None) -> None:
        self.name = name
        self.internal_id = internal_id
        if seq_motif_ids is None:
            self.seq_motif_ids = []
        else:
            self.seq_motif_ids = seq_motif_ids
        if str_motif_ids is None:
            self.str_motif_ids = []
        else:
            self.str_motif_ids = str_motif_ids
        if seq_motif_hits is None:
            self.seq_motif_hits = []
        else:
            self.seq_motif_hits = seq_motif_hits
        if str_motif_hits is None:
            self.str_motif_hits = []
        else:
            self.str_motif_hits = str_motif_hits
        self.c_hit_reg = c_hit_reg
        self.c_no_hit_reg = c_no_hit_reg
        self.perc_hit_reg = perc_hit_reg
        self.c_motif_hits = c_motif_hits
        self.c_uniq_motif_hits = c_uniq_motif_hits
        self.c_uniq_motif_nts = c_uniq_motif_nts
        self.perc_uniq_motif_nts_eff_reg = perc_uniq_motif_nts_eff_reg
        self.perc_uniq_motif_nts_cal_reg = perc_uniq_motif_nts_cal_reg
        self.uniq_motif_hits_eff_1000nt = uniq_motif_hits_eff_1000nt
        self.uniq_motif_hits_cal_1000nt = uniq_motif_hits_cal_1000nt
        # self.ks_pval = ks_pval
        self.wc_pval = wc_pval
        self.wc_pval_less = wc_pval_less
        self.organism = organism


################################################################################

class GenomicRegion:
    """
    Genomic region class, with 1-based start+end coordinates.

    """
    def __init__(self, chr_id: str, start: int, end: int, strand: str,
                 score: Optional[float] = None, 
                 name: Optional[str] = None, 
                 genome: Optional[str] = None) -> None:
        self.chr_id = chr_id
        self.start = start
        self.end = end
        self.name = name
        self.score = score
        self.strand = strand
        self.genome = genome

    def __repr__(self) -> str:
        return f"{self.chr_id}:{self.start}-{self.end}({self.strand})"

    def __eq__(self, other) -> bool:
        if not isinstance(other, GenomicRegion):
            return False
        return (self.chr_id == other.chr_id and
                self.start == other.start and
                self.end == other.end and
                self.strand == other.strand and
                self.genome == other.genome)

    def __len__(self):
        return self.end - self.start + 1

    def overlaps(self, other) -> bool:
        if self.chr_id != other.chr_id:
            return False
        if self.strand != other.strand:
            return False
        return self.start <= other.end and self.end >= other.start


################################################################################

class FimoHit(GenomicRegion):
    """
    Fimo motif hit class, with 1-based start+end coordinates.
    
    """
    def __init__(self, chr_id, start, end, strand, score,
                 motif_id: str, seq_name: str, pval: float, 
                 qval: Optional[float] = None,
                 matched_seq: Optional[str] = None,
                 seq_s: Optional[int] = None,
                 seq_e: Optional[int] = None,
                 genome: Optional[str] = None) -> None:
        super().__init__(chr_id, start, end, strand, score, genome)
        self.motif_id = motif_id
        self.seq_name = seq_name
        self.pval = pval
        self.qval = qval
        self.matched_seq = matched_seq
        self.seq_s = seq_s
        self.seq_e = seq_e

    def __eq__(self, other) -> bool:
        if not isinstance(other, FimoHit):
            return False
        return (self.chr_id == other.chr_id and
                self.start == other.start and
                self.end == other.end and
                self.strand == other.strand and
                self.motif_id == other.motif_id and
                self.genome == other.genome)

    def __repr__(self) -> str:
        return f"{self.chr_id}:{self.start}-{self.end}({self.strand}),{self.motif_id}"


################################################################################

def read_in_fimo_results(fimo_tsv,
                         fast_fimo=True):
    """
    Read in FIMO motif finding results (TSV file).
    
    Note: sequence_name coordinates BED format, i.e. start 0-based, end 1-based.

    Example output "Nothing found":
    # FIMO (Find Individual Motif Occurrences): Version 5.4.1 compiled on Aug  1 2022 at 17:19:30
    # The format of this file is described at https://meme-suite.org/meme/doc/fimo-output-format.html.
    # fimo --norc --verbosity 1 -oc test_out test_out/513dbc80-cef8-11ed-95a7-901b0eb924fa.meme_motifs.xml test_out/513dbc44-cef8-11ed-95a7-901b0eb924fa.in_sites.fa

    Example output "Something found":

    motif_id	motif_alt_id	sequence_name	start	stop	strand	score	p-value	q-value	matched_sequence
    AGGF1_2		chr6:35575787-35575923(-)	6	12	+	10.4412	0.000093	0.00632	AAGAAGA
    AGGF1_2		chr6:35575788-35575924(-)	7	13	+	10.4412	0.000093	0.00632	AAGAAGA
    AGGF1_2		chr16:69653258-69653362(+)	9	15	+	10.4412	0.000093	0.00632	AAGAAGA
    AGGF1_2		chr16:69653257-69653361(+)	10	16	+	10.4412	0.000093	0.00632	AAGAAGA
    AGGF1_2		chr2:24067815-24067851(-)	22	28	+	10.4412	0.000093	0.00632	AAGAAGA
    AGGF1_2		chr6:35575787-35575923(-)	61	67	+	10.4412	0.000093	0.00632	AAGAAGA
    AGGF1_2		chr6:35575788-35575924(-)	62	68	+	10.4412	0.000093	0.00632	AAGAAGA

    # FIMO (Find Individual Motif Occurrences): Version 5.4.1 compiled on Aug  1 2022 at 17:19:30
    # The format of this file is described at https://meme-suite.org/meme/doc/fimo-output-format.html.
    # fimo --norc --verbosity 1 -oc test_out test_out/980739c0-cef8-11ed-8c3c-901b0eb924fa.meme_motifs.xml test_out/9807398e-cef8-11ed-8c3c-901b0eb924fa.in_sites.fa
 
    fimo --norc --verbosity 1 --skip-matched-sequence --text test_pum2.meme test_pum2.fa > check.tsv
    motif_id	motif_alt_id	sequence_name	start	stop	strand	score	p-value	q-value	matched_sequence
    PUM2_2		chr19:100000-100085(+)	21	28	+	15.7959	0.0000171		
    PUM2_2		chr18:100000-100096(+)	17	24	+	15.7959	0.0000171		
    PUM2_2		chr18:100000-100096(+)	75	82	+	15.7959	0.0000171    

    UPDATES:
    Do not read in q-value or matched_sequence, as these are not produced when 
    using run_fast_fimo().

    """

    fimo_hits_list = []

    with open(fimo_tsv) as f:
        for line in f:
            if re.search("^#", line):
                continue
            cols = line.strip().split("\t")
            if cols[0] == "motif_id" or cols[0] == "":
                continue

            motif_id = cols[0]
            seq_name = cols[2]
            motif_s = int(cols[3])-1 # make 0-based.
            motif_e = int(cols[4])
            score = float(cols[6])
            pval = float(cols[7])
            qval = None
            matched_seq = None
            if not fast_fimo:
                qval = float(cols[8])
                matched_seq = cols[9]

            gen_motif_coords = get_genomic_coords_from_seq_name(seq_name, motif_s, motif_e,
                                                                one_based_start=True)
            
            fimo_hit = FimoHit(chr_id=gen_motif_coords[0], 
                               start=gen_motif_coords[1], 
                               end=gen_motif_coords[2],
                               strand=gen_motif_coords[3], 
                               score=score, 
                               motif_id=motif_id, 
                               seq_name=seq_name, 
                               pval=pval, 
                               qval=qval,
                               seq_s=motif_s+1,
                               seq_e=motif_e,
                               matched_seq=matched_seq)

            fimo_hits_list.append(fimo_hit)

    f.closed

    return fimo_hits_list


################################################################################

def output_motif_hits_to_bed(rbp_id, unique_motifs_dic, out_bed,
                             one_based_start=True):
    """
    Output unique motif hit regions to out_bed BED, stored in 
    unique_motifs_dic[rbp_id] for given RBP ID rbp_id.

    one_based_start:
        Assumes FIMO hit string contains 1-based start.

    Format of fimo_hit_string:
    chr_id:start-end(strand),motif_id

    """
    
    assert rbp_id in unique_motifs_dic, "rbp_id %s not in unique_motifs_dic" %(rbp_id)

    OUTMRBED = open(out_bed, "w")

    for fh_str in unique_motifs_dic[rbp_id]:
        if re.search("\w+:\d+-\d+\([+|-]\),.+", fh_str):
            m = re.search("(\w+):(\d+)-(\d+)\(([+|-])\),(.+)", fh_str)
            chr_id = m.group(1)
            reg_s = int(m.group(2))
            reg_e = int(m.group(3))
            strand = m.group(4)
            motif_id = m.group(5)

            if one_based_start:
                reg_s -= 1

            OUTMRBED.write("%s\t%i\t%i\t%s\t0\t%s\n" % (chr_id, reg_s, reg_e, motif_id, strand))
        else:
            assert False, "invalid fh_str format given (%s)" %(fh_str)

    OUTMRBED.close()


################################################################################

def batch_output_motif_hits_to_bed(unique_motifs_dic, out_bed,
                                   one_based_start=True):
    """
    Output unique motif hit regions to out_bed BED, stored in 
    unique_motifs_dic.

    one_based_start:
        Assumes FIMO hit string contains 1-based start.

    Format of fimo_hit_string:
    chr_id:start-end(strand),motif_id

    """
    
    OUTMRBED = open(out_bed, "w")

    for fh_str in unique_motifs_dic:
        if re.search("\w+:\d+-\d+\([+|-]\),.+", fh_str):
            m = re.search("(\w+):(\d+)-(\d+)\(([+|-])\),(.+)", fh_str)
            chr_id = m.group(1)
            reg_s = int(m.group(2))
            reg_e = int(m.group(3))
            strand = m.group(4)
            motif_id = m.group(5)

            if one_based_start:
                reg_s -= 1

            OUTMRBED.write("%s\t%i\t%i\t%s\t0\t%s\n" % (chr_id, reg_s, reg_e, motif_id, strand))
        else:
            assert False, "invalid fh_str format given (%s)" %(fh_str)

    OUTMRBED.close()


################################################################################

def read_in_cmsearch_results(in_tab,
                             check=True,
                             hits_list=None):
    """
    Read in cmsearch motif finding results file.
    
    Output columns:
    1 target name
    2 accession (taget)
    3 query name
    4 accession (query) -> RFAM ID, use this as motif identifer.
    5 mdl
    6 mdl from
    7 mdl to
    8 seq from
    9 seq to
    10 strand (for --toponly should always be +)
    11 trunc
    12 pass
    13 gc
    14 bias
    15 score
    16 E-value
    17 inc
    18 description of target

    """

    if check:
        assert os.path.exists(in_tab), "cmsearch output file %s does not exist! Check cmsearch output for errors (likely a failed cmsearch call due to invalid parameter settings)" % (in_tab)

    if hits_list is None:
        hits_list = []

    c_hits = 0

    with open(in_tab) as f:
        for line in f:
            if re.search("^#", line):
                continue
            cols_pre = line.strip().split(" ")
            # Remove empty column values.
            cols = []
            for c in cols_pre:
                if c:
                    cols.append(c)
            c_cols = len(cols)
            assert c_cols == 18, "invalid cmsearch results table row encountered (# columns != 18 (%i)). Problem row is:\n%s" %(c_cols, str(cols))

            seq_name = cols[0] # FASTA sequence header.
            query_name = cols[2]
            motif_id = cols[3] # use query accession ID (RFAM ID) as motif ID.
            model = cols[4]
            model_s = int(cols[5])
            model_e = int(cols[6])
            seq_s = int(cols[7]) - 1 # make 0-based.
            seq_e = int(cols[8])
            strand = cols[9]
            score = float(cols[14])
            e_value = float(cols[15])

            assert strand == "+", "invalid cmsearch results table row encountered. strand (column 10) expected to be +, but found:\n%s" %(str(cols))

            gen_motif_coords = get_genomic_coords_from_seq_name(seq_name, seq_s, seq_e,
                                                                one_based_start=True)

            # print("cols:", cols)
            # print("e_value:", e_value)

            cmsearch_hit = CmsearchHit(chr_id=gen_motif_coords[0], 
                                start=gen_motif_coords[1], 
                                end=gen_motif_coords[2],
                                strand=gen_motif_coords[3], 
                                score=score,
                                motif_id=motif_id,
                                query_name=query_name,
                                seq_name=seq_name,
                                e_value=e_value,
                                model=model,
                                seq_s=seq_s+1,
                                seq_e=seq_e,
                                model_s=model_s,
                                model_e=model_e)

            # print(cmsearch_hit.__dict__)
            c_hits += 1
            hits_list.append(cmsearch_hit)

    f.closed

    return hits_list, c_hits


################################################################################

class CmsearchHit(GenomicRegion):
    """
    cmsearch motif hit class, with 1-based start+end coordinates.
    
    """
    def __init__(self, chr_id, start, end, strand, score,
                 motif_id: str, seq_name: str, e_value: float, 
                 query_name: Optional[str] = None,
                 model: Optional[str] = None,
                 model_s: Optional[int] = None,
                 model_e: Optional[int] = None,
                 seq_s: Optional[int] = None,
                 seq_e: Optional[int] = None,
                 genome: Optional[str] = None) -> None:
        super().__init__(chr_id, start, end, strand, score, genome)
        self.motif_id = motif_id
        self.seq_name = seq_name
        self.e_value = e_value
        self.model = model
        self.model_s = model_s
        self.model_e = model_e
        self.seq_s = seq_s
        self.seq_e = seq_e

    def __eq__(self, other) -> bool:
        if not isinstance(other, CmsearchHit):
            return False
        return (self.chr_id == other.chr_id and
                self.start == other.start and
                self.end == other.end and
                self.strand == other.strand and
                self.motif_id == other.motif_id and
                self.genome == other.genome)

    def __repr__(self) -> str:
        return f"{self.chr_id}:{self.start}-{self.end}({self.strand}),{self.motif_id}"


################################################################################

def remove_special_chars_from_str(check_str,
                                  reg_ex='[^A-Za-z0-9_-]+'):
    """
    Remove special characters from string.

    reg_ex:
        Regular expression defining what to keep.

    >>> check_str = "{_}[-](_)\V/"
    >>> remove_special_chars_from_str(check_str)
    '_-_V'
    >>> check_str = ""
    >>> remove_special_chars_from_str(check_str)
    ''

    """
    clean_string = re.sub(reg_ex, '', check_str)
    return clean_string


################################################################################

def get_motif_id_from_str_repr(hit_str_repr):
    """
    From motif string representation:
    chr1:100-200(-),motif_id
    return motif_id

    >>> hit_str_repr = "chr6:66-666(-),satan6666"
    >>> get_motif_id_from_str_repr(hit_str_repr)
    'satan6666'

    """

    if re.search("\w+:\d+-\d+\([+|-]\),.+", hit_str_repr):
        m = re.search(".+,(.+)", hit_str_repr)
        motif_id = m.group(1)
        return motif_id
    else:
        assert False, "motif_id extraction failed due to invalid hit_str_repr format given (%s)" %(hit_str_repr)


################################################################################

def search_generate_html_report(df_corr, df_pval, pval_cont_lll,
                                search_rbps_dic,
                                fimo_hits_list, cmsearch_hits_list,
                                id2name_dic, out_folder, 
                                benchlib_path,
                                rbp2regidx_dic,
                                c_regions,
                                upset_plot_min_degree=2,
                                upset_plot_min_subset_size=10,
                                html_report_out="report.rbpbench_search.html",
                                plot_abs_paths=False,
                                plots_subfolder="html_report_plots"):
    """
    Create additional hit statistics for selected RBPs, 
    e.g. correlation / co-occurrence between RBPs.

    """
    # Use absolute paths?
    if plot_abs_paths:
        out_folder = os.path.abspath(out_folder)
    
    plots_folder = plots_subfolder
    plots_out_folder = out_folder + "/" + plots_folder
    if plot_abs_paths:
        plots_folder = plots_out_folder

    if not os.path.exists(plots_out_folder):
        os.makedirs(plots_out_folder)
    html_out = out_folder + "/" + "report.rbpbench_search.html"
    md_out = out_folder + "/" + "report.rbpbench_search.md"
    if html_report_out:
        html_out = html_report_out

    # Makes tables sortable.
    sorttable_js_path = benchlib_path + "/content/sorttable.js"
    # plotly js path.
    # https://plotly.com/javascript/getting-started/#download
    # plotly-2.20.0.min.js
    # plotly-latest.min.js
    plotly_js_path = benchlib_path + "/content/plotly-2.20.0.min.js"
    assert os.path.exists(plotly_js_path), "plotly js %s not found" %(plotly_js_path)

    # Create theme-specific HTML header.
    mdtext = """
<head>
<title>RBPBench - Search Report</title>
<script src="%s" type="text/javascript"></script>
</head>
""" %(sorttable_js_path)

    # Add first section markdown.
    mdtext += """

# Search report

List of available statistics and plots generated
by RBPBench (rbpbench search --report):

- [RBP motif enrichment statistics](#rbp-enrich-stats)
- [RBP co-occurrences heat map](#cooc-heat-map)
- [RBP correlations heat map](#corr-heat-map)
- [RBP combinations upset plot](#rbp-comb-upset-plot)"""

    mdtext += "\n&nbsp;\n"

    """
    RBP motif enrichment statistics

    """

    c_in_regions = 0
    for rbp_id in search_rbps_dic:
        c_in_regions += search_rbps_dic[rbp_id].c_hit_reg
        c_in_regions += search_rbps_dic[rbp_id].c_no_hit_reg      
        break

    mdtext += """
## RBP motif enrichment statistics ### {#rbp-enrich-stats}

**Table:** RBP motif enrichment statistics. Given a score for each genomic region (# input regions = %i), 
RBPbench checks whether motifs are enriched 
in higher-scoring regions (using Wilcoxon rank-sum test). A low Wilcoxon rank-sum test p-value for a given RBP thus indicates 
that higher-scoring regions are more likely to contain motif hits of the respective RBP. NOTE that if scores associated to 
input genomic regions are all the same, p-values become meaningless (i.e., they result in p-values of 1.0).


""" %(c_in_regions)
    mdtext += '| RBP ID | # hit regions | % hit regions | # motif hits | p-value |' + " \n"
    mdtext += "| :-: | :-: | :-: | :-: | :-: |\n"
    for rbp_id in search_rbps_dic:
        wc_pval = search_rbps_dic[rbp_id].wc_pval
        c_hit_reg = search_rbps_dic[rbp_id].c_hit_reg
        perc_hit_reg = search_rbps_dic[rbp_id].perc_hit_reg
        c_uniq_motif_hits = search_rbps_dic[rbp_id].c_uniq_motif_hits
        mdtext += "| %s | %i | %.2f | %i | %s |\n" %(rbp_id, c_hit_reg, perc_hit_reg, c_uniq_motif_hits, str(wc_pval))
    mdtext += "\n&nbsp;\n&nbsp;\n"
    mdtext += "\nColumn IDs have the following meanings: "
    mdtext += "**RBP ID** -> RBP ID from database or user-defined (typically RBP name), "
    mdtext += '**# hit regions** -> number of input genomic regions with motif hits (after filtering and optional extension), '
    mdtext += '**% hit regions** -> percentage of hit regions over all regions (i.e. how many input regions contain >= 1 RBP binding motif), '
    mdtext += '**# motif hits** -> number of unique motif hits in input regions (removed double counts), '
    mdtext += '**p-value** -> Wilcoxon rank-sum test p-value.' + "\n"

    """
    Co-occurrence heat map.

    """

    cooc_plot_plotly =  "co-occurrence_plot.plotly.html"
    cooc_plot_plotly_out = plots_out_folder + "/" + cooc_plot_plotly

    create_cooc_plot_plotly(df_pval, pval_cont_lll,
                            plotly_js_path, cooc_plot_plotly_out)

    plot_path = plots_folder + "/" + cooc_plot_plotly

    mdtext += """
## RBP co-occurrences heat map ### {#cooc-heat-map}

RBP co-occurrences heat map.

"""
    mdtext += '<div class=class="container-fluid" style="margin-top:40px">' + "\n"
    mdtext += '<iframe src="' + plot_path + '" width="1200" height="1200"></iframe>' + "\n"
    mdtext += '</div>'
    mdtext += """

**Figure:** Heat map of co-occurrences (Fisher's exact test p-values) between RBPs. 
Legend color: negative logarithm (base 10) of Fisher's exact test p-value.
Hover box: 1) RBP1. 2) RBP2. 3) p-value: Fisher's exact test p-value (calculated based on contingency table between RBP1 and RBP2). 
4) RBPs compaired. 5) Counts[]: Contingency table of co-occurrence counts (i.e., number of genomic regions with/without shared motif hits) between compaired RBPs, 
with format [[A, B], [C, D]], where 
A: RBP1 AND RBP2, 
B: NOT RBP1 AND RBP2
C: RBP1 AND NOT RBP2
D: NOT RBP1 AND NOT RBP2. 

&nbsp;

"""

    """
    Correlation heat map.

    """

    corr_plot_plotly =  "correlation_plot.plotly.html"
    corr_plot_plotly_out = plots_out_folder + "/" + corr_plot_plotly

    create_corr_plot_plotly(df_corr, pval_cont_lll,
                            plotly_js_path, corr_plot_plotly_out)

    plot_path = plots_folder + "/" + corr_plot_plotly

    mdtext += """
## RBP correlations heat map ### {#corr-heat-map}

RBP correlations heat map.

"""

    mdtext += '<div class=class="container-fluid" style="margin-top:40px">' + "\n"
    mdtext += '<iframe src="' + plot_path + '" width="1200" height="1200"></iframe>' + "\n"
    mdtext += '</div>'
    mdtext += """

**Figure:** Heat map of correlations (Pearson correlation coefficients) between RBPs. 
Genomic regions are labelled 1 or 0 (RBP motif present or not), resulting in a vector of 1s and 0s for each RBP.
Correlations are then calculated by comparing vectors for every pair of RBPs.
Legend color: Pearson correlation coefficient. 
Hover box: 1) RBP1. 2) RBP2.
3) RBPs compaired. 5) Counts[]: Contingency table of co-occurrence counts (i.e., number of genomic regions with/without shared motif hits) between compaired RBPs, 
with format [[A, B], [C, D]], where 
A: RBP1 AND RBP2, 
B: NOT RBP1 AND RBP2
C: RBP1 AND NOT RBP2
D: NOT RBP1 AND NOT RBP2. 

&nbsp;

"""



    """
    RBP region occupancies upset plot.

    """

    rbp_reg_occ_upset_plot =  "rbp_region_occupancies.upset_plot.png"
    rbp_reg_occ_upset_plot_out = plots_out_folder + "/" + rbp_reg_occ_upset_plot

    plotted, reason = create_rbp_reg_occ_upset_plot(rbp2regidx_dic, c_regions,
                                  min_degree=upset_plot_min_degree,
                                  min_subset_size=upset_plot_min_subset_size,
                                  plot_out=rbp_reg_occ_upset_plot_out)


    plot_path = plots_folder + "/" + rbp_reg_occ_upset_plot

    mdtext += """
## RBP combinations upset plot ### {#rbp-comb-upset-plot}

RBP combinations upset plot.

"""

    if plotted:
        mdtext += '<img src="' + plot_path + '" alt="RBP region occupancies upset plot"' + "\n"
        mdtext += 'title="RBP region occupancies upset plot" width="600" />' + "\n"
        mdtext += """

**Figure:** Upset plot of RBP combinations found in the given set of genomic regions (# of regions = %i). 
Intersection size == how often a specific RBP combination is found in the regions dataset.
For example, if two regions in the input set contain motif hits for RBP1, RBP3, and RBP5, then the RBP combination RBP1,RBP3,RBP5 will get a count (== Intersection size) of 2.
Minimum occurrence number for a combination to be reported = %i (command line parameter: --upset-plot-min-subset-size). 
How many RBPs a combination must contain to be reported (--up) = %i (command line parameter: --upset-plot-min-degree). 
Number of

&nbsp;

""" %(c_regions, upset_plot_min_subset_size, upset_plot_min_degree)

    else:

        if reason == "min_degree":

            mdtext += """

No upset plot generated since set --upset-plot-min-degree > maximum degree found in the RBP combination set. Please use lower number for --upset-plot-min-degree parameter.

&nbsp;

"""

        elif reason == "min_degree":

            mdtext += """

No upset plot generated since set --upset-plot-min-subset-size > maximum subset size found in the RBP combination set. Please use lower number for --upset-plot-min-subset-size parameter.

&nbsp;

"""
        else:
            assert False, "invalid reason given for not plotting upset plot"


    # ALAMO

    # Convert mdtext to html.
    md2html = markdown(mdtext, extensions=['attr_list', 'tables'])

    OUTMD = open(md_out,"w")
    OUTMD.write("%s\n" %(mdtext))
    OUTMD.close()
    OUTHTML = open(html_out,"w")
    OUTHTML.write("%s\n" %(md2html))
    OUTHTML.close()

    # change <table> to sortable.
    check_cmd = "sed -i 's/<table>/<table class=" + '"sortable"' + ">/g' " + html_out
    output = subprocess.getoutput(check_cmd)
    error = False
    if output:
        error = True
    assert error == False, "sed command returned error:\n%s" %(output)


################################################################################

def create_rbp_reg_occ_upset_plot(rbp2regidx_dic, c_regions,
                                  min_degree=2,
                                  min_subset_size=10,
                                  plot_out="rbp_region_occupancies.upset_plot.png"):
    """
    Create upset plot for RBP region occupancies.

    conda install -c conda-forge upsetplot
    import pandas as pd
    from upsetplot import plot
    from matplotlib import pyplot

    Return False if plotting not possible (due to set min_degree or min_subset_size too high).

    """


    # All regions (# regions).
    assert c_regions > 0, "c_regions must be >= 0"
    rbp_idx_list = []
    for i in range(c_regions):
        rbp_idx_list.append(i)

    # Convert to sets.
    rbp_id_list = []
    for rbp_id, rbp_list in sorted(rbp2regidx_dic.items()):
        rbp_set = set(rbp_list)
        rbp2regidx_dic[rbp_id] = rbp_set
        rbp_id_list.append(rbp_id)
    rbp_idx_list = set(rbp_idx_list)

    # Create boolean list of lists to convert to pandas DataFrame.
    bool_ll = []
    for idx in rbp_idx_list:
        bool_l = []
        for rbp_id, rbp_set in sorted(rbp2regidx_dic.items()):
            if idx in rbp_set:
                bool_l.append(True)
            else:
                bool_l.append(False)        
        bool_ll.append(bool_l)

    df = pd.DataFrame(bool_ll, columns=rbp_id_list)

    df_up = df.groupby(rbp_id_list).size()

    # Check if set min_degree and min_subset_size are too high (to prevent plotting error).
    df_up_check = df_up[df_up.index.map(sum) >= min_degree]
    if df_up_check.empty:
        # Set min_degree too high.
        return False, "min_degree"
    max_subset_size = 0
    for elem_c in df_up_check:
        if elem_c > max_subset_size:
            max_subset_size = elem_c

    if max_subset_size < min_subset_size:
        # Set min_subset_size too high.
        return False, "min_subset_size"

    # Move on to plotting.
    upset_plot = plot(df_up, orientation='horizontal', 
                    min_degree=min_degree,  # number of RBPs in set (e.g. 2 -> at least 2 RBP pairs, not single RBPs)
                    min_subset_size=min_subset_size,  # min size of a set to be reported.
                    sort_by="cardinality")

    plt.savefig(plot_out, dpi=150)
    plt.close()
    return True, "yowza"


################################################################################

def create_corr_plot_plotly(df, pval_cont_lll, plotly_js_path, plot_out):
    """
    Plot correlations as heat map with plotly.


    ax_labels_dic = {rank_df : x_label, ws_sc_df : y_label}

    plot = px.scatter(data_frame=df, x=rank_df, y=ws_sc_df, hover_name=seq_id_df,
                      labels=ax_labels_dic,
                      hover_data=[sal_peak_pos_df, sal_peak_sc_df, sal_win_coords_df, sal_win_seq_df, seq_df],
                      color_discrete_sequence=[dot_col])
    """

    plot = px.imshow(df)
    # plot.show()
    # plot.layout.template = 'seaborn'
    # plot.update_layout(hoverlabel=dict(font_size=11))
    # plot.update_traces(marker=dict(size=3))
    # plot.update_traces(marker=dict(size=3))
    #fig.update_traces(hovertemplate='GDP: %{x} <br>Life Expectancy: %{y}')

    plot.update(data=[{'customdata': pval_cont_lll,
                    'hovertemplate': 'RBP1: %{x}<br>RBP2: %{y}<br>RBPs: %{customdata[1]}<br>Counts: %{customdata[2]}<br>Correlation: %{z}<extra></extra>'}])
    # plot.update(data=[{'customdata': pval_cont_lll,
    #                 'hovertemplate': 'RBP1: %{x}<br>RBP2: %{y}<br>RBPs: %{customdata[1]}<br>Counts: %{customdata[2]}<br>Color: %{z}<extra></extra>'}])
    # plot.update_layout(plot_bgcolor='white', legend_title_text='Correlation')
    # plot.update_traces(colorbar_title='Correlation')
    # plot.update_traces(marker=dict(colorbar={"title": "Correlation"}))
    plot.update_layout(plot_bgcolor='white')
    plot.write_html(plot_out,
                    full_html=False,
                    include_plotlyjs=plotly_js_path)
    

################################################################################

def create_cooc_plot_plotly(df, pval_cont_lll, plotly_js_path, plot_out):
    """
    Plot co-occurrences as heat map with plotly.

    """

    plot = px.imshow(df)
    plot.update(data=[{'customdata': pval_cont_lll,
                    'hovertemplate': 'RBP1: %{x}<br>RBP2: %{y}<br>p-value: %{customdata[0]}<br>RBPs: %{customdata[1]}<br>Counts: %{customdata[2]}<br>-log10(p-value): %{z}<extra></extra>'}])
    plot.update_layout(plot_bgcolor='white')
    plot.write_html(plot_out,
                    full_html=False,
                    include_plotlyjs=plotly_js_path)


################################################################################

def log_tf_pval(pval):
    """

    Log transform p-values:
    Returns negative logarithm (base 10) of given p-value pval.

    >>> pv = 0.1
    >>> log_tf_pval(pv)
    1.0
    
    """

    neg_log_pval = -log10(pval)
    # return abs(neg_log_pval)
    return neg_log_pval


################################################################################

def log_tf_df(df,
              min_pv=2.2e-308,
              convert_zero_pv=False,
              rbp_list=None):
    """
    Log transform quadratic dataframe of p-values.

    convert_zero_pv:
        If true, converts p-values of 0.0 to smallest float() class number 
        possible (i.e., 2.2e-308). Otherwise -log10(0.0) will cause
        math domain error.

    """

    if rbp_list is None:
        # Old .iloc way of addressing (gives warning from pandas 2.1.0 on).
        for i in range(len(df)):
            for j in range(len(df)):
                if df.iloc[i][j] is not None:
                    pv = df.iloc[i][j]
                    if convert_zero_pv:
                        if pv == 0:
                            pv = min_pv
                    ltf_pval = log_tf_pval(pv)
                    df.iloc[i][j] = ltf_pval
    else:
        # Way to deal with pandas 2.1.0 deprecation warning "treating keys as positions is deprecated ...".
        assert len(df) == len(rbp_list), "len(df) != len(rbp_list) (%i != %i)" %(len(df), len(rbp_list)) 
        for i,rbp_i in enumerate(rbp_list):
            for j,rbp_j in enumerate(rbp_list):
                if df.loc[rbp_i][rbp_j] is not None:
                    pv = df.loc[rbp_i][rbp_j]
                    if convert_zero_pv:
                        if pv == 0:
                            pv = min_pv
                    ltf_pval = log_tf_pval(pv)
                    df.loc[rbp_i][rbp_j] = ltf_pval


################################################################################

def make_pos_freq_matrix(seqs_dic, 
                         exp_len=21,
                         report=True,
                         to_ppm=False):
    """
    Make position frequency matrix (PFM, DNA alphabet with N), based on given 
    dictionary of sequences.

    exp_len:
        Only process sequences with length exp_len.
    to_ppm:
        Normalize frequencies to probabilities, creating a position 
        probability matrix (PPM) out of PFM.

    >>> seqs_dic = {'s1': 'ACT', 's2': 'ACA'}
    >>> make_pos_freq_matrix(seqs_dic, exp_len=3, report=False)
    {'A': [2, 0, 1], 'C': [0, 2, 0], 'G': [0, 0, 0], 'T': [0, 0, 1], 'N': [0, 0, 0]}
    >>> make_pos_freq_matrix(seqs_dic, exp_len=3, report=False, to_ppm=True)
    {'A': [1.0, 0.0, 0.5], 'C': [0.0, 1.0, 0.0], 'G': [0.0, 0.0, 0.0], 'T': [0.0, 0.0, 0.5], 'N': [0.0, 0.0, 0.0]}

    """

    # Init matrix.
    pos_freq_matrix = {
        'A' : [0]*exp_len,
        'C' : [0]*exp_len,
        'G' : [0]*exp_len,
        'T' : [0]*exp_len,
        'N' : [0]*exp_len
    }

    # Record sequences.
    c_trunc_filter = 0
    c_full_len = 0
    for seq_id in seqs_dic:
        seq = seqs_dic[seq_id]
        reg_len = len(seq)
        if reg_len != exp_len:
            c_trunc_filter += 1
            continue
        c_full_len += 1
        for idx, nt in enumerate(seq):
            assert nt in pos_freq_matrix, "nucleotide \"%s\" not in pos_freq_matrix. Please provide valid sequences"
            pos_freq_matrix[nt][idx] += 1

    assert c_full_len, "no full-length sequences found in seqs_dic (given expected length: %i)" %(exp_len)

    if report:
        print("# sequences with exp_len:  %i" %(c_full_len))
        print("# truncated sequences:     %i" %(c_trunc_filter))

    # Normalize.
    if to_ppm:
        for nt in pos_freq_matrix:
            for idx, pos_freq in enumerate(pos_freq_matrix[nt]):
                pos_prob = pos_freq / c_full_len
                pos_freq_matrix[nt][idx] = pos_prob

    return pos_freq_matrix


################################################################################

def plot_nt_distribution_zero_pos(ppm, ext_lr,
                                  plot_out="nt_dist_zero_pos.png"):
    """
    Plot nucleotide distribution around zero position.
    
    """

    exp_len = (ext_lr*2) + 1

    # Plot x-tick labels.
    x_idx = []
    x_idx_s = ext_lr*-1
    for idx in range(exp_len):
        x_idx.append(x_idx_s + idx)

    # Color scheme for plotting.
    color_scheme = {
        'A' : [0, .5, 0],
        'C' : [0, 0, 1],
        'G' : [1, .65, 0],
        'T' : [1, 0, 0],
        'N': 'gray'
    }

    fig_x = ext_lr
    fix_y = 2.25

    ppm_df = pd.DataFrame(ppm)
    if fig_x < 5:
        logo = Logo(ppm_df, color_scheme=color_scheme)
    else:
        logo = Logo(ppm_df, color_scheme=color_scheme, figsize=[fig_x, fix_y])
    logo.style_spines(visible=False)
    logo.ax.set_xticks(range(exp_len))
    logo.ax.set_xticklabels('%+d'%x for x in x_idx)

    # Draw boxes for certain regions.
    # logo.highlight_position_range(pmin=0, pmax=2, color='lightcyan')
    # logo.highlight_position_range(pmin=3, pmax=4, color='honeydew')
    # logo.highlight_position_range(pmin=4, pmax=7, color='lavenderblush')

    logo.ax.set_ylabel("probability", labelpad=5, fontsize=10)
    plt.savefig(plot_out, dpi=125)
    plt.close()


################################################################################

def search_generate_html_motif_plots(search_rbps_dic, 
                                     seq_motif_blocks_dic, str_motif_blocks_dic,
                                     out_folder, benchlib_path, motif2db_dic,
                                     html_report_out="motif_plots.rbpbench_search.html",
                                     plot_abs_paths=False,
                                     plots_subfolder="html_motif_plots"):
    """
    Create motif plots for selected RBPs.

    """
    # Use absolute paths?
    if plot_abs_paths:
        out_folder = os.path.abspath(out_folder)
    
    plots_folder = plots_subfolder
    plots_out_folder = out_folder + "/" + plots_folder
    if plot_abs_paths:
        plots_folder = plots_out_folder

    if not os.path.exists(plots_out_folder):
        os.makedirs(plots_out_folder)
    html_out = out_folder + "/" + "motif_plots.rbpbench_search.html"
    md_out = out_folder + "/" + "motif_plots.rbpbench_search.md"
    if html_report_out:
        html_out = html_report_out

    sorttable_js_path = benchlib_path + "/content/sorttable.js"

    # Create theme-specific HTML header.
    mdtext = """
<head>
<title>RBPBench - Motif Plots and Hit Statistics</title>
<script src="%s" type="text/javascript"></script>
</head>
""" %(sorttable_js_path)

    # Add first section markdown.
    mdtext += """

# Motif Plots and Hit Statistics

List of available motif hit statistics and motif plots generated
by RBPBench (rbpbench search --plot-motifs):

- [Motif hit statistics](#motif-hit-stats)
"""

    motif_plot_ids_dic = {}
    idx = 0
    for rbp_id, rbp in sorted(search_rbps_dic.items()):
        idx += 1
        tab_id = "plot-%i" %(idx)
        motif_plot_ids_dic[rbp_id] = tab_id
        mdtext += "- [%s motifs](#%s)\n" %(rbp_id, tab_id)
    mdtext += "\n&nbsp;\n"


    """
    Motif hit statistics table.

    """
    mdtext += """
## Motif hit statistics ### {#motif-hit-stats}

**Table:** RBP motif hit statistics with RBP ID, motif ID, motif database ID (set to "user" if user-supplied motif, otherwise internal database ID), and respective number of motif hits found in supplied genomic regions.

"""

    mdtext += '| &nbsp; RBP ID &nbsp; | &nbsp; Motif ID &nbsp; | Motif database | # motif hits |' + " \n"
    mdtext += "| :-: | :-: | :-: | :-: |\n"
    for rbp_id, rbp in sorted(search_rbps_dic.items()):
        for idx, motif_id in enumerate(rbp.seq_motif_ids):
            c_motif_hits = rbp.seq_motif_hits[idx]
            motif_db = motif2db_dic[motif_id]
            mdtext += "| %s | %s | %s | %i |\n" %(rbp_id, motif_id, motif_db, c_motif_hits)
        for idx, motif_id in enumerate(rbp.str_motif_ids):
            c_motif_hits = rbp.str_motif_hits[idx]
            motif_db = motif2db_dic[motif_id]
            mdtext += "| %s | %s | %s | %i |\n" %(rbp_id, motif_id, motif_db, c_motif_hits)
    mdtext += "\n&nbsp;\n&nbsp;\n"
    mdtext += "\nColumn IDs have the following meanings: "
    mdtext += "**RBP ID** -> RBP ID from database or user-defined (typically RBP name), "
    mdtext += "**Motif ID** -> Motif ID from database or user-defined, "
    mdtext += "**Motif database** -> Motif database used for search run, "
    mdtext += '**# motif hits** -> number of individual motif hits (i.e., hits for motif with motif ID).' + "\n"

    """
    Motif plots.

    """

    for rbp_id, rbp in sorted(search_rbps_dic.items()):
        tab_id = motif_plot_ids_dic[rbp_id]

        # RBP has sequence motifs?
        if rbp.seq_motif_ids:
            mdtext += """
## %s motifs ### {#%s}

RBP "%s" sequence motif plots.

""" %(rbp_id, tab_id, rbp_id)

        else:
            mdtext += """
## %s motifs ### {#%s}

RBP "%s" only contains structure motifs, which are currently not available for plotting.

""" %(rbp_id, tab_id, rbp_id)

        for idx, motif_id in enumerate(rbp.seq_motif_ids):
            c_motif_hits = rbp.seq_motif_hits[idx]
            motif_db = motif2db_dic[motif_id]
            motif_plot = "%s.%s.png" %(rbp_id, motif_id)
            motif_plot_out = plots_out_folder + "/" + motif_plot
            plot_path = plots_folder + "/" + motif_plot

            create_motif_plot(motif_id, seq_motif_blocks_dic,
                              motif_plot_out)

            mdtext += '<img src="' + plot_path + '" alt="' + "sequence motif plot %s" %(motif_id) + "\n"
            mdtext += 'title="' + "sequence motif plot %s" %(motif_id) + '" width="500" />' + "\n"
            mdtext += """

**Figure:** Sequence motif plot for motif ID "%s" (RBP ID: %s, motif database ID: %s). X-axis: motif position. Y-axis: nucleotide probability. Number of %s motif hits in supplied genomic regions: %i.

&nbsp;

""" %(motif_id, rbp_id, motif_db, motif_id, c_motif_hits)


        for idx, motif_id in enumerate(rbp.str_motif_ids):
            c_motif_hits = rbp.str_motif_hits[idx]
            # NO STRUCTURE MOTIF PLOTTING IMPLEMENTED YET.
            # TO DO ...

    # print("Generate motif plots HTML ... ")

    # Convert mdtext to html.
    md2html = markdown(mdtext, extensions=['attr_list', 'tables'])

    OUTMD = open(md_out,"w")
    OUTMD.write("%s\n" %(mdtext))
    OUTMD.close()
    OUTHTML = open(html_out,"w")
    OUTHTML.write("%s\n" %(md2html))
    OUTHTML.close()

    # change <table> to sortable.
    check_cmd = "sed -i 's/<table>/<table class=" + '"sortable"' + ">/g' " + html_out
    output = subprocess.getoutput(check_cmd)
    error = False
    if output:
        error = True
    assert error == False, "sed command returned error:\n%s" %(output)

"""
                 internal_id: str,
                 seq_motif_ids = None,
                 str_motif_ids = None,
                 seq_motif_hits = None,
                 str_motif_hits = None,
                 c_hit_reg = 0, # # regions with motif hits.
                 perc_hit_reg = 0.0, # % hit regions over all regions (i.e. how many input regions contain >= 1 RBP motif).
                 c_motif_hits = 0, # # motif hits.
                 c_uniq_motif_hits = 0, # # unique motif hits.
                 c_uniq_motif_nts = 0, # # unique motif nucleotides.
                 perc_uniq_motif_nts_eff_reg = 0.0, # % unique motif nts over effective region length.
                 perc_uniq_motif_nts_cal_reg = 0.0, # % unique motif nts over called region length.
                 uniq_motif_hits_eff_1000nt = 0.0, # unique motif hits per effective 1000 nt.
                 uniq_motif_hits_cal_1000nt = 0.0, # unique motif hits per called 1000 nt.
                 # ks_pval = 1.0, # Kolmogorov-Smirnov (KS) statistic p-value (are higher scoring sites enriched with motifs).
                 wc_pval = 1.0, # Wilcoxon rank-sum test (Mann-Whitney U test) statistic p-value (are higher scoring sites enriched with motifs).
                 wc_pval_less = 1.0, # Wilcoxon rank-sum test (Mann-Whitney U test) statistic p-value (are lower scoring sites enriched with motifs).
                 organism: Optional[str] = None) -> None:
"""


################################################################################

def create_motif_plot(motif_id, 
                      seq_motif_blocks_dic,
                      motif_plot_out):
    """
    Create sequence motif plot from MEME XML motif block.

    Block format:
    block = ['letter-probability matrix: alength= 4 w= 9 nsites= 20 E= 0', 
    ' 0.000000  0.000000  0.000000  1.000000 ', ' 0.000000  0.000000  0.500000  0.500000 ', 
    ' 0.000000  0.000000  0.000000  1.000000 ', ' 1.000000  0.000000  0.000000  0.000000 ', 
    ' 1.000000  0.000000  0.000000  0.000000 ', ' 0.000000  0.000000  0.000000  1.000000 ', 
    ' 0.500000  0.000000  0.500000  0.000000 ', ' 0.000000  0.000000  0.000000  1.000000 ', 
    ' 0.000000  0.000000  0.000000  1.000000 ']

    """

    assert motif_id in seq_motif_blocks_dic, "motif ID %s not found in seq_motif_blocks_dic" %(motif_id)
    motif_rows = seq_motif_blocks_dic[motif_id]

    motif_freqs_dic = {}
    motif_freqs_dic["A"] = []
    motif_freqs_dic["C"] = []
    motif_freqs_dic["G"] = []
    motif_freqs_dic["T"] = []

    for idx, row in enumerate(motif_rows):
        if idx == 0:
            continue
        freqs_pre = row.split(" ")
        freqs = []
        for f in freqs_pre:
            if f:
                freqs.append(float(f))
        
        assert len(freqs) == 4, "invalid number of nucleotide frequencies (expected 4 but found %i)" %(len(freqs))
        motif_freqs_dic["A"].append(freqs[0])
        motif_freqs_dic["C"].append(freqs[1])
        motif_freqs_dic["G"].append(freqs[2])
        motif_freqs_dic["T"].append(freqs[3])

    df = pd.DataFrame(motif_freqs_dic)
    # df = pd.DataFrame(motif_freqs_dic, columns=rbp_list)
    logo = Logo(df)
    logo.style_spines(visible=False)
    #logo.style_spines(spines=['left'], visible=False, bounds=[0, 1])
    #logo.style_spines(spines=['bottom'], visible=False)
    #logo.ax.set_xticks([])
    #logo.ax.set_yticks([])
    #plt.yticks(fontsize=7)
    logo.ax.set_ylabel("probability", labelpad=10, fontsize=10)
    plt.savefig(motif_plot_out, dpi=100)
    plt.close()


################################################################################

def compare_generate_html_report(compare_methods_dic, compare_datasets_dic,
                                 rbp_stats_dic, motif_stats_dic,
                                 out_folder, benchlib_path,
                                 html_report_out="report.rbpbench_compare.html",
                                 plot_abs_paths=False,
                                 plots_subfolder="html_plots"):
    """
    Create comparison statistics and HTML report.

    """
    # Use absolute paths?
    if plot_abs_paths:
        out_folder = os.path.abspath(out_folder)

    plots_folder = plots_subfolder
    plots_out_folder = out_folder + "/" + plots_folder
    if plot_abs_paths:
        plots_folder = plots_out_folder

    if not os.path.exists(plots_out_folder):
        os.makedirs(plots_out_folder)
    html_out = out_folder + "/" + "report.rbpbench_compare.html"
    md_out = out_folder + "/" + "report.rbpbench_compare.md"
    if html_report_out:
        html_out = html_report_out

    # # Logo paths.
    # logo1_path = benchlib_path + "/content/logo1.png"
    # logo2_path = benchlib_path + "/content/logo2.png"
    sorttable_js_path = benchlib_path + "/content/sorttable.js"
    # # plotly js path.
    # plotly_js_path = benchlib_path + "/content/plotly-latest.min.js"
    # assert os.path.exists(plotly_js_path), "plotly js %s not found" %(plotly_js_path)


    # Create theme-specific HTML header.
    mdtext = """
<head>
<title>RBPBench - Motif Search Comparison Report</title>
<script src="%s" type="text/javascript"></script>
</head>
""" %(sorttable_js_path)

    # Add first section markdown.
    mdtext += """

# Comparison Report

List of available comparison statistics generated
by RBPBench (rbpbench compare):

"""
    # Comparison based on method ID: tables.
    method_table_ids_dic = {}
    idx = 0
    for comp_id, data in sorted(compare_methods_dic.items()):
        # mdtext += "\n"
        if len(data) < 2:
            continue
        idx += 1
        tab_id = "method-tab-%i" %(idx)
        method_table_ids_dic[comp_id] = tab_id
        mdtext += "- [%s method comparison table](#%s)\n" %(comp_id, tab_id)

    # Comparison based on data ID: tables.
    data_table_ids_dic = {}
    idx = 0
    for comp_id, data in sorted(compare_datasets_dic.items()):
        # mdtext += "\n"
        if len(data) < 2:
            continue
        idx += 1
        tab_id = "data-tab-%i" %(idx)
        data_table_ids_dic[comp_id] = tab_id
        mdtext += "- [%s dataset comparison table](#%s)\n" %(comp_id, tab_id)

    # Comparison based on method ID: Venn diagrams.
    method_plot_ids_dic = {}
    idx = 0
    for comp_id, data in sorted(compare_methods_dic.items()):
        if len(data) < 2:
            continue
        idx += 1
        tab_id = "method-venn-%i" %(idx)
        method_plot_ids_dic[comp_id] = tab_id
        mdtext += "- [%s method comparison plot](#%s)\n" %(comp_id, tab_id)

    # Comparison based on data ID: Venn diagrams.
    data_plot_ids_dic = {}
    idx = 0
    for comp_id, data in sorted(compare_datasets_dic.items()):
        # mdtext += "\n"
        if len(data) < 2:
            continue
        idx += 1
        tab_id = "data-venn-%i" %(idx)
        data_plot_ids_dic[comp_id] = tab_id
        mdtext += "- [%s dataset comparison plot](#%s)\n" %(comp_id, tab_id)
    mdtext += "\n&nbsp;\n"


    # Method comparison tables.
    for comp_id, data in sorted(compare_methods_dic.items()):
        if len(data) < 2:
            continue
        tab_id = method_table_ids_dic[comp_id]
        mdtext += """
## %s method comparison statistics ### {#%s}

**Table:** RBP motif hit statistics for combined ID "%s" (data ID, motif database ID, RBP ID) over different methods (method ID column).

""" %(comp_id, tab_id, comp_id)
        mdtext += '| Method ID | # regions | # motif hits | % regions with motifs | % motif nucleotides | # motif hits per 1000 nt |' + " \n"
        mdtext += "| :-: | :-: | :-: | :-: | :-: | :-: |\n"
        for method_id in data:
            int_id = data[method_id]
            c_regions = rbp_stats_dic[int_id].c_regions
            c_uniq_motif_hits = rbp_stats_dic[int_id].c_uniq_motif_hits
            perc_reg_with_hits = rbp_stats_dic[int_id].perc_reg_with_hits
            perc_uniq_motif_nts_eff_reg = rbp_stats_dic[int_id].perc_uniq_motif_nts_eff_reg
            uniq_motif_hits_cal_1000nt = rbp_stats_dic[int_id].uniq_motif_hits_cal_1000nt
            mdtext += "| %s | %i | %i | %.2f | %.2f | %.2f |\n" %(method_id, c_regions, c_uniq_motif_hits, perc_reg_with_hits, perc_uniq_motif_nts_eff_reg, uniq_motif_hits_cal_1000nt)
        mdtext += "\n&nbsp;\n&nbsp;\n"
        mdtext += "\nColumn IDs have the following meanings: "
        mdtext += "**Method ID** -> method ID set for dataset (typically peak calling method ID), "
        mdtext += '**# regions** -> number of peak regions used for motif search, '
        mdtext += '**# motif hits** -> number of unique motif hits in peak regions (removed double counts), '
        mdtext += '**% regions with motifs** -> percentage of peak regions with motif hits, '
        mdtext += '**% motif nucleotides** -> percentage of unique motif nucleotides over effective peak region size (overlapping regions merged), '
        mdtext += '**# motif hits per 1000 nt** -> number of motif hits over 1000 nt of called peak region size (overlapping regions NOT merged).' + "\n"

    # Dataset comparison tables.
    for comp_id, data in sorted(compare_datasets_dic.items()):
        if len(data) < 2:
            continue
        tab_id = data_table_ids_dic[comp_id]
        mdtext += """
## %s dataset comparison statistics ### {#%s}

**Table:** RBP motif hit statistics for combined ID "%s" (method ID, motif database ID, RBP ID) over different datasets (data ID column).

""" %(comp_id, tab_id, comp_id)
        mdtext += '| Data ID | # regions | # motif hits | % regions with motifs | % motif nucleotides | # motif hits per 1000 nt |' + " \n"
        mdtext += "| :-: | :-: | :-: | :-: | :-: | :-: |\n"
        for dataset_id in data:
            int_id = data[dataset_id]
            c_regions = rbp_stats_dic[int_id].c_regions
            c_uniq_motif_hits = rbp_stats_dic[int_id].c_uniq_motif_hits
            perc_reg_with_hits = rbp_stats_dic[int_id].perc_reg_with_hits
            perc_uniq_motif_nts_eff_reg = rbp_stats_dic[int_id].perc_uniq_motif_nts_eff_reg
            uniq_motif_hits_cal_1000nt = rbp_stats_dic[int_id].uniq_motif_hits_cal_1000nt
            mdtext += "| %s | %i | %i | %.2f | %.2f | %.2f |\n" %(dataset_id, c_regions, c_uniq_motif_hits, perc_reg_with_hits, perc_uniq_motif_nts_eff_reg, uniq_motif_hits_cal_1000nt)
        mdtext += "\n&nbsp;\n&nbsp;\n"
        mdtext += "\nColumn IDs have the following meanings: "
        mdtext += "**Data ID** -> data ID set for dataset (typically describing CLIP data, e.g. CLIP method + cell type combination), "
        mdtext += '**# regions** -> number of peak regions used for motif search, '
        mdtext += '**# motif hits** -> number of unique motif hits in peak regions (removed double counts), '
        mdtext += '**% regions with motifs** -> percentage of peak regions with motif hits, '
        mdtext += '**% motif nucleotides** -> percentage of unique motif nucleotides over effective peak region size (overlapping regions merged), '
        mdtext += '**# motif hits per 1000 nt** -> number of motif hits over 1000 nt of called peak region size (overlapping regions NOT merged).' + "\n"

    """
    Venn diagrams for method ID comparisons.

    """
    for comp_id, method_dic in sorted(compare_methods_dic.items()):
        if len(method_dic) < 2:
            continue
        tab_id = method_plot_ids_dic[comp_id]

        int_ids = []
        method_ids = []
        for method_id, int_id in sorted(method_dic.items()):
            int_ids.append(int_id)
            method_ids.append(method_id)

        venn_plot = "venn_diagram.method_comp.%s.png" %(comp_id)
        venn_plot_out = plots_out_folder + "/" + venn_plot
        plot_path = plots_folder + "/" + venn_plot

        c_methods = len(method_ids)
        method_ids_str = ",".join(method_ids)

        if len(method_ids) == 2:
            create_venn2_diagram(int_ids[0], int_ids[1],
                            motif_stats_dic,
                            venn_plot_out,
                            set1_label=method_ids[0],
                            set2_label=method_ids[1])
        elif len(method_ids) == 3:
            create_venn3_diagram(int_ids[0], int_ids[1], int_ids[2],
                            motif_stats_dic,
                            venn_plot_out,
                            set1_label=method_ids[0],
                            set2_label=method_ids[1],
                            set3_label=method_ids[2])
        elif len(method_ids) > 3 and len(method_ids) <= 24:
            create_vennx_diagram(int_ids, method_ids,
                            motif_stats_dic, venn_plot_out)
        else:
            assert False, "two many methods to compare (comp_id: %s). Please use less methods for plotting (current limit: 24)" %(comp_id)


        mdtext += """
## %s method comparison plot ### {#%s}

Based on the same combined ID "%s" (data ID, motif database ID, RBP ID), motif hit occurrences for %i different methods (%s) are compared via Venn diagram.
Any given motif hit can either be found only by one method, or be identified by any set (>=2) of methods (intersection areas).

""" %(comp_id, tab_id, comp_id, c_methods, method_ids_str)
        mdtext += '<img src="' + plot_path + '" alt="' + "dataset comparison plot %s" %(comp_id) + "\n"
        mdtext += 'title="' + "dataset comparison plot %s" %(comp_id) + '" width="700" />' + "\n"
        mdtext += """

**Figure:** Venn diagram of motif hit occurrences for %i different methods (%s) with identical combined ID (%s) + corresponding percentages of total motif hits for each region (method exclusive and intersection(s)).

&nbsp;

""" %(c_methods, method_ids_str, comp_id)

    """
    Venn diagrams for data ID comparisons.

    """
    for comp_id, data_dic in sorted(compare_datasets_dic.items()):
        if len(data_dic) < 2:
            continue
        tab_id = data_plot_ids_dic[comp_id]

        int_ids = []
        data_ids = []
        for data_id, int_id in sorted(data_dic.items()):
            int_ids.append(int_id)
            data_ids.append(data_id)

        venn_plot = "venn_diagram.data_comp.%s.png" %(comp_id)
        venn_plot_out = plots_out_folder + "/" + venn_plot
        plot_path = plots_folder + "/" + venn_plot

        c_data_ids = len(data_ids)
        data_ids_str = ",".join(data_ids)

        if len(data_ids) == 2:
            create_venn2_diagram(int_ids[0], int_ids[1],
                            motif_stats_dic,
                            venn_plot_out,
                            set1_label=data_ids[0],
                            set2_label=data_ids[1])
        elif len(data_ids) == 3:
            create_venn3_diagram(int_ids[0], int_ids[1], int_ids[2],
                            motif_stats_dic,
                            venn_plot_out,
                            set1_label=data_ids[0],
                            set2_label=data_ids[1],
                            set3_label=data_ids[2])
        elif len(data_ids) > 3 and len(data_ids) <= 24:
            create_vennx_diagram(int_ids, data_ids,
                            motif_stats_dic, venn_plot_out)
        else:
            assert False, "two many datasets to compare (comp_id: %s). Please use less datasets for plotting (current limit: 24)" %(comp_id)

        method_id = comp_id.split(",")[0]

        mdtext += """
## %s dataset comparison plot ### {#%s}

Based on the same method combination ID "%s" (method ID, motif database ID, RBP ID), motif hit occurrences in %i different datasets (%s) are compared via Venn diagram.
Any given motif hit can either be found only in one dataset, or be common to any >= 2 datasets (intersection areas).

""" %(comp_id, tab_id, comp_id, c_data_ids, data_ids_str)
        mdtext += '<img src="' + plot_path + '" alt="' + "dataset comparison plot %s" %(comp_id) + "\n"
        mdtext += 'title="' + "dataset comparison plot %s" %(comp_id) + '" width="700" />' + "\n"
        mdtext += """

**Figure:** Venn diagram of motif hit occurrences for %i different datasets (%s) with identical combined ID (%s) + corresponding percentages of total motif hits for each region (method exclusive and intersection(s)).

&nbsp;

""" %(c_data_ids, data_ids_str, comp_id)

    """
    Method comparisons.

    Two comparisons here:
    compare_methods_dic: {
    'k562_eclip,human_v0.1,PUM1': {'clipper_idr': ['GeR1WaFK'], 'dewseq_w100_s5': ['aErgL-gv']}, 
    'k562_eclip,human_v0.1,PUM2': {'clipper_idr': ['fEag63xM'], 'dewseq_w100_s5': ['ZzLuEZR9']}
    }
    rbp_stats_dic
    motif_stats_dic
    
    """

    print("Generate HTML report ... ")

    # Convert mdtext to html.
    md2html = markdown(mdtext, extensions=['attr_list', 'tables'])

    OUTMD = open(md_out,"w")
    OUTMD.write("%s\n" %(mdtext))
    OUTMD.close()
    OUTHTML = open(html_out,"w")
    OUTHTML.write("%s\n" %(md2html))
    OUTHTML.close()

    # change <table> to sortable.
    check_cmd = "sed -i 's/<table>/<table class=" + '"sortable"' + ">/g' " + html_out
    output = subprocess.getoutput(check_cmd)
    error = False
    if output:
        error = True
    assert error == False, "sed command returned error:\n%s" %(output)


################################################################################

def create_venn2_diagram(int_id1, int_id2,
                         motif_stats_dic,
                         out_plot,
                         alpha=0.5,
                         set1_label="Set1",
                         set2_label="Set2"):
    """
    Create Venn Diagram for two sets.

    """
    assert int_id1 in motif_stats_dic, "given internal_id %s not in motif_stats_dic" %(int_id1)
    assert int_id2 in motif_stats_dic, "given internal_id %s not in motif_stats_dic" %(int_id2)

    # Generate sets.
    set1_list = []
    set2_list = []
    for stats in motif_stats_dic[int_id1]:
        set1_list.append(stats.hit_id)
    for stats in motif_stats_dic[int_id2]:
        set2_list.append(stats.hit_id)

    set1 = set(set1_list)
    set2 = set(set2_list)

    # # Intersection.
    # intersection_set = set1 & set2
    # # Difference.
    # only_set1 = set1.difference(set2)
    # only_set2 = set2.difference(set1)

    total = len(set1.union(set2))
    venn2([set1, set2], #set_colors=('skyblue', 'salmon'), 
          alpha=alpha, set_labels = (set1_label, set2_label),
          subset_label_formatter=lambda x: str(x) + "\n(" + f"{(x/total):1.0%}" + ")")
    plt.savefig(out_plot, dpi=150)
    plt.clf()


################################################################################

def create_venn3_diagram(int_id1, int_id2, int_id3,
                         motif_stats_dic,
                         out_plot,
                         alpha=0.5,
                         set1_label="Set1",
                         set2_label="Set2",
                         set3_label="Set3"):
    """
    Create Venn Diagram for three sets.

    #36e6e6
    #ee2a9a
    #f8e318

    """
    assert int_id1 in motif_stats_dic, "given internal_id %s not in motif_stats_dic" %(int_id1)
    assert int_id2 in motif_stats_dic, "given internal_id %s not in motif_stats_dic" %(int_id2)
    assert int_id3 in motif_stats_dic, "given internal_id %s not in motif_stats_dic" %(int_id3)

    # Generate sets.
    set1_list = []
    set2_list = []
    set3_list = []
    for stats in motif_stats_dic[int_id1]:
        set1_list.append(stats.hit_id)
    for stats in motif_stats_dic[int_id2]:
        set2_list.append(stats.hit_id)
    for stats in motif_stats_dic[int_id3]:
        set3_list.append(stats.hit_id)

    """
    Set union operator: |
    Set intersection operator: &
    """

    set1 = set(set1_list)
    set2 = set(set2_list)
    set3 = set(set3_list)

    total = len(set1 | set2 | set3) 

    venn3([set1, set2, set3], # set_colors=('#36e6e6', '#ee2a9a', '#f8e318'), 
          alpha=alpha, set_labels = (set1_label, set2_label, set3_label),
          subset_label_formatter=lambda x: str(x) + "\n(" + f"{(x/total):1.0%}" + ")")
    plt.savefig(out_plot, dpi=150)
    plt.clf()


################################################################################

def create_vennx_diagram(int_ids, set_labels,
                         motif_stats_dic, out_plot,
                         cmap="jet"):
    """
    Create Venn Diagram for x sets.

    Some matplotlib color map styles to try:
    jet, nipy_spectral, tab20, brg, turbo, gnuplot
    https://matplotlib.org/stable/tutorials/colors/colormaps.html


    """

    for int_id in int_ids:
        assert int_id in motif_stats_dic, "given internal_id %s not in motif_stats_dic" %(int_id)
    assert len(int_ids) == len(set_labels), "len(int_ids) != len(set_labels)"

    set_dic = {}
    for idx, int_id in enumerate(int_ids):
        set_label = set_labels[idx]
        id_list = []
        for stats in motif_stats_dic[int_id]:
            id_list.append(stats.hit_id)
        id_set = set(id_list)
        set_dic[set_label] = id_set

    venn(set_dic,
         cmap=cmap,
         # alpha=0.6,
         fmt="{size}"+"\n"+"{percentage:.1f}%")
    plt.savefig(out_plot, dpi=100)
    plt.clf()


################################################################################

def print_some_banner():
    """
    Print some banner.

    """

    a = """
              ######  ######  ######             
              #     # #     # #     #            
              #     # #     # #     #            
              ######  ######  ######             
              #   #   #     # #                  
              #    #  #     # #
   ##         #     # ######  #              ##
 ####                                        ####
##################################################
 ####                                        ####
   ##   #####  ###### #    #  ####  #    #   ##
        #    # #      ##   # #    # #    #            
        #####  #####  # #  # #      ###### 
        #    # #      #  # # #      #    # 
        #    # #      #   ## #    # #    # 
        #####  ###### #    #  ####  #    #

"""

    return(a)

################################################################################
