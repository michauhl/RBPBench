from distutils.spawn import find_executable
from typing import Optional
import os
import re
import subprocess
import gzip
import shutil
import matplotlib.pyplot as plt
from matplotlib_venn import venn2, venn3
from venn import venn
from logomaker import Logo
from markdown import markdown
import pandas as pd
import plotly.express as px
from math import log10
import textdistance
import numpy as np
# from upsetplot import plot
from upsetplot import UpSet
# from matplotlib import cm
from packaging import version


"""

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
~~~~~~ OPEN FOR BUSINESS ~~~~~~
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


~~~~~~~~~~~~~
Run doctests
~~~~~~~~~~~~~

python3 -m doctest benchlib.py
python3 -m doctest -v benchlib.py



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

def check_tool_version(tool_call, min_version):
    """
    Check if version of a command line tool is >= min_version.
    Use tool_call to provide command line call which produces version number.
    E.g. check_tool_version("meme -version", "5.3.0")

    Using:
    from packaging import version

    """

    output = subprocess.getoutput(tool_call)
    # version_string = output.split('\n')[0]
    version_string = output.strip()
    tool_version = version.parse(version_string)
    min_version = version.parse(min_version)
    check = tool_version >= min_version
    return check, tool_version


################################################################################

def get_colormap_hex_colors(cm_id):
    """
    Get all colors of a matplotlib colormap (provided ID cm_id), as HEX values.
    E.g. "Pastel1", "Spectral", ...
    Return list of HEX color values.
    
    Also see:
    https://matplotlib.org/stable/gallery/color/colormap_reference.html

    """
    import matplotlib

    cmap = plt.cm.get_cmap(cm_id)
    colors = cmap(np.linspace(0,1,cmap.N))

    hex_colors = []
    for color in colors:
        hex_color = matplotlib.colors.rgb2hex(color)
        hex_colors.append(hex_color)

    return(hex_colors)


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

def make_contingency_table_max_dist_2x2(region_rbp_motif_pos_dic, rbp_id1, rbp_id2,
                                        motif_max_dist=50,
                                        count_to_d=True):

    """
    Make a contingency table 2x2, but only count regions as occupied by both 
    RBPs if any of their motifs are <= motif_max_dist away from each other.
    Count regions with > motif_max_dist as c_out, but also add them to D (otherwise
    we get different total counts for each RBP pair).
    
    """

    cont_a = 0
    cont_b = 0
    cont_c = 0
    cont_d = 0
    c_out = 0

    for reg_id in region_rbp_motif_pos_dic:
        if rbp_id1 in region_rbp_motif_pos_dic[reg_id]:
            if rbp_id2 in region_rbp_motif_pos_dic[reg_id]:
                # Compare motif distances between rbp_id1 + rbp_id2.
                found_pair = False
                for p1 in region_rbp_motif_pos_dic[reg_id][rbp_id1]:
                    for p2 in region_rbp_motif_pos_dic[reg_id][rbp_id2]:
                        if abs(p1 - p2) <= motif_max_dist:
                            found_pair = True
                            break
                if found_pair:
                    cont_a += 1
                else:
                    c_out += 1
                    if count_to_d:
                        cont_d += 1
            else:
                cont_c += 1
        else:
            if rbp_id2 in region_rbp_motif_pos_dic[reg_id]:
                cont_b += 1
            else:
                cont_d += 1

    table = [[cont_a, cont_b], [cont_c, cont_d]]
    return table, c_out


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

def read_in_cm_blocks(cm_file, 
                      empty_check=True):
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
                # Remove special characters from motif/accession ID.
                new_acc_id = remove_special_chars_from_str(acc_id)
                assert new_acc_id, "no characters left after removal of special characters from covariance model accession ID \"%s\". Please use valid covariance model accession IDs (i.e., modify ACC column strings in .cm file)" %(acc_id)
                acc_list.append(new_acc_id)
                blocks_list[idx] += line
            else:
                blocks_list[idx] += line
    f.closed

    if empty_check:
        assert blocks_list, "no covariance models read in from %s (blocks_list empty). Please provide valid covariance model(s) (.cm) file" %(cm_file)
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
                  empty_check=True):
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

    if empty_check:
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

def get_fasta_headers(in_fa,
                      full_header=False):
    """
    Get FASTA header IDs.
    This grep version appears to be much faster than reading in file via 
    Python line by line.

    full_header:
        If true, use whole header (after >) as ID. By default, use ID up to 
        first space character.

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
            if full_header:
                m = re.search(">(.+)", line)
                seq_id = m.group(1)
                seq_ids_dic[seq_id] = 1
            else:
                m = re.search(">(\S+)", line)
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
                        full_header=False,
                        report=1,
                        all_uc=False,
                        empty_check=True,
                        id_check=True,
                        skip_data_id="set",
                        skip_n_seqs=True):
    """
    Read in FASTA sequences, store in dictionary and return dictionary.
    FASTA file can be plain text or gzipped (watch out for .gz ending).

    full_header:
        If true, use whole header (after >) as ID. By default, use ID up to 
        first space character.

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
            m = False
            if full_header:
                m = re.search(">(.+)", line)
            else:
                m = re.search(">(\S+)", line)

            assert m, "header ID extraction failed for FASTA header line \"%s\"" %(line)
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
                       empty_check=True):
    """
    Read in XML motifs, store in blocks dictionary.

    """

    assert os.path.exists(meme_xml_file), "meme_xml_file %s not found" % (meme_xml_file)

    raw_text = ""
    with open(meme_xml_file, "r") as f:
        raw_text = f.read()
    assert raw_text, "nothing read in from MEME XML file %s. Please provide valid MEME/DREME XML sequence motifs file" %(meme_xml_file)

    # Get motif blocks.
    motif_blocks_dic = extract_motif_blocks(raw_text)

    if empty_check:
        assert motif_blocks_dic, "motif_blocks_dic empty (malformatted MEME/DREME XML file provided?)"

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
            new_motif_id = remove_special_chars_from_str(motif_id)
            assert new_motif_id, "no characters left after removal of special characters from motif ID \"%s\". Please use valid MEME XML motif IDs (i.e., modify MOTIF column strings in motifs xml file)" %(motif_id)
            motif_id = new_motif_id
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

def check_convert_chr_id(chr_id,
                         id_style=1):
    """
    Check and convert chromosome IDs to format:
    chr1, chr2, chrX, ...
    If chromosome IDs like 1,2,X, .. given, convert to chr1, chr2, chrX ..
    Return False if given chr_id not standard and not convertable.

    Filter out scaffold IDs like:
    GL000009.2, KI270442.1, chr14_GL000009v2_random
    chrUn_KI270442v1 ...

    id_style:
        Defines to which style chromosome IDs should be converted to.
        0: Do not convert at all, just return chr_id.
        1: to chr1,chr2,...,chrM style.
        2: to 1,2,...,MT style.

    >>> chr_id = "chrX"
    >>> check_convert_chr_id(chr_id)
    'chrX'
    >>> chr_id = "4"
    >>> check_convert_chr_id(chr_id)
    'chr4'
    >>> chr_id = "MT"
    >>> check_convert_chr_id(chr_id)
    'chrM'
    >>> chr_id = "GL000009.2"
    >>> check_convert_chr_id(chr_id)
    False
    >>> chr_id = "chrUn_KI270442v1"
    >>> check_convert_chr_id(chr_id)
    False
    >>> chr_id = "chr2R"
    >>> check_convert_chr_id(chr_id)
    'chr2R'
    >>> chr_id = "3L"
    >>> check_convert_chr_id(chr_id)
    'chr3L'
    >>> chr_id = "4L"
    >>> check_convert_chr_id(chr_id)
    False
    >>> chr_id = "chrM"
    >>> check_convert_chr_id(chr_id, id_style=2)
    'MT'
    >>> chr_id = "chr2R"
    >>> check_convert_chr_id(chr_id, id_style=2)
    '2R'
    >>> chr_id = "5"
    >>> check_convert_chr_id(chr_id, id_style=2)
    '5'
    >>> chr_id = "chrA"
    >>> check_convert_chr_id(chr_id, id_style=2)
    False
    >>> chr_id = "chrA"
    >>> check_convert_chr_id(chr_id, id_style=0)
    'chrA'


    """
    assert chr_id, "given chr_id empty"

    if not id_style: # If id_style == 0 or False.
        return chr_id

    elif id_style == 1:
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

    elif id_style == 2:

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
        assert False, "invalid id_style set"


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

def gtf_read_in_gene_infos(in_gtf,
                           tr2gid_dic=None,
                           tr_types_dic=None,
                           check_chr_ids_dic=None,
                           chr_style=0,
                           empty_check=False):
    """
    Read in gene infos into GeneInfo objects, including information on 
    transcript isoforms for the gene. Note that only features on standard 
    chromosomes (1,2,...,X Y MT) are currently used.

    Assuming gtf file with order: gene,transcript(s),exon(s) ...

    tr_types_dic:
        Store transcript biotype IDs and number of appearances.
        transcript biotype ID -> # appearances
    
    chr_style:
        0: do not change
        1: change to chr1, chr2 ...
        2: change to 1, 2, 3, ...

    """

    # Transcript ID to exonic length dictionary.
    tr2len_dic = {}
    # Gene info objects dictionary (gene_id -> gene info object).
    gid2gio_dic = {}

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

        chr_id = check_convert_chr_id(chr_id, id_style=chr_style)
        # If not one of standard chromosomes, continue.
        if not chr_id:
            continue

        assert feat_e >= feat_s, "feature end < feature start in GTF file \"%s\", line \"%s\". Since both coordinates are expected to have 1-based index, this should not happen" %(in_gtf, line)

        if check_chr_ids_dic is not None:
            assert chr_id in check_chr_ids_dic, "chromosme ID \"%s\" from GTF file \"%s\" not found in --genome FASTA file. Make sure GTF + FASTA + BED files use compatible chromosome IDs" %(chr_id, in_gtf)

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
            # assert m, "gene_biotype / gene_type entry missing in GTF file \"%s\", line \"%s\"" %(in_gtf, line)
            # gene_biotype = m.group(1)

            gene_infos = GeneInfo(gene_id, gene_name, gene_biotype, chr_id, feat_s, feat_e, feat_pol)

            assert gene_id not in gid2gio_dic, "gene feature with gene ID %s already encountered in GTF file \"%s\"" %(gene_id, in_gtf)
            gid2gio_dic[gene_id] = gene_infos

        elif feature == "transcript":

            m = re.search('gene_id "(.+?)"', infos)
            assert m, "gene_id entry missing in GTF file \"%s\", line \"%s\"" %(in_gtf, line)
            gene_id = m.group(1)
            m = re.search('transcript_id "(.+?)"', infos)
            assert m, "transcript_id entry missing in GTF file \"%s\", line \"%s\"" %(in_gtf, line)
            tr_id = m.group(1)
            assert gene_id in gid2gio_dic, "gene_id %s belonging to transcript ID %s not (yet) encountered. Gene feature expected to come before transcript and exon features in GTF file \"%s\"" %(gene_id, tr_id, in_gtf)

            if tr2gid_dic is not None:
                tr2gid_dic[tr_id] = gene_id

            tr_biotype = "-"  # # optional.
            m = re.search('transcript_biotype "(.+?)"', infos)
            if not m:
                m = re.search('transcript_type "(.+?)"', infos)
            # assert m, "transcript_biotype / transcript_type entry missing in GTF file \"%s\", line \"%s\"" %(in_gtf, line)
            if m:
                tr_biotype = m.group(1)

            if tr_types_dic is not None:
                if tr_biotype in tr_types_dic:
                    tr_types_dic[tr_biotype] += 1
                else:
                    tr_types_dic[tr_biotype] = 1

            # Basic tag.
            basic_tag = 0
            m = re.search('tag "basic"', infos)
            if m:
                basic_tag = 1
            # Ensembl canonical.
            ensembl_canonical = 0
            m = re.search('tag "Ensembl_canonical"', infos)
            if m:
                ensembl_canonical = 1
            # Transcript support level (TSL).
            m = re.search('transcript_support_level "(.+?)"', infos)
            tsl_id = "NA"
            if m:
                tsl_id = m.group(1)
                if re.search("assigned to previous", tsl_id):
                    m = re.search("(.+?) \(", tsl_id)
                    tsl_id = m.group(1)
            # Dummy length for now.
            tr_length = 0

            # Update gene infos.
            gid2gio_dic[gene_id].tr_ids.append(tr_id)
            gid2gio_dic[gene_id].tr_biotypes.append(tr_biotype)
            gid2gio_dic[gene_id].tr_basic_tags.append(basic_tag)
            gid2gio_dic[gene_id].tr_ensembl_canonical_tags.append(ensembl_canonical)
            gid2gio_dic[gene_id].tr_tsls.append(tsl_id)
            gid2gio_dic[gene_id].tr_lengths.append(tr_length)

        elif feature == "exon":
            # Extract transcript ID.
            m = re.search('transcript_id "(.+?)"', infos)
            assert m, "transcript_id entry missing in GTF file \"%s\", line \"%s\"" %(in_gtf, line)    
            tr_id = m.group(1)
            # Sum up length.
            ex_len = feat_e - feat_s + 1
            if not tr_id in tr2len_dic:
                tr2len_dic[tr_id] = ex_len
            else:
                tr2len_dic[tr_id] += ex_len

    f.close()

    if empty_check:
        assert gid2gio_dic, "no gene infos read in from GTF file \"%s\"" %(in_gtf)

    # Update transcript lengths.
    for gene_id in gid2gio_dic:
        for idx, tr_id in enumerate(gid2gio_dic[gene_id].tr_ids):
            gid2gio_dic[gene_id].tr_lengths[idx] = tr2len_dic[tr_id]

    return gid2gio_dic


################################################################################

def gtf_read_in_transcript_infos(in_gtf,
                                 tr_ids_dic=None,
                                 tr_types_dic=None,
                                 correct_min_ex_order=True,
                                 chr_style=0,
                                 empty_check=True):
    """
    Read in transcript infos into TransriptInfo objects. Note that only 
    features on standard chromosomes (1,2,...,X Y MT) are currently used.

    tr_ids_dic:
        Transcript IDs dictionary to define which transcript IDs to read in.
    correct_min_ex_order:
        If True, assume minus strand exon numbering is correct, i.e., most 
        downstream exon gets exon number 1. In some GTF files, most upstream 
        exon get exon number 1, which is incorrect as it is the last exon 
        of the minus strand transcript (in this case set to False).
    chr_style:
        0: do not change
        1: change to chr1, chr2 ...
        2: change to 1, 2, 3, ...

    Assuming gtf file with order: gene,transcript(s),exon(s) ...

    """

    # Transcript info objects dictionary (transcript ID -> transcript info object).
    tid2tio_dic = {}
    tr_ids_seen_dic = {}

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

        chr_id = check_convert_chr_id(chr_id, id_style=chr_style)
        # If not one of standard chromosomes, continue.
        if not chr_id:
            continue

        if feature == "transcript":

            m = re.search('gene_id "(.+?)"', infos)
            assert m, "gene_id entry missing in GTF file \"%s\", line \"%s\"" %(in_gtf, line)
            gene_id = m.group(1)
            m = re.search('transcript_id "(.+?)"', infos)
            assert m, "transcript_id entry missing in GTF file \"%s\", line \"%s\"" %(in_gtf, line)
            tr_id = m.group(1)

            tr_ids_seen_dic[tr_id] = 1

            if tr_ids_dic is not None:
                if tr_id not in tr_ids_dic:
                    continue

            tr_biotype = "-"  # optional.
            m = re.search('transcript_biotype "(.+?)"', infos)
            if not m:
                m = re.search('transcript_type "(.+?)"', infos)
            if m:
                tr_biotype = m.group(1)

            if tr_types_dic is not None:
                if tr_biotype not in tr_types_dic:
                    tr_types_dic[tr_biotype] = 1
                else:
                    tr_types_dic[tr_biotype] += 1

            tr_infos = TranscriptInfo(tr_id, tr_biotype, chr_id, feat_s, feat_e, feat_pol, gene_id,
                                      tr_length=0,
                                      exon_c=0)
            assert tr_id not in tid2tio_dic, "transcript feature with transcript ID %s already encountered in GTF file \"%s\"" %(tr_id, in_gtf)
            tid2tio_dic[tr_id] = tr_infos

        elif feature == "exon":

            m = re.search('transcript_id "(.+?)"', infos)
            assert m, "transcript_id entry missing in GTF file \"%s\", line \"%s\"" %(in_gtf, line)    
            tr_id = m.group(1)
            assert tr_id in tr_ids_seen_dic, "transcript ID %s in exon feature not (yet) encountered. Transcript feature expected to come before exon features in GTF file \"%s\"" %(tr_id, in_gtf)
    
            if tr_ids_dic is not None:
                if tr_id not in tr_ids_dic:
                    continue

            m = re.search('exon_number "(\d+?)"', infos)
            if not m:
                m = re.search('exon_number (\d+?);', infos)  # GENCODE encoding.
            assert m, "exon_number entry missing in GTF file \"%s\", line \"%s\"" %(in_gtf, line)
            exon_nr = int(m.group(1))

            tid2tio_dic[tr_id].exon_c += 1

            assert exon_nr == tid2tio_dic[tr_id].exon_c, "ascending exon numbers expected in GTF file \"%s\" but instead found: %i (exon number) %i (expected number) in line:\n%s" %(in_gtf, exon_nr, tid2tio_dic[tr_id].exon_c, line)

            if correct_min_ex_order:
                # If correct minus trand exon order, add exon coordinates to end of list.
                tid2tio_dic[tr_id].exon_coords.append([feat_s, feat_e])
            else:
                # If reverse minus strand exon order, insert exon coordinates at beginning.
                tid2tio_dic[tr_id].exon_coords.insert(0, [feat_s, feat_e])

            # Record exon length.
            exon_length = feat_e - feat_s + 1
            tid2tio_dic[tr_id].tr_length += exon_length

        elif feature == "CDS":

            m = re.search('transcript_id "(.+?)"', infos)
            assert m, "transcript_id entry missing in GTF file \"%s\", line \"%s\"" %(in_gtf, line)    
            tr_id = m.group(1)
            assert tr_id in tr_ids_seen_dic, "transcript ID %s in exon feature not (yet) encountered. Transcript feature expected to come before CDS features in GTF file \"%s\"" %(tr_id, in_gtf)
        
            if tr_ids_dic is not None:
                if tr_id not in tr_ids_dic:
                    continue

            m = re.search('transcript_id "(.+?)"', infos)
            assert m, "transcript_id entry missing in GTF file \"%s\", line \"%s\"" %(in_gtf, line)    
            tr_id = m.group(1)

            assert tr_id in tr_ids_seen_dic, "transcript ID %s in exon feature not (yet) encountered. Transcript feature expected to come before exon features in GTF file \"%s\"" %(tr_id, in_gtf)

            # Store/update CDS genomic start and end coordinates.
            if tid2tio_dic[tr_id].cds_s is None:
                tid2tio_dic[tr_id].cds_s = feat_s
                tid2tio_dic[tr_id].cds_e = feat_e
            else:
                if tid2tio_dic[tr_id].cds_s > feat_s:
                    tid2tio_dic[tr_id].cds_s = feat_s
                if tid2tio_dic[tr_id].cds_e < feat_e:
                    tid2tio_dic[tr_id].cds_e = feat_e

    f.close()

    if empty_check:
        assert tid2tio_dic, "no transcript infos read in from GTF file \"%s\"" %(in_gtf)

    return tid2tio_dic


################################################################################

def output_transcript_info_intron_exon_to_bed(tid2tio_dic, out_bed,
                                              output_mode=1,
                                              report_counts=True,
                                              add_tr_id=True,
                                              empty_check=False):
    """
    Output exon and/or intron regions to BED given a dictionary of TranscriptInfo 
    objects (tid2tio_dic), with key transcript ID and value TranscriptInfo object.
    
    output_mode:
        1: output exon + intron regions to BED.
        2: output exon regions only.
        3: output intron regions only.
    add_tr_id:
        If True, add transcript ID after exon intron label (BED column 4), 
        format: intron;tr_id
        
    """
    assert tid2tio_dic, "given tid2tio_dic empty"

    OUTBED = open(out_bed, "w")

    c_exon_out = 0
    c_intron_out = 0

    for tr_id in tid2tio_dic:
        
        chr_id = tid2tio_dic[tr_id].chr_id
        tr_pol = tid2tio_dic[tr_id].tr_pol
        exon_coords = tid2tio_dic[tr_id].exon_coords
        assert exon_coords is not None, "exon coordinates list not set for transcript ID %s" %(tr_id)
    
        # Get intron coordinates.
        intron_coords = []
        if tr_pol == "+":
            for i in range(len(exon_coords) - 1):
                intron_coords.append([exon_coords[i][1]+1, exon_coords[i+1][0]-1])
        elif tr_pol == "-":
            for i in range(len(exon_coords) - 1):
                intron_coords.append([exon_coords[i+1][1]+1,exon_coords[i][0]-1])
        else:
            assert False, "invalid strand given (%s) for transcript ID %s" %(tr_pol, tr_id)

        if output_mode == 1:
            for exon in exon_coords:
                c_exon_out += 1
                if add_tr_id:
                    OUTBED.write("%s\t%i\t%i\texon;%s\t0\t%s\n" % (chr_id, exon[0]-1, exon[1], tr_id, tr_pol))
                else:
                    OUTBED.write("%s\t%i\t%i\texon\t0\t%s\n" % (chr_id, exon[0]-1, exon[1], tr_pol))
            for intron in intron_coords:
                c_intron_out += 1
                if add_tr_id:
                    OUTBED.write("%s\t%i\t%i\tintron;%s\t0\t%s\n" % (chr_id, intron[0]-1, intron[1], tr_id, tr_pol))
                else:
                    OUTBED.write("%s\t%i\t%i\tintron\t0\t%s\n" % (chr_id, intron[0]-1, intron[1], tr_pol))
        elif output_mode == 2:
            for exon in exon_coords:
                c_exon_out += 1
                if add_tr_id:
                    OUTBED.write("%s\t%i\t%i\texon;%s\t0\t%s\n" % (chr_id, exon[0]-1, exon[1], tr_id, tr_pol))
                else:
                    OUTBED.write("%s\t%i\t%i\texon\t0\t%s\n" % (chr_id, exon[0]-1, exon[1], tr_pol))
        elif output_mode == 3:
            for intron in intron_coords:
                c_intron_out += 1
                if add_tr_id:
                    OUTBED.write("%s\t%i\t%i\tintron;%s\t0\t%s\n" % (chr_id, intron[0]-1, intron[1], tr_id, tr_pol))
                else:
                    OUTBED.write("%s\t%i\t%i\tintron\t0\t%s\n" % (chr_id, intron[0]-1, intron[1], tr_pol))

    OUTBED.close()

    if report_counts:
        print("# of exon features output to BED:   ", c_exon_out)
        print("# of intron features output to BED: ", c_intron_out)
    if empty_check:
        assert (c_exon_out+c_intron_out) > 0, "no exon/intron features output to BED"


################################################################################

def bed_intersect_files(a_bed, b_bed, out_bed,
                                    params="-s -u"):

    """
    Intersect two files using bedtools intersect.

    """
    assert os.path.exists(a_bed), "a_bed does not exist"
    assert os.path.exists(b_bed), "b_bed does not exist"

    check_cmd = "intersectBed -a " + a_bed + " -b " + b_bed + " " + params + " > " + out_bed
    output = subprocess.getoutput(check_cmd)
    error = False
    if output:
        error = True
    assert error == False, "intersectBed has problems with your input:\n%s\n%s" %(check_cmd, output)


################################################################################

class TranscriptInfo:
    """
    Store transcript infos.

    Store exons in order, i.e., for plus strand exon 1 is most upstream,
    while for minus strand exon 1 is most downstream.


    """
    def __init__(self,
                 tr_id: str,
                 tr_biotype: str,
                 chr_id: str,
                 tr_s: int,  # 1-based index (GTF).
                 tr_e: int,  # 1-based index.
                 tr_pol: str,
                 gene_id: str,
                 tr_length: Optional[int] = None,
                 exon_c: Optional[int] = None,
                 cds_s: Optional[int] = None,
                 cds_e: Optional[int] = None,
                 exon_coords=None) -> None:

        self.tr_id = tr_id
        self.tr_biotype = tr_biotype
        self.chr_id = chr_id
        self.tr_s = tr_s
        self.tr_e = tr_e
        self.tr_pol = tr_pol
        self.gene_id = gene_id
        self.tr_length = tr_length
        self.exon_c = exon_c
        self.cds_s = cds_s
        self.cds_e = cds_e
        if exon_coords is None:
            self.exon_coords = []
        else:
            self.exon_coords = exon_coords


################################################################################

class GeneInfo:
    """
    Stores gene infos.
    
    GENCODE example (gencode.v44.annotation.gtf.gz):
chr1	HAVANA	gene	11869	14409	.	+	.	gene_id "ENSG00000290825.1"; gene_type "lncRNA"; gene_name "DDX11L2"; level 2; tag "overlaps_pseudogene";
chr1	HAVANA	transcript	11869	14409	.	+	.	gene_id "ENSG00000290825.1"; transcript_id "ENST00000456328.2"; gene_type "lncRNA"; gene_name "DDX11L2"; transcript_type "lncRNA"; transcript_name "DDX11L2-202"; level 2; transcript_support_level "1"; tag "basic"; tag "Ensembl_canonical"; havana_transcript "OTTHUMT00000362751.1";
chr1	HAVANA	exon	11869	12227	.	+	.	gene_id "ENSG00000290825.1"; transcript_id "ENST00000456328.2"; gene_type "lncRNA"; gene_name "DDX11L2"; transcript_type "lncRNA"; transcript_name "DDX11L2-202"; exon_number 1; exon_id "ENSE00002234944.1"; level 2; transcript_support_level "1"; tag "basic"; tag "Ensembl_canonical"; havana_transcript "OTTHUMT00000362751.1";

    Ensembl example (Homo_sapiens.GRCh38.110.gtf.gz):
1	havana	gene	11869	14409	.	+	.	gene_id "ENSG00000290825"; gene_version "1"; gene_name "DDX11L2"; gene_source "havana"; gene_biotype "lncRNA";
1	havana	transcript	11869	14409	.	+	.	gene_id "ENSG00000290825"; gene_version "1"; transcript_id "ENST00000456328"; transcript_version "2"; gene_name "DDX11L2"; gene_source "havana"; gene_biotype "lncRNA"; transcript_name "DDX11L2-202"; transcript_source "havana"; transcript_biotype "lncRNA"; tag "basic"; tag "Ensembl_canonical"; transcript_support_level "1";
1	havana	exon	11869	12227	.	+	.	gene_id "ENSG00000290825"; gene_version "1"; transcript_id "ENST00000456328"; transcript_version "2"; exon_number "1"; gene_name "DDX11L2"; gene_source "havana"; gene_biotype "lncRNA"; transcript_name "DDX11L2-202"; transcript_source "havana"; transcript_biotype "lncRNA"; exon_id "ENSE00002234944"; exon_version "1"; tag "basic"; tag "Ensembl_canonical"; transcript_support_level "1";

    """

    def __init__(self,
                 gene_id: str,
                 gene_name: str,
                 gene_biotype: str,
                 chr_id: str,
                 gene_s: int,  # 1-based index (GTF).
                 gene_e: int,  # 1-based index.
                 gene_pol: str,
                 tr_ids=None,
                 tr_biotypes=None,
                 tr_basic_tags=None,
                 tr_ensembl_canonical_tags=None,
                 tr_lengths=None,
                 tr_tsls=None) -> None:
        self.gene_id = gene_id
        self.gene_name = gene_name
        self.gene_biotype = gene_biotype
        self.chr_id = chr_id
        self.gene_s = gene_s
        self.gene_e = gene_e
        self.gene_pol = gene_pol
        if tr_ids is None:
            self.tr_ids = []
        else:
            self.tr_ids = tr_ids
        if tr_biotypes is None:
            self.tr_biotypes = []
        else:
            self.tr_biotypes = tr_biotypes
        if tr_basic_tags is None:
            self.tr_basic_tags = []
        else:
            self.tr_basic_tags = tr_basic_tags
        if tr_ensembl_canonical_tags is None:
            self.tr_ensembl_canonical_tags = []
        else:
            self.tr_ensembl_canonical_tags = tr_ensembl_canonical_tags
        if tr_lengths is None:
            self.tr_lengths = []
        else:
            self.tr_lengths = tr_lengths
        if tr_tsls is None:
            self.tr_tsls = []
        else:
            self.tr_tsls = tr_tsls


################################################################################

def get_cds_exon_overlap(cds_s, cds_e, exon_s, exon_e, 
                         strand="+"):
    """
    Based on CDS and exon coordinates, get CDS, 5'UTR, and 3'UTR parts of 
    the exon (coordinates).
    
    cds_s:
        CDS start (1-based)
    cds_e:
        CDS end (1-based)
    exon_s:
        Exon start (1-based)
    exon_e:
        Exon end (1-based)
    strand:
        Set strand of exon/cds coordinates. "+" (plus) or "-" (minus) strand
    
    Returns:
    CDS overlap coordinates (False if no overlap with exon),
    Upstream (5'UTR) overlap coordinates (False if no 5'UTR part),
    Upstream (3'UTR) overlap coordinates (False if no 3'UTR part),
    CDS label (False if not present in exon),
    5'UTR label (False if not present in exon),
    3'UTR label (False if not present in exon)

    >>> cds_s, cds_e = 200, 300
    >>> exon_s, exon_e = 150, 250
    >>> get_cds_exon_overlap(cds_s, cds_e, exon_s, exon_e)
    ((200, 250), (150, 199), False, 'CDS', "5'UTR", False)
    >>> cds_s, cds_e = 200, 220
    >>> exon_s, exon_e = 150, 250
    >>> get_cds_exon_overlap(cds_s, cds_e, exon_s, exon_e)
    ((200, 220), (150, 199), (221, 250), 'CDS', "5'UTR", "3'UTR")
    >>> cds_s, cds_e = 50, 100
    >>> exon_s, exon_e = 150, 250
    >>> get_cds_exon_overlap(cds_s, cds_e, exon_s, exon_e)
    (False, False, (150, 250), False, False, "3'UTR")
    >>> cds_s, cds_e = 300, 400
    >>> exon_s, exon_e = 150, 250
    >>> get_cds_exon_overlap(cds_s, cds_e, exon_s, exon_e)
    (False, (150, 250), False, False, "5'UTR", False)
    >>> cds_s, cds_e = 100, 200
    >>> exon_s, exon_e = 150, 250
    >>> get_cds_exon_overlap(cds_s, cds_e, exon_s, exon_e, strand="-")
    ((150, 200), (201, 250), False, 'CDS', "5'UTR", False)
    >>> cds_s, cds_e = 200, 300
    >>> exon_s, exon_e = 150, 250
    >>> get_cds_exon_overlap(cds_s, cds_e, exon_s, exon_e, strand="-")
    ((200, 250), False, (150, 199), 'CDS', False, "3'UTR")
    >>> cds_s, cds_e = 300, 400
    >>> exon_s, exon_e = 150, 250
    >>> get_cds_exon_overlap(cds_s, cds_e, exon_s, exon_e, strand="-")
    (False, False, (150, 250), False, False, "3'UTR")
    >>> cds_s, cds_e = 50, 100
    >>> exon_s, exon_e = 150, 250
    >>> get_cds_exon_overlap(cds_s, cds_e, exon_s, exon_e, strand="-")
    (False, (150, 250), False, False, "5'UTR", False)
    >>> cds_s, cds_e = 200, 220
    >>> exon_s, exon_e = 150, 250
    >>> get_cds_exon_overlap(cds_s, cds_e, exon_s, exon_e, strand="-")
    ((200, 220), (221, 250), (150, 199), 'CDS', "5'UTR", "3'UTR")


    """
    assert strand == "+" or strand == "-", "invalid strand set"

    overlap_start = max(cds_s, exon_s)
    overlap_stop = min(cds_e, exon_e)
    
    # If there is CDS-exon overlap.
    if overlap_start <= overlap_stop:
        non_overlap1 = (exon_s, overlap_start-1) if cds_s > exon_s else False
        non_overlap2 = (overlap_stop+1, exon_e) if exon_e > cds_e else False
        
        # Negative strand, just switch 5'UTR and 3'UTR parts.
        if strand == "-":
            no1 = non_overlap1
            non_overlap1 = non_overlap2
            non_overlap2 = no1

        if non_overlap1 and non_overlap2:
            return (overlap_start, overlap_stop), non_overlap1, non_overlap2, "CDS", "5'UTR", "3'UTR"
        elif non_overlap1 and not non_overlap2:
            return (overlap_start, overlap_stop), non_overlap1, non_overlap2, "CDS", "5'UTR", False
        elif not non_overlap1 and non_overlap2:
            return (overlap_start, overlap_stop), non_overlap1, non_overlap2, "CDS", False, "3'UTR"
        else:
            return (overlap_start, overlap_stop), non_overlap1, non_overlap2, "CDS", False, False

    else:  # If no CDS-exon overlap.
        if cds_e <= exon_s:
            if strand == "+":
                return False, False, (exon_s, exon_e), False, False, "3'UTR"
            else:
                return False, (exon_s, exon_e), False, False, "5'UTR", False
        else:
            if strand == "+":
                return False, (exon_s, exon_e), False, False, "5'UTR", False
            else:
                return False, False, (exon_s, exon_e), False, False, "3'UTR"


################################################################################

def output_exon_annotations(tid2tio_dic, out_bed,
                            custom_annot_dic=None,
                            append=False):
    """
    Output detailed exon annotations to BED file (out_bed.
    If transcript has CDS, output CDS / 5'UTR / 3'UTR annotations.
    If not, use biotype from a dictionary of valid types. If biotype not 
    present in dictionary, use other_annot label.

    append:
        If set append to set BED file, do not overwrite it.
    custom_annot_dic:
        Custom annotation dictionary, overwrites valid_annot_dic. Keys are
        expected to be valid transcript biotype strings found in GTF file. 
        
    """

    # Valid RNA annotations.
    valid_annot_dic = {
        "lncRNA" : "lncRNA",
        "snRNA" : "snRNA",
        "miRNA" : "miRNA",
        "snoRNA" : "snoRNA",
        "rRNA" : "rRNA",
        "Mt_rRNA" : "rRNA",
        "scaRNA" : "scaRNA",
        "processed_pseudogene" : "pseudogene",
        "unprocessed_pseudogene" : "pseudogene",
        "transcribed_unprocessed_pseudogene" : "pseudogene",
        "transcribed_processed_pseudogene" : "pseudogene",
        "rRNA_pseudogene" : "pseudogene",
        "IG_V_pseudogene" : "pseudogene",
        "transcribed_unitary_pseudogene" : "pseudogene",
        "unitary_pseudogene" : "pseudogene",
        "polymorphic_pseudogene" : "pseudogene",
        "TR_V_pseudogene" : "pseudogene",
        "pseudogene" : "pseudogene",
        "IG_C_pseudogene" : "pseudogene",
        "TR_J_pseudogene" : "pseudogene",
        "translated_processed_pseudogene" : "pseudogene",
        "translated_unprocessed_pseudogene" : "pseudogene",
        "IG_pseudogene" : "pseudogene",
        "IG_C_gene" : "IG gene",
        "IG_D_gene" : "IG gene",
        "IG_J_gene" : "IG gene",
        "IG_V_gene" : "IG gene",
        "TR_C_gene" : "TR gene",
        "TR_D_gene" : "TR gene",
        "TR_J_gene" : "TR gene",
        "TR_V_gene" : "TR gene"
    }
    if custom_annot_dic is not None:
        valid_annot_dic = custom_annot_dic

    # If transcript biotype is not in valid_annot_dic, use this annotation/ label:
    other_annot = "other (nc)RNA"

    write_op = "w"
    if append:
        write_op = "a"

    OUTEXAN = open(out_bed, write_op)

    for tid in tid2tio_dic:
        tio = tid2tio_dic[tid]

        if tio.cds_s is not None:
            for exon in tio.exon_coords:
                cds_se, utr5_se, utr3_se, cds, utr5, utr3 = get_cds_exon_overlap(tio.cds_s, tio.cds_e, exon[0], exon[1], strand=tio.tr_pol)
                if cds:
                    OUTEXAN.write("%s\t%i\t%i\tCDS;%s\t0\t%s\n" %(tio.chr_id, cds_se[0]-1, cds_se[1], tio.tr_id, tio.tr_pol))
                if utr5:
                    OUTEXAN.write("%s\t%i\t%i\t5'UTR;%s\t0\t%s\n" %(tio.chr_id, utr5_se[0]-1, utr5_se[1], tio.tr_id, tio.tr_pol))
                if utr3:
                    OUTEXAN.write("%s\t%i\t%i\t3'UTR;%s\t0\t%s\n" %(tio.chr_id, utr3_se[0]-1, utr3_se[1], tio.tr_id, tio.tr_pol))
        else:
            label = other_annot
            if tio.tr_biotype in valid_annot_dic:
                label = valid_annot_dic[tio.tr_biotype]
            for exon in tio.exon_coords:
                OUTEXAN.write("%s\t%i\t%i\t%s;%s\t0\t%s\n" %(tio.chr_id, exon[0]-1, exon[1], label, tio.tr_id, tio.tr_pol))

    OUTEXAN.close()


"""

$ gtf_extract_transcript_biotypes.py Homo_sapiens.GRCh38.106.gtf.gz 

protein_coding	87777
lncRNA	51298
retained_intron	32801
processed_transcript	30341
nonsense_mediated_decay	20247
processed_pseudogene	10151
unprocessed_pseudogene	2609
misc_RNA	2220
snRNA	1910
miRNA	1877
TEC	1147
transcribed_unprocessed_pseudogene	953
snoRNA	943
transcribed_processed_pseudogene	505
rRNA_pseudogene	496
IG_V_pseudogene	187
transcribed_unitary_pseudogene	151
IG_V_gene	145
TR_V_gene	106
non_stop_decay	99
unitary_pseudogene	96
TR_J_gene	79
polymorphic_pseudogene	71
rRNA	53
scaRNA	49
IG_D_gene	37
TR_V_pseudogene	33
IG_C_gene	23
Mt_tRNA	22
pseudogene	19
IG_J_gene	18
IG_C_pseudogene	9
ribozyme	8
TR_C_gene	6
sRNA	5
TR_D_gene	4
TR_J_pseudogene	4
IG_J_pseudogene	3
translated_processed_pseudogene	2
translated_unprocessed_pseudogene	2
Mt_rRNA	2
scRNA	1
vault_RNA	1
IG_pseudogene	1


$ gtf_extract_transcript_biotypes.py gencode.v44.annotation.gtf.gz 
protein_coding	89067
lncRNA	56357
retained_intron	34128
protein_coding_CDS_not_defined	26483
nonsense_mediated_decay	21384
processed_pseudogene	10147
unprocessed_pseudogene	2610
misc_RNA	2212
snRNA	1901
miRNA	1879
processed_transcript	1352
TEC	1145
transcribed_unprocessed_pseudogene	961
snoRNA	942
transcribed_processed_pseudogene	512
rRNA_pseudogene	497
IG_V_pseudogene	187
transcribed_unitary_pseudogene	155
IG_V_gene	145
TR_V_gene	107
non_stop_decay	103
unitary_pseudogene	99
TR_J_gene	79
protein_coding_LoF	73
scaRNA	49
rRNA	47
IG_D_gene	37
TR_V_pseudogene	33
IG_C_gene	23
Mt_tRNA	22
artifact	19
IG_J_gene	18
pseudogene	15
IG_C_pseudogene	9
ribozyme	8
TR_C_gene	6
sRNA	5
TR_D_gene	5
TR_J_pseudogene	4
IG_J_pseudogene	3
translated_processed_pseudogene	2
Mt_rRNA	2
scRNA	1
vault_RNA	1
IG_pseudogene	1

"""

################################################################################

def get_region_annotations(overlap_annotations_bed,
                           reg_ids_dic=None):
    """
    Get region annotations from overlapping genomic regions with exon / intron 
    / transript biotype annotations.

    reg_ids_dic:
        If set, compare genomic region IDs with IDs in dictionary. If region ID 
        from dictionary not in overlap_annotations_bed, set label "intergenic".

    """

    assert os.path.exists(overlap_annotations_bed), "file %s does not exist" %(overlap_annotations_bed)

    reg2maxol_dic = {}
    reg2annot_dic = {}

    with open(overlap_annotations_bed) as f:
        for line in f:
            cols = line.strip().split("\t")
            reg_id = cols[3]
            annot_ids = cols[9].split(";")
            assert len(annot_ids) == 2, "len(annot_ids) != 2 (expected ; separated string, but got: \"%s\")" %(cols[9])
            annot_id = annot_ids[0]
            tr_id = annot_ids[1]
            c_overlap_nts = cols[12]
            if reg_id not in reg2maxol_dic:
                reg2maxol_dic[reg_id] = c_overlap_nts
                reg2annot_dic[reg_id] = [annot_id, tr_id]
            else:
                if c_overlap_nts > reg2maxol_dic[reg_id]:
                    reg2maxol_dic[reg_id] = c_overlap_nts
                    reg2annot_dic[reg_id][0] = annot_id
                    reg2annot_dic[reg_id][1] = tr_id
    f.closed

    if reg_ids_dic is not None:
        for reg_id in reg_ids_dic:
            if reg_id not in reg2annot_dic:
                reg2annot_dic[reg_id] = ["intergenic", False]

    return reg2annot_dic


"""
$ intersectBed -a regions.bed -b introns.bed -s -wo -f 0.5
chr1	1980	2020	s1	0	+	chr1	1000	2000	exon	0	+	20
chr1	1980	2020	s1	0	+	chr1	2000	3000	intron	0	+	20

"""


################################################################################

def get_hex_colors_list(min_len=0):
    """
    Get list of distinct hex color values for plotting.

    matplotlib color maps:
    https://matplotlib.org/stable/gallery/color/colormap_reference.html

    How to get from colormap to hex colors:

    import matplotlib.pyplot as plt
    import matplotlib
    import numpy as np

    cmap = plt.cm.get_cmap('tab20b') # colormap tab20b
    colors = cmap(np.linspace(0,1,cmap.N))

    hex_colors = []
    for color in colors:
        hex_color = matplotlib.colors.rgb2hex(color)
        hex_colors.append(hex_color)

    """

    assert min_len < 1000, "provide reasonable min_len"

    hex_colors_tab20_list = ['#1f77b4', '#aec7e8', '#ff7f0e', '#ffbb78', '#2ca02c', '#98df8a', '#d62728', '#ff9896', '#9467bd', '#c5b0d5', '#8c564b', '#c49c94', '#e377c2', '#f7b6d2', '#7f7f7f', '#c7c7c7', '#bcbd22', '#dbdb8d', '#17becf', '#9edae5']
    hex_colors_tab20b_list = ['#393b79', '#5254a3', '#6b6ecf', '#9c9ede', '#637939', '#8ca252', '#b5cf6b', '#cedb9c', '#8c6d31', '#bd9e39', '#e7ba52', '#e7cb94', '#843c39', '#ad494a', '#d6616b', '#e7969c', '#7b4173', '#a55194', '#ce6dbd', '#de9ed6']
    hex_colors_Pastel1_list = ['#fbb4ae', '#b3cde3', '#ccebc5', '#decbe4', '#fed9a6', '#ffffcc', '#e5d8bd', '#fddaec', '#f2f2f2']

    hex_colors_list = []
    for hc in hex_colors_tab20_list:
        hex_colors_list.append(hc)
    for hc in hex_colors_tab20b_list:
        hex_colors_list.append(hc)
    for hc in hex_colors_Pastel1_list:
        hex_colors_list.append(hc)

    while len(hex_colors_list) < min_len:
        hex_colors_list = hex_colors_list*2

    return hex_colors_list


################################################################################

def read_ids_into_dic(ids_file,
                      check_dic=True,
                      ids_dic=False):
    """
    Read in IDs file, where each line stores one ID.

    >>> test_ids_file = "test_data/test.ids"
    >>> ids_dic = read_ids_into_dic(test_ids_file)
    >>> print(ids_dic)
    {'clip1': 1, 'clip2': 1, 'clip3': 1}

    """
    if not ids_dic:
        ids_dic = {}
    # Read in file content.
    with open(ids_file) as f:
        for line in f:
            row_id = line.strip()
            ids_dic[row_id] = 1
    f.closed
    if check_dic:
        assert ids_dic, "IDs dictionary ids_dic empty"
    return ids_dic


################################################################################

def gtf_check_exon_order(in_gtf):
    """
    Check exon_number ordering. Return True if ordering for minus strand
    and plus strand is different (i.e. for minus exon_number 1 is most downstream).
    Return False, if upstream to downstream order numbering is used for both
    plus and minus strand transcripts.

    >>> test_true_gtf = "test_data/test_order_true.gtf"
    >>> test_false_gtf = "test_data/test_order_false.gtf"
    >>> gtf_check_exon_order(test_true_gtf)
    True
    >>> gtf_check_exon_order(test_false_gtf)
    False

    """
    tr2exon_nr_dic = {}
    tr2exon_s_dic = {}

    check = 6666

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
        feat_s = int(cols[3])
        feat_e = int(cols[4])
        feat_pol = cols[6]
        infos = cols[8]

        if feature != "exon":
            continue

        # Transcript ID.
        m = re.search('transcript_id "(.+?)"', infos)
        assert m, "transcript_id entry missing in GTF file \"%s\", line \"%s\"" %(in_gtf, line)
        transcript_id = m.group(1)
        # Exon number.
        m = re.search('exon_number "(\d+?)"', infos)
        # Try GENCODE encoding.
        if not m:
            m = re.search('exon_number (\d+?);', infos)
        assert m, "exon_number entry missing in GTF file \"%s\", line \"%s\"" %(in_gtf, line)
        exon_nr = int(m.group(1))

        # Check whether exon numbers are incrementing for each transcript ID.
        if transcript_id not in tr2exon_nr_dic:
            tr2exon_nr_dic[transcript_id] = exon_nr
        else:
            assert tr2exon_nr_dic[transcript_id] < exon_nr, "transcript ID \"%s\" without increasing exon number order in GTF file \"%s\"" %(transcript_id, in_gtf)
            tr2exon_nr_dic[transcript_id] = exon_nr

        # Check ordering of exons for minus strand transcripts.
        if transcript_id not in tr2exon_s_dic:
            tr2exon_s_dic[transcript_id] = feat_s
        else:
            if feat_pol == "-":
                if tr2exon_s_dic[transcript_id] > feat_s:
                    check = True
                else:
                    check = False
                break
            elif feat_pol == "+":
                assert tr2exon_s_dic[transcript_id] < feat_s, "transcript ID \"%s\" on plus strand but exon region coordinates are decreasing" %(transcript_id)
    f.close()

    assert check != 6666, "no minus strand exon regions found in GTF file %s" %(in_gtf)
    return check


################################################################################

def select_mpts_from_gene_infos(gid2gio_dic,
                                basic_tag=True,
                                ensembl_canonical_tag=False,
                                only_tsl=False,
                                tr_min_len=False
                                ):
    """
    Select most prominent transcripts from GeneInfo objects. 
    Selection is based on transcript support level (TSL), 
    with lower value == better, as well as preference for basic tag + 
    Ensembl canonical.
    If several transcripts with lowest TSL / same tags, select longest transcript.
    Depending on filter settings, there might be genes for which no TSL 
    is reported.
    
    basic_tag:
        If True only report transcripts with "basic" tag.
    ensembl_canonical_tag:
        If True only report transcripts with "Ensembl_canonical" tag.
    only_tsl:
        If True only report transcripts with TSL 1-5 (excluding "NA").
    tr_min_len:
        If length set, only report transcripts with length >= tr_min_len.


    >>> test_gtf = "test_data/test_mpt_selection.gtf"
    >>> gid2gio_dic = gtf_read_in_gene_infos(test_gtf)
    >>> tr_ids_dic = select_mpts_from_gene_infos(gid2gio_dic, basic_tag=False, ensembl_canonical_tag=False, only_tsl=False)
    >>> tr_ids_dic
    {'ENST00000357266': 'ENSG00000096060'}

    """

    # Comparison dictionary.
    id2sc = {}
    for i in range(5):
        pos = i + 1
        pos_str = "%i" %(pos)
        id2sc[pos_str] = pos
    id2sc["NA"] = 6

    mpt2gid_dic = {}

    for gene_id in gid2gio_dic:
        gene_info = gid2gio_dic[gene_id]
        mpt_id = "-"
        mpt_tsl = "NA"
        mpt_len = 0
        mpt_bt = 0
        mpt_ec = 0

        for idx, tr_id in enumerate(gene_info.tr_ids):
            # print("mpt_id:", mpt_id, "tr_id:", tr_id)
            tr_tsl = gene_info.tr_tsls[idx]  # 1-5 or NA
            tr_bt = gene_info.tr_basic_tags[idx]  # 0 or 1
            tr_ec = gene_info.tr_ensembl_canonical_tags[idx]  # 0 or 1
            tr_length = gene_info.tr_lengths[idx]
            # print(tr_id, tr_bt, tr_ec, tr_length, tr_tsl)
            if basic_tag:
                if not tr_bt:
                    continue
            if ensembl_canonical_tag:
                if not tr_ec:
                    continue
            if only_tsl:
                if tr_tsl == "NA":
                    continue
            if tr_min_len:
                if tr_length < tr_min_len:
                    continue

            if id2sc[tr_tsl] < id2sc[mpt_tsl]:
                mpt_id = tr_id
                mpt_tsl = tr_tsl
                mpt_len = tr_length
                mpt_bt = tr_bt
                mpt_ec = tr_ec
            elif id2sc[tr_tsl] == id2sc[mpt_tsl]:
                # print("Now equal, comparing tr_id %s with mpt_id %s" %(tr_id, mpt_id))
                # If transcript has basic tag, use this.
                if tr_bt > mpt_bt:
                    mpt_id = tr_id
                    mpt_tsl = tr_tsl
                    mpt_len = tr_length
                    mpt_bt = tr_bt
                    mpt_ec = tr_ec
                    continue
                # If transcript has Ensembl canonical tag, use this.
                if tr_ec > mpt_ec:
                    mpt_id = tr_id
                    mpt_tsl = tr_tsl
                    mpt_len = tr_length
                    mpt_bt = tr_bt
                    mpt_ec = tr_ec
                    continue
                # If same basic/Ensembl canonical tag combination.
                if tr_ec == mpt_ec and tr_bt == mpt_bt:
                    if tr_length > mpt_len:
                        mpt_id = tr_id
                        mpt_tsl = tr_tsl
                        mpt_len = tr_length
                        mpt_bt = tr_bt
                        mpt_ec = tr_ec

        if not mpt_len:
            continue

        mpt2gid_dic[mpt_id] = gene_id

    return mpt2gid_dic


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

    RBPBench v0.2:
    Currently ignore Organism column / do not use this information.

    """
    name2ids_dic = {}
    id2type_dic = {}
    # id2org_dic = {}

    with open(rbp2ids_file) as f:
        for line in f:
            if re.search("^RBP_motif_id", line) or re.search("^#", line):
                continue
            cols = line.strip().split("\t")
            motif_id = cols[0]
            rbp_name = cols[1]
            motif_type = cols[2]
            # organism = cols[3]            
            # if organism != "human":
            #     rbp_name = rbp_name + "_" + organism
            # id2org_dic[motif_id] = organism

            id2type_dic[motif_id] = motif_type
            if rbp_name in name2ids_dic:
                name2ids_dic[rbp_name].append(motif_id)
            else:
                name2ids_dic[rbp_name] = [motif_id]
    f.closed

    assert name2ids_dic, "no RBP IDs read in from %s" %(rbp2ids_file)

    return name2ids_dic, id2type_dic  #, id2org_dic


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

def extract_pol_from_seq_ids(out_seqs_dic):
    """
    Extract strand/polarity info from FASTA sequence header ID.

    out_seqs_dic:
        FASTA header ID -> sequence.
        Header ID format:
        chr20:62139082-62139128(-)

    """
    reg2pol_dic = {}
    for seq_id in out_seqs_dic:
        if re.search("\w+:\d+-\d+\([+|-]\)", seq_id):
            m = re.search("\w+:\d+-\d+\(([+|-])\)", seq_id)
            reg2pol_dic[seq_id] = m.group(1)
        else:
            assert False, "region ID has invalid format (%s). Please contact developers" %(reg_id)
    return reg2pol_dic


################################################################################

def bed_filter_extend_bed(in_bed, out_bed,
                          ext_up=0,
                          ext_down=0,
                          remove_dupl=True,
                          reg2sc_dic=None,
                          score_col=5,
                          chr_ids_dic=None,
                          bed_chr_ids_dic=None,
                          use_region_ids=False,
                          unstranded=False):
    """
    Filter / extend in_bed BED file, output to out_bed.

    remove_dupl:
        Remove duplicated regions, i.e., regions with same 
        chr_id, reg_s, reg_e, reg_pol
    use_region_ids:
        If True set column 4 region ID to location ID, in style: "chr20:62139082-62139128(-)".
        This way we get FASTA headers == BED col4 IDs. 
    unstranded:
        If True, output both strands of each region (+ and -), so given one 
        input region, two output regions are generated.
    reg2sc_dic:
        region ID -> column score_col BED score

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

            # Record present IDs.
            if bed_chr_ids_dic is not None:
                bed_chr_ids_dic[chr_id] = 1

            # Use "chr20:62139082-62139128(-)" style.
            # if use_region_ids:
            #     if unstranded:
            #         reg_id = "%s:%i-%i" %(chr_id, reg_s, reg_e)
            #     else:
            #         reg_id = "%s:%i-%i(%s)" %(chr_id, reg_s, reg_e, reg_pol)

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

                if use_region_ids:
                    OUTBED.write("%s\t%i\t%i\t%s:%i-%i(+)\t%s\t+\n" % (chr_id, new_s, new_e, chr_id, new_s, new_e, reg_sc))
                    OUTBED.write("%s\t%i\t%i\t%s:%i-%i(-)\t%s\t-\n" % (chr_id, new_s, new_e, chr_id, new_s, new_e, reg_sc))
                else:
                    OUTBED.write("%s\t%i\t%i\t%s\t%s\t+\n" % (chr_id, new_s, new_e, reg_id, reg_sc))
                    OUTBED.write("%s\t%i\t%i\t%s\t%s\t-\n" % (chr_id, new_s, new_e, reg_id, reg_sc)) 

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

                if use_region_ids:
                    OUTBED.write("%s\t%i\t%i\t%s:%i-%i(%s)\t%s\t%s\n" % (chr_id, new_s, new_e, chr_id, new_s, new_e, reg_pol, reg_sc, reg_pol))
                else:
                    OUTBED.write("%s\t%i\t%i\t%s\t%s\t%s\n" % (chr_id, new_s, new_e, reg_id, reg_sc, reg_pol))

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

def convert_sci_not_to_decimal(sn_value):
    """
    Convert scientific notation value to decimal string.

    >>> convert_sci_not_to_decimal(6.879779511865709E-8)
    '0.00000006879779511865709'

    """
    dec_str = np.format_float_positional(sn_value, trim='-')
    return dec_str


################################################################################

def search_generate_html_report(df_pval, pval_cont_lll,
                                search_rbps_dic,
                                id2name_dic, name2ids_dic,
                                region_rbp_motif_pos_dic,
                                reg2pol_dic,
                                out_folder, 
                                benchlib_path,
                                rbp2regidx_dic,
                                reg_ids_list,
                                set_rbp_id=None,
                                motif_db_str=False,
                                seq_motif_blocks_dic=None,
                                motif_distance_plot_range=50,
                                motif_min_pair_count=10,
                                rbp_min_pair_count=10,
                                reg2annot_dic=None,
                                upset_plot_min_degree=2,
                                upset_plot_max_degree=None,
                                upset_plot_min_subset_size=10,
                                html_report_out="report.rbpbench_search.html",
                                plot_abs_paths=False,
                                plotly_js_source=1,
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

    # Number of genomic regions.
    c_regions = len(reg_ids_list)

    # Makes tables sortable.
    sorttable_js_path = benchlib_path + "/content/sorttable.js"
    # plotly js path.
    # https://plotly.com/javascript/getting-started/#download
    # plotly-2.20.0.min.js
    # plotly-latest.min.js
    plotly_js_path = "cdn"
    if plotly_js_source == 1:
       plotly_js_path = "cdn"
    elif plotly_js_source == 2:
        plotly_js_path = benchlib_path + "/content/plotly-2.20.0.min.js"
        assert os.path.exists(plotly_js_path), "plotly js %s not found" %(plotly_js_path)
    else:
        assert False, "invalid plotly_js_path set"

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
- [RBP combinations upset plot](#rbp-comb-upset-plot)"""

    # Additional plot if GTF annotations given.
    if reg2annot_dic is not None:
        mdtext += "\n"
        mdtext += "- [Region annotations per RBP](#annot-rbp-plot)\n"
    # If --set-rbp-id given.
    if set_rbp_id is not None:
        if reg2annot_dic is None:
            mdtext += "\n"
        mdtext += "- [Set RBP %s motifs distance statistics](#rbp-motif-dist-stats)\n" %(set_rbp_id)
        # mdtext += "- [Set RBP %s motif distance plot](#rbp-motif-dist-plot)\n" %(set_rbp_id)
        for idx, motif_id in enumerate(name2ids_dic[set_rbp_id]):
            mdtext += "    - [Motif %s distance statistics](#single-motif-%i-dist-stats)\n" %(motif_id, idx)
            # mdtext += "    - [Single motif %s distance plot](#single-motif-%i-dist-plot)\n" %(motif_id, idx)

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
        wc_pval_str = convert_sci_not_to_decimal(wc_pval)  # Convert scientific notation to decimal string for sorting to work.
        c_hit_reg = search_rbps_dic[rbp_id].c_hit_reg
        perc_hit_reg = search_rbps_dic[rbp_id].perc_hit_reg
        c_uniq_motif_hits = search_rbps_dic[rbp_id].c_uniq_motif_hits
        mdtext += "| %s | %i | %.2f | %i | %s |\n" %(rbp_id, c_hit_reg, perc_hit_reg, c_uniq_motif_hits, wc_pval_str)
    mdtext += "\n&nbsp;\n&nbsp;\n"
    mdtext += "\nColumn IDs have the following meanings: "
    mdtext += "**RBP ID** -> RBP ID from database or user-defined (typically RBP name), "
    mdtext += '**# hit regions** -> number of input genomic regions with motif hits (after filtering and optional extension), '
    mdtext += '**% hit regions** -> percentage of hit regions over all regions (i.e. how many input regions contain >= 1 RBP binding motif), '
    mdtext += '**# motif hits** -> number of unique motif hits in input regions (removed double counts), '
    mdtext += '**p-value** -> Wilcoxon rank-sum test p-value.' + "\n"
    mdtext += "\n&nbsp;\n"

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

"""
    mdtext += '<div class=class="container-fluid" style="margin-top:40px">' + "\n"
    mdtext += '<iframe src="' + plot_path + '" width="1200" height="1200"></iframe>' + "\n"
    mdtext += '</div>'
    mdtext += """

**Figure:** Heat map of co-occurrences (Fisher's exact test p-values) between RBPs. 
Legend color: negative logarithm (base 10) of Fisher's exact test p-value.
Hover box: 1) RBP1. 2) RBP2. 3) p-value: Fisher's exact test p-value (calculated based on contingency table between RBP1 and RBP2). 
4) RBPs compaired. 5) Counts[]: contingency table of co-occurrence counts (i.e., number of genomic regions with/without shared motif hits) between compaired RBPs, 
with format [[A, B], [C, D]], where 
A: RBP1 AND RBP2, 
B: NOT RBP1 AND RBP2
C: RBP1 AND NOT RBP2
D: NOT RBP1 AND NOT RBP2.
6) Correlation: Pearson correlation coefficients between RBPs. 
Genomic regions are labelled 1 or 0 (RBP motif present or not), resulting in a vector of 1s and 0s for each RBP.
Correlations are then calculated by comparing vectors for every pair of RBPs.
7) -log10 of p-value, used as color values.

&nbsp;

"""

    """
    RBP region occupancies upset plot.

    """

    rbp_reg_occ_upset_plot =  "rbp_region_occupancies.upset_plot.png"
    rbp_reg_occ_upset_plot_out = plots_out_folder + "/" + rbp_reg_occ_upset_plot

    plotted, reason, count = create_rbp_reg_occ_upset_plot(rbp2regidx_dic, reg_ids_list, 
                                  reg2annot_dic=reg2annot_dic,
                                  min_degree=upset_plot_min_degree,
                                  max_degree=upset_plot_max_degree,
                                  min_subset_size=upset_plot_min_subset_size,
                                  plot_out=rbp_reg_occ_upset_plot_out)


    plot_path = plots_folder + "/" + rbp_reg_occ_upset_plot

    mdtext += """
## RBP combinations upset plot ### {#rbp-comb-upset-plot}

"""

    if plotted:
        mdtext += '<img src="' + plot_path + '" alt="RBP region occupancies upset plot"' + "\n"
        mdtext += 'title="RBP region occupancies upset plot" />' + "\n"
        mdtext += """

**Figure:** Upset plot of RBP combinations found in the given set of genomic regions (# of regions = %i). 
Intersection size == how often a specific RBP combination is found in the regions dataset.
For example, if two regions in the input set contain motif hits for RBP1, RBP3, and RBP5, then the RBP combination RBP1,RBP3,RBP5 will get a count (== Intersection size) of 2.
Minimum occurrence number for a combination to be reported = %i (command line parameter: --upset-plot-min-subset-size). 
How many RBPs a combination must contain to be reported = %i (command line parameter: --upset-plot-min-degree).
The numbers on the left side for each RBP tell how many genomic regions have motif hits (1 or more) of the respective RBP. 
If a GTF file was given, bar charts become stacked bar charts, showing what GTF annotations the regions overlap with (see legend for region types).

&nbsp;

""" %(c_regions, upset_plot_min_subset_size, upset_plot_min_degree)

    else:

        if reason == "min_degree":

            mdtext += """

No upset plot generated since set --upset-plot-min-degree > maximum degree found in the RBP combination set. Please use lower number for --upset-plot-min-degree parameter.

&nbsp;

"""

        elif reason == "min_subset_size":

            mdtext += """

No upset plot generated since set --upset-plot-min-subset-size (%i) > maximum subset size (%i) found in the RBP combination set. Please use lower number for --upset-plot-min-subset-size parameter.

&nbsp;

""" %(upset_plot_min_subset_size, count)

        elif reason == "len(rbp_id_list) == 1":

            mdtext += """

No upset plot generated since number of selected RBPs == 1.

&nbsp;

"""

        else:
            assert False, "invalid reason given for not plotting upset plot"



    if reg2annot_dic is not None:

        annot_stacked_bars_plot =  "annotation_stacked_bars_plot.png"
        annot_stacked_bars_plot_out = plots_out_folder + "/" + annot_stacked_bars_plot


        mdtext += """
## Region annotations per RBP ### {#annot-rbp-plot}

"""
        create_annotation_stacked_bars_plot(rbp2regidx_dic, reg_ids_list, reg2annot_dic,
                                            plot_out=annot_stacked_bars_plot_out)

        plot_path = plots_folder + "/" + annot_stacked_bars_plot

        mdtext += '<img src="' + plot_path + '" alt="Annotation stacked bars plot"' + "\n"
        # mdtext += 'title="Annotation stacked bars plot" width="800" />' + "\n"
        mdtext += 'title="Annotation stacked bars plot" />' + "\n"
        mdtext += """
**Figure:** For each RBP, a stacked bar shows the corresponding region annotations 
(from --gtf GTF file, see legend for region types) for the genomic regions 
with motif hits for the respective RBP. 
Total bar height equals to the number of genomic regions with >= 1 motif hit for the RBP.

&nbsp;

"""

    """
    Set RBP motif distance stats + plot.

    """
    if set_rbp_id is not None:

        rbp_motif_dist_plot_plotly =  "%s.motif_dist.plotly.html" %(set_rbp_id)
        rbp_motif_dist_plot_plotly_out = plots_out_folder + "/" + rbp_motif_dist_plot_plotly

        plotted, pc_dic, in_dic, out_dic = plot_motif_dist_rbp_level(set_rbp_id,
            region_rbp_motif_pos_dic,
            id2name_dic,
            name2ids_dic,
            reg2pol_dic,
            html_out=rbp_motif_dist_plot_plotly_out,
            plotly_js_path=plotly_js_path,
            line_plot_range=motif_distance_plot_range,
            min_pair_count=rbp_min_pair_count)

        mdtext += """
## Set RBP %s motifs distance statistics ### {#rbp-motif-dist-stats}

**Table:** Motif distance statistics between set RBP ID (%s) motifs and other RBP ID motifs (including itself). 
Note that the statistics are generated by focussing (i.e., centering) on the highest-scoring motif of the set RBP 
(lowest p-value for sequence motifs, highest bit score for structure motifs) in each input region.
In case of an empty table, try to lower --rbp-min-pair-count (current value: %i).

""" %(set_rbp_id, set_rbp_id, rbp_min_pair_count)
        mdtext += '| Set RBP ID | Other RBP ID | Pair count | # near motifs | # distant motifs |' + " \n"
        mdtext += "| :-: | :-: | :-: | :-: | :-: |\n"

        top_c = 0
        for rbp_id, pair_c in sorted(pc_dic.items(), key=lambda item: item[1], reverse=True):
            if pair_c >= rbp_min_pair_count:
                top_c += 1
                mdtext += "| %s | %s | %i | %i | %i |\n" %(set_rbp_id, rbp_id, pair_c, in_dic[rbp_id], out_dic[rbp_id])

        mdtext += "\n&nbsp;\n&nbsp;\n"
        mdtext += "\nColumn IDs have the following meanings: "
        mdtext += "**Set RBP ID** -> set RBP ID (specified via --set-rbp-id), "
        mdtext += "**Other RBP ID** -> other RBP ID that set RBP ID is compared to, "
        mdtext += "**Pair count** -> number of input regions with motif hits from both RBP IDs (minimum pair count to be reported: %i (set by --rbp-min-pair-count)), " %(rbp_min_pair_count)
        mdtext += '**# near motifs** -> ' + "number of other RBP ID motifs within specified distance (----motif-distance-plot-range %i) of the centered set RBP ID motifs, " %(motif_distance_plot_range)
        mdtext += '**# distant motifs** -> ' + "number of other RBP ID motifs outside specified distance (----motif-distance-plot-range %i) of the centered set RBP ID motifs." %(motif_distance_plot_range)
        # mdtext += "\n"
        mdtext += "\n&nbsp;\n"

#         mdtext += """
# ## %s motif distance plot ### {#rbp-motif-dist-plot}

# """ %(set_rbp_id)

        if plotted:

            plot_path = plots_folder + "/" + rbp_motif_dist_plot_plotly

            mdtext += '<div class=class="container-fluid" style="margin-top:40px">' + "\n"
            mdtext += '<iframe src="' + plot_path + '" width="1200" height="700"></iframe>' + "\n"
            mdtext += '</div>'
            mdtext += """

**Figure:** Line plot showing motif distances between set RBP ID (%s) motifs (using highest-scoring %s motif for each input region) and other RBP ID motifs (including itself).
Coverage corresponds to positions occupied by RBP motifs. 
The plot is centered on the highest-scoring %s motif (lowest p-value for sequence or highest bit score) for each region containing >= 1 %s motif.
Each RBP with a pair count (definition see table above) of >= %i is shown, and the coverage of the RBP is the accumulation of its individual motif coverages.

&nbsp;

""" %(set_rbp_id, set_rbp_id, set_rbp_id, set_rbp_id, rbp_min_pair_count)

        else:

            mdtext += """

<em>No motif distance plot generated for set RBP %s. Try to lower --rbp-min-pair-count (current value: %i).</em>

&nbsp;

""" %(set_rbp_id, rbp_min_pair_count)

        """
        Set RBP single motif distance stats + plots.
        
        """

        for idx, motif_id in enumerate(name2ids_dic[set_rbp_id]):

            single_motif_dist_plot_plotly =  "%s.motif_dist.plotly.html" %(motif_id)
            single_motif_dist_plot_plotly_out = plots_out_folder + "/" + single_motif_dist_plot_plotly

            plotted, pc_dic, in_dic, out_dic = plot_motif_dist_motif_level(motif_id,
                region_rbp_motif_pos_dic,
                name2ids_dic,
                reg2pol_dic,
                html_out=single_motif_dist_plot_plotly_out,
                plotly_js_path=plotly_js_path,
                line_plot_range=motif_distance_plot_range,
                min_pair_count=motif_min_pair_count)


            mdtext += """
### Motif %s distance statistics ### {#single-motif-%i-dist-stats}

""" %(motif_id, idx)

            # Plot motif (if sequence motif).
            if motif_id in seq_motif_blocks_dic:

                motif_plot = "%s.%s.png" %(set_rbp_id, motif_id)
                motif_plot_out = plots_out_folder + "/" + motif_plot
                plot_path = plots_folder + "/" + motif_plot

                # Check if motif in motif database folder.
                if motif_db_str:
                    db_motif_path = benchlib_path + "/content/%s_motif_plots/%s" %(motif_db_str, motif_plot)
                    if os.path.exists(db_motif_path):
                        shutil.copy(db_motif_path, motif_plot_out)

                if not os.path.exists(motif_plot_out):
                    create_motif_plot(motif_id, seq_motif_blocks_dic,
                                      motif_plot_out)

                mdtext += '<img src="' + plot_path + '" alt="' + "sequence motif plot %s" %(motif_id) + "\n"
                mdtext += 'title="' + "sequence motif plot %s" %(motif_id) + '" width="500" />' + "\n"
                mdtext += """

**Figure:** Sequence motif %s.

&nbsp;

""" %(motif_id)

            mdtext += """

**Table:** Motif distance statistics between motif %s (RBP: %s) and other motifs (from any RBP). 
Note that the statistics are generated by focussing (i.e., centering) on the highest-scoring hit of motif %s 
(lowest p-value for sequence motifs, highest bit score for structure motifs) in each input region.
In case of an empty table, try to lower --motif-min-pair-count (current value: %i).

""" %(motif_id, set_rbp_id, motif_id, motif_min_pair_count)

            mdtext += '| Motif ID | Other motif ID | &nbsp; Other motif ID plot &nbsp; | Pair count | # near motifs | # distant motifs |' + " \n"
            mdtext += "| :-: | :-: | :-: | :-: | :-: | :-: |\n"

            for other_motif_id, pair_c in sorted(pc_dic.items(), key=lambda item: item[1], reverse=True):

                if pair_c >= motif_min_pair_count:

                    # Plot motif (if sequence motif).
                    plot_str = "-"
                    if other_motif_id in seq_motif_blocks_dic:

                        other_rbp_id = id2name_dic[other_motif_id]
                        motif_plot = "%s.%s.png" %(other_rbp_id, other_motif_id)
                        motif_plot_out = plots_out_folder + "/" + motif_plot
                        plot_path = plots_folder + "/" + motif_plot

                        # Check if motif in motif database folder.
                        if motif_db_str:
                            db_motif_path = benchlib_path + "/content/%s_motif_plots/%s" %(motif_db_str, motif_plot)
                            if os.path.exists(db_motif_path):
                                shutil.copy(db_motif_path, motif_plot_out)

                        if not os.path.exists(motif_plot_out):
                            create_motif_plot(other_motif_id, seq_motif_blocks_dic,
                                              motif_plot_out)

                        plot_str = '<image src = "' + plot_path + '" width="300px"></image>'

                    else:
                        print("Motif ID %s not in seq_motif_blocks_dic ... " %(seq_motif_blocks_dic))

                    mdtext += "| %s | %s | %s | %i | %i | %i |\n" %(motif_id, other_motif_id, plot_str, pair_c, in_dic[other_motif_id], out_dic[other_motif_id])

            mdtext += "\n&nbsp;\n&nbsp;\n"
            mdtext += "\nColumn IDs have the following meanings: "
            mdtext += "**Motif ID** -> motif ID (motif belonging to set RBP), "
            mdtext += "**Other motif ID** -> other motif ID that set RBP motif ID is compared to, "
            mdtext += "**Other motif ID plot** -> other motif ID sequence motif plot (if motif is sequence motif), "
            mdtext += "**Pair count** -> number of input regions containing hits for both motifs (minimum pair count to be reported: %i (set by --motif-min-pair-count)), " %(motif_min_pair_count)
            mdtext += '**# near motifs** -> ' + "number of other motifs within specified distance (----motif-distance-plot-range %i) of the centered motif belonging to set RBP, " %(motif_distance_plot_range)
            mdtext += '**# distant motifs** -> ' + "number of other motifs outside specified distance (----motif-distance-plot-range %i) of the centered motif belonging to set RBP." %(motif_distance_plot_range)
            # mdtext += "\n"
            mdtext += "\n&nbsp;\n"

            if plotted:

                plot_path = plots_folder + "/" + single_motif_dist_plot_plotly

                mdtext += '<div class=class="container-fluid" style="margin-top:40px">' + "\n"
                mdtext += '<iframe src="' + plot_path + '" width="1200" height="700"></iframe>' + "\n"
                mdtext += '</div>'
                mdtext += """

**Figure:** Line plot showing motif distances between motif %s (using highest-scoring hit of motif %s for each input region) and other motifs (from same and different RBPs).
Coverage corresponds to positions occupied by motifs.
The plot is centered on the highest-scoring hit of motif %s (lowest p-value for sequence or highest bit score) for each region containing >= 1 hit of motif %s.
Only motifs with a pair count of >= %i appear in the plot.

&nbsp;

""" %(motif_id, motif_id, motif_id, motif_id, motif_min_pair_count)


            else:

                mdtext += """

<em>No motif distance plot generated for motif %s. Try to lower --motif-min-pair-count (current value: %i).</em>

&nbsp;

""" %(motif_id, motif_min_pair_count)


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

def plot_motif_dist_rbp_level(set_rbp_id,
                              region_rbp_motif_pos_dic,
                              id2name_dic,
                              rbp2mids_dic,
                              reg2pol_dic,
                              html_out="rbp_motif_distances.plotly.html",
                              plotly_js_path="cdn",
                              line_plot_range=30,
                              min_pair_count=1):
    """
    Make motif distances plot on RBP level, relative to given RBP ID set_rbp.

    """

    # Initialize pair counts for each RBP ID, relative to set_rbp.
    set_rbp_other_rbp_pair_count_dic = {}
    rbp_out_range_c_dic = {}
    rbp_in_range_c_dic = {}

    for rbp_id in rbp2mids_dic:
        set_rbp_other_rbp_pair_count_dic[rbp_id] = 0
        rbp_out_range_c_dic[rbp_id] = 0
        rbp_in_range_c_dic[rbp_id] = 0

    # Plot position counts per RBP ID.
    set_rbp_counts_dic = {}
    for rbp_id in rbp2mids_dic:
        set_rbp_counts_dic[rbp_id] = [0]*(2*line_plot_range+1)

    for reg_id in region_rbp_motif_pos_dic:

        # Choose best set_rbp motif in region as center.
        best_motif_id = ""
        best_motif_str = ""
        best_motif_s = 0
        best_motif_e = 0
        best_motif_pval = 1000

        for motif_str in region_rbp_motif_pos_dic[reg_id]:
            motif_id, s, e, p = motif_str.split(",")
            if id2name_dic[motif_id] != set_rbp_id:
                continue
            pval = float(p)
            if pval < best_motif_pval:
                best_motif_id = motif_id
                best_motif_str = motif_str
                best_motif_s = int(s)
                best_motif_e = int(e)
                best_motif_pval = pval

        # If no set_rbp motif hit in region, continue.
        if not best_motif_id:
            continue
        # If set_rbp motif hit is only hit.
        if len(region_rbp_motif_pos_dic[reg_id]) == 1:
            continue

        # Center position of best motif used for centering plot.
        best_motif_cp = get_center_position(best_motif_s-1, best_motif_e)

        # print("best_motif_id:", best_motif_id)
        # print("best_motif_cp:", best_motif_cp)

        # Record counts.
        reg_pol = reg2pol_dic[reg_id]

        paired_rbp_ids_dic = {}

        for motif_str in region_rbp_motif_pos_dic[reg_id]:
            if motif_str == best_motif_str:
                continue
            motif_id, s, e, p = motif_str.split(",")
            rbp_id = id2name_dic[motif_id]
            paired_rbp_ids_dic[rbp_id] = 1

            s = int(s)
            e = int(e)

            # Record for every motif position.
            inside_plot_range = False

            for pos in range(s, e+1):
                dist = pos - best_motif_cp
                # If on minus strand, reverse distance.
                if reg_pol == "-":
                    dist = -1*dist
                if abs(dist) <= line_plot_range:
                    set_rbp_counts_dic[rbp_id][line_plot_range+dist] += 1
                    inside_plot_range = True
                # else:
                #     print("distance %i > line_plot_range")

            if inside_plot_range:
                rbp_in_range_c_dic[rbp_id] += 1
            else:
                rbp_out_range_c_dic[rbp_id] += 1

        for rbp_id in paired_rbp_ids_dic:
            set_rbp_other_rbp_pair_count_dic[rbp_id] += 1


    # Plot top x pair RBPs (pair with set_rbp) only (can also be set_rbp).
    for rbp_id, pair_c in sorted(set_rbp_other_rbp_pair_count_dic.items(), key=lambda item: item[1], reverse=True):
        # print(rbp_id, pair_c)
        if pair_c < min_pair_count:
            # print("Remove %s from plotting ... " %(rbp_id))
            del set_rbp_counts_dic[rbp_id]

    # If not keys remain after min_pair_count filtering.
    if not set_rbp_counts_dic:
        return False, set_rbp_other_rbp_pair_count_dic, rbp_in_range_c_dic, rbp_out_range_c_dic

    line_plot_index = list(range(-line_plot_range, line_plot_range+1))

    # print("set_rbp_counts_dic:", set_rbp_counts_dic)

    df = pd.DataFrame(set_rbp_counts_dic, index=line_plot_index)


    df_reset = df.reset_index()
    df_melt = df_reset.melt(id_vars='index')

    # print("df_melt:", df_melt)
    # index variable  value

    df_melt.columns = ['position', 'rbp_id', 'coverage']

    fig = px.line(df_melt, x='position', y='coverage', color='rbp_id')
    fig.update_layout(
        xaxis_title='RBP ID motif positions relative to RBP ID "%s" motif centers' %(set_rbp_id),
        yaxis_title='Coverage',
        autosize=True,
        legend_title_text='RBP IDs')

    fig.write_html(html_out,
        full_html=False,
        include_plotlyjs=plotly_js_path)

    return True, set_rbp_other_rbp_pair_count_dic, rbp_in_range_c_dic, rbp_out_range_c_dic


################################################################################

def plot_motif_dist_motif_level(set_motif_id,
                              region_rbp_motif_pos_dic,
                              rbp2mids_dic,
                              reg2pol_dic,
                              html_out="motif_distances.plotly.html",
                              plotly_js_path="cdn",
                              line_plot_range=30,
                              min_pair_count=1,
                              mode=1):
    """
    Make motif distances plot on motif level, using set_motif_id motif ID.
    
    mode:
        1: if several set_motif_id motif hits in one region, choose one with 
        best p-value as center.
        2: if several set_motif_id motif hits in one region, use every motif 
        as center once.
    
    """

    # Initialize pair counts for each motif ID, relative best set_rbp motif.
    set_motif_other_motif_pair_count_dic = {}
    motif_out_range_c_dic = {}
    motif_in_range_c_dic = {}

    for rbp_id in rbp2mids_dic:
        for motif_id in rbp2mids_dic[rbp_id]:
            set_motif_other_motif_pair_count_dic[motif_id] = 0
            motif_out_range_c_dic[motif_id] = 0
            motif_in_range_c_dic[motif_id] = 0

    # Plot position counts per motif ID.
    set_motif_counts_dic = {}
    for rbp_id in rbp2mids_dic:
        for motif_id in rbp2mids_dic[rbp_id]:
            set_motif_counts_dic[motif_id] = [0]*(2*line_plot_range+1)


    for reg_id in region_rbp_motif_pos_dic:

        # I case several set_motif_id motif hits in region, select best.
        best_motif_id = ""
        best_motif_str = ""
        best_motif_s = 0
        best_motif_e = 0
        best_motif_pval = 1000

        for motif_str in region_rbp_motif_pos_dic[reg_id]:
            motif_id, s, e, p = motif_str.split(",")
            if motif_id != set_motif_id:
                continue
            pval = float(p)
            if pval < best_motif_pval:
                best_motif_id = motif_id
                best_motif_str = motif_str
                best_motif_s = int(s)
                best_motif_e = int(e)
                best_motif_pval = pval

        # If no set_motif_id motif hit in region, go to next region.
        if not best_motif_id:
            continue
        # If set_motif_id motif hit is only hit.
        if len(region_rbp_motif_pos_dic[reg_id]) == 1:
            continue

        # Center position of best motif used for centering plot.
        best_motif_cp = get_center_position(best_motif_s-1, best_motif_e)

        # print("reg_id:", reg_id)
        # print("best_motif_id:", best_motif_id)
        # print("best_motif_cp:", best_motif_cp)
        # print("best_motif_str:", best_motif_str)

        # Record counts.
        reg_pol = reg2pol_dic[reg_id]
        paired_motif_ids_dic = {}

        for motif_str in region_rbp_motif_pos_dic[reg_id]:
            # Do not compare with motif hit itself.
            if motif_str == best_motif_str:
                continue
            motif_id, s, e, p = motif_str.split(",")
            paired_motif_ids_dic[motif_id] = 1

            s = int(s)
            e = int(e)

            # Record for every motif position.
            inside_plot_range = False

            for pos in range(s, e+1):
                dist = pos - best_motif_cp
                # If on minus strand, reverse distance.
                if reg_pol == "-":
                    dist = -1*dist
                if abs(dist) <= line_plot_range:
                    set_motif_counts_dic[motif_id][line_plot_range+dist] += 1
                    inside_plot_range = True
                # else:
                #     print("distance %i > line_plot_range" %(dist))

            if inside_plot_range:
                motif_in_range_c_dic[motif_id] += 1
            else:
                motif_out_range_c_dic[motif_id] += 1

        for motif_id in paired_motif_ids_dic:
            set_motif_other_motif_pair_count_dic[motif_id] += 1


    # Plot top x pair motifs (pair with set_rbp) only (can also be set_rbp).
    for motif_id, pair_c in sorted(set_motif_other_motif_pair_count_dic.items(), key=lambda item: item[1], reverse=True):
        # print(rbp_id, pair_c)
        if pair_c < min_pair_count:
            # print("Remove %s from plotting (pair_c = %i) ... " %(motif_id, pair_c))
            del set_motif_counts_dic[motif_id]

    # If not keys remain after min_pair_count filtering.
    if not set_motif_counts_dic:
        return False, set_motif_other_motif_pair_count_dic, motif_in_range_c_dic, motif_out_range_c_dic

    line_plot_index = list(range(-line_plot_range, line_plot_range+1))

    # print("set_rbp_counts_dic:", set_rbp_counts_dic)

    df = pd.DataFrame(set_motif_counts_dic, index=line_plot_index)

    df_reset = df.reset_index()
    df_melt = df_reset.melt(id_vars='index')

    # print("df_melt:", df_melt)
    # index variable  value

    df_melt.columns = ['position', 'motif_id', 'coverage']

    fig = px.line(df_melt, x='position', y='coverage', color='motif_id')
    fig.update_layout(
        xaxis_title='Motif positions relative to "%s" motif centers' %(set_motif_id), 
        yaxis_title='Coverage',
        autosize=True,
        legend_title_text='Motif IDs')

    fig.write_html(html_out,
        full_html=False,
        include_plotlyjs=plotly_js_path)

    return True, set_motif_other_motif_pair_count_dic, motif_in_range_c_dic, motif_out_range_c_dic


################################################################################

def create_rbp_reg_occ_upset_plot(rbp2regidx_dic, reg_ids_list,
                                  reg2annot_dic=None,
                                  min_degree=2,
                                  max_degree=None,
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
    c_regions = len(reg_ids_list)
    assert c_regions > 0, "c_regions must be >= 0"
    rbp_idx_list = []  # region indices list [0, 1, 2, ... c_regions-1]
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
    idx_with_hits_dic = {}

    for idx in rbp_idx_list:
        bool_l = []
        for rbp_id, rbp_set in sorted(rbp2regidx_dic.items()):
            if idx in rbp_set:
                idx_with_hits_dic[idx] = 1
                bool_l.append(True)
            else:
                bool_l.append(False)        
        bool_ll.append(bool_l)

    df = pd.DataFrame(bool_ll, columns=rbp_id_list)


    """
    Some checks about degrees and subset sizes.
    """

    df_up = df.groupby(rbp_id_list).size()

    # Check if set min_degree and min_subset_size are too high (to prevent plotting error).
    if len(rbp_id_list) > 1:
        df_up_check = df_up[df_up.index.map(sum) >= min_degree]
        if df_up_check.empty:
            # Set min_degree too high.
            return False, "min_degree", 0
        max_subset_size = 0
        for elem_c in df_up_check:
            if elem_c > max_subset_size:
                max_subset_size = elem_c
        if max_subset_size < min_subset_size:
            # Set min_subset_size too high.
            return False, "min_subset_size", max_subset_size
    else:
        # if min_degree > 1:
        #     return False, "min_degree", 0
        # max_subset_size = len(rbp2regidx_dic[rbp_id_list[0]])
        # if max_subset_size < min_subset_size:
        #     return False, "min_subset_size", max_subset_size
        return False, "len(rbp_id_list) == 1", 0

    """
    Plot with GTF annotations (if reg2annot_dic given),
    which includes adding a stacked bar plot,
    or without.

    https://upsetplot.readthedocs.io/en/latest/auto_examples/plot_discrete.html
    """

    if reg2annot_dic is not None:
        assert len(reg2annot_dic) == c_regions, "len(reg2annot_dic) != c_regions"
        # List of annotations.
        annot_ids_list = []
        annot_with_hits_dic = {}

        for idx, reg_id in enumerate(reg_ids_list):
            annot_id = reg2annot_dic[reg_id][0]
            annot_ids_list.append(annot_id)
            if idx in idx_with_hits_dic:
                annot_with_hits_dic[annot_id] = 1
        df['annot'] = annot_ids_list

        # df = df.sort_values('annot', ascending=False) # Sorting does not change plot.

        # Get annotation ID -> hex color dictionary.
        annot2color_dic = {}
        hex_colors = get_hex_colors_list(min_len=len(annot_with_hits_dic))
        idx = 0
        for annot in sorted(annot_with_hits_dic, reverse=False):
            hc = hex_colors[idx]
            # print("Assigning hex color %s to annotation %s ... " %(hc, annot))
            annot2color_dic[annot] = hex_colors[idx]
            idx += 1
        
        # Set indices (== RBP ID columns).
        df = df.set_index(rbp_id_list)
        # Move on to plotting.
        print("Plotting upset plot with GTF annotations ... ")
        upset = UpSet(df, orientation='horizontal', 
                        min_degree=min_degree,  # number of RBPs in set (e.g. 2 -> at least 2 RBP pairs, not single RBPs)
                        max_degree=max_degree,
                        min_subset_size=min_subset_size,  # min size of a set to be reported.
                        show_counts=True,
                        intersection_plot_elements=0,  # disable the default bar chart.
                        sort_by="cardinality")
        upset.add_stacked_bars(by="annot", colors=annot2color_dic,
                        title="Intersection size", elements=10)  # elements (ticks) of plotting axis it seems.
        upset.plot()
        plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.5)

    else:
        df = df.set_index(rbp_id_list)
        print("Plotting upset plot ... ")
        upset = UpSet(df, orientation='horizontal', 
                      min_degree=min_degree,  # number of RBPs in set (e.g. 2 -> at least 2 RBP pairs, not single RBPs)
                      max_degree=max_degree,
                      min_subset_size=min_subset_size,  # min size of a set to be reported.
                      show_counts=True,
                      sort_by="cardinality")
        upset.plot()

    plt.savefig(plot_out, dpi=125, bbox_inches='tight')
    plt.close()
    return True, "yowza", 0


################################################################################

def create_annotation_stacked_bars_plot(rbp2regidx_dic, reg_ids_list, reg2annot_dic,
                                        plot_out="annotation_stacked_bars_plot.png"):
    """
    Create a stacked bars plot, with each bar showing the annotations for one RBPs,
    i.e. with which GTF annotations the sites with motifs of the RBP overlap.

    rbp2regidx_dic:
        RBP ID -> [indices of regions with RBP motif hits]
    reg_ids_list:
        Region IDs list.
    reg2annot_dic:
        Region ID -> [assigned annotation, transcript ID]

    """
    assert rbp2regidx_dic, "given rbp2regidx_dic empty"
    assert reg_ids_list, "given reg_ids_list empty"
    assert reg2annot_dic, "given reg2annot_dic empty"

    # Scale plot height depending on # of features.
    c_ids = len(rbp2regidx_dic)
    fheight = 1 * c_ids
    fwidth = 10

    # Get all annotation IDs in dataset.
    annot_dic = {}
    for reg_id in reg2annot_dic:
        annot = reg2annot_dic[reg_id][0]
        if annot not in annot_dic:
            annot_dic[annot] = 1
        else:
            annot_dic[annot] += 1

    rbp_id_list = []
    rbp2idx_dic = {}
    idx = 0
    for rbp_id, rbp_list in sorted(rbp2regidx_dic.items(), reverse=True):
        rbp_id_list.append(rbp_id)
        rbp2idx_dic[rbp_id] = idx
        idx += 1

    data_dic = {}
    for annot in sorted(annot_dic, reverse=True):  # make reverse sort to have same color coding as upset plot.
        data_dic[annot] = []
        for rbp_id in rbp_id_list:
            data_dic[annot].append(0)

    data_dic["rbp_id"] = []
    for rbp_id in rbp_id_list:
        data_dic["rbp_id"].append(rbp_id)

    annot_with_hits_dic = {}

    for rbp_id in rbp2regidx_dic:
        for hit_idx in rbp2regidx_dic[rbp_id]:
            reg_id = reg_ids_list[hit_idx]
            annot = reg2annot_dic[reg_id][0]
            annot_with_hits_dic[annot] = 1
            rbp_idx = rbp2idx_dic[rbp_id]
            data_dic[annot][rbp_idx] += 1

    df = pd.DataFrame(data_dic)
    # Remove annotation columns with no counts.
    df = df.loc[:, (df != 0).any(axis=0)]

    # Get annotation ID -> hex color mapping for plotting.
    annot2color_dic = {}

    hex_colors = get_hex_colors_list(min_len=len(annot_with_hits_dic))

    idx = 0
    for annot in sorted(annot_with_hits_dic, reverse=False):
        hc = hex_colors[idx]
        # print("Assigning hex color %s to annotation %s ... " %(hc, annot))
        annot2color_dic[annot] = hex_colors[idx]
        idx += 1

    ax = df.set_index('rbp_id').plot(kind='barh', stacked=True, legend=False, color=annot2color_dic, edgecolor="none", figsize=(fwidth, fheight))
    # ax = df.set_index('rbp_id').plot(kind='barh', stacked=True, legend=False, edgecolor="lightgrey", figsize=(fwidth, fheight))

    plt.xlabel('Annotation overlap')
    ax.set_ylabel('')
    ax.yaxis.grid(False)
    ax.xaxis.grid(True)
    ax.set_axisbelow(True)

    # Remove border lines.
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)

    # plt.legend(loc=(1.01, 0.4), fontsize=12, framealpha=0)
    # plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.5)
    plt.legend(bbox_to_anchor=(1.01, 1), loc='upper left')

    plt.savefig(plot_out, dpi=125, bbox_inches='tight')
    plt.close()


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
                    'hovertemplate': 'RBP1: %{x}<br>RBP2: %{y}<br>p-value: %{customdata[0]}<br>RBPs: %{customdata[1]}<br>Counts: %{customdata[2]}<br>Correlation: %{customdata[3]}<br>-log10(p-value): %{z}<extra></extra>'}])
    plot.update_layout(plot_bgcolor='white')
    plot.write_html(plot_out,
                    full_html=False,
                    include_plotlyjs=plotly_js_path)


################################################################################

def boundary_check(value, min, max):
    """
    Return if value is within integer/float boundaries.

    >>> value = 0.5
    >>> min = 1E-9
    >>> max = 1.0
    >>> boundary_check(value, min, max)
    True
    >>> value = 0.0
    >>> min = 1E-9
    >>> max = 1.0
    >>> boundary_check(value, min, max)
    False

    """
    if value >= min and value <= max:
        return True
    else:
        return False


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
                                     motif_db_str=False,
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

    # Delete plots if already present.
    if os.path.exists(plots_out_folder):
        shutil.rmtree(plots_out_folder)
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

            # Check if motif in motif database folder.
            if motif_db_str:
                db_motif_path = benchlib_path + "/content/%s_motif_plots/%s" %(motif_db_str, motif_plot)
                if os.path.exists(db_motif_path):
                    shutil.copy(db_motif_path, motif_plot_out)

            if not os.path.exists(motif_plot_out):
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
