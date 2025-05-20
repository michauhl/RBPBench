from distutils.spawn import find_executable
from typing import Optional
import os
import re
import subprocess
import gzip
import shutil
import statistics
from itertools import combinations, product
from collections import defaultdict
import matplotlib.pyplot as plt
from matplotlib_venn import venn2, venn3
from venn import venn
from logomaker import Logo
from markdown import markdown
import pandas as pd
import plotly.express as px
from math import floor, log10, ceil, log
import textdistance
import numpy as np
from upsetplot import UpSet
from packaging import version
from sklearn.decomposition import PCA
from sklearn.preprocessing import MinMaxScaler
# from sklearn.decomposition import SparsePCA
from itertools import product
from scipy.stats import fisher_exact
from scipy.stats import false_discovery_control  # Benjamini-Hochberg correction.
from scipy.stats import mannwhitneyu
from goatools.obo_parser import GODag
from goatools.goea.go_enrichment_ns import GOEnrichmentStudy
from decimal import Decimal, getcontext

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

def bed_read_in_regions(in_bed,
                        no_id_check=False,
                        use_regions_as_ids=False):
    """
    Read in BED regions into a dictionary.

    """

    bed_regions_dic = {}

    with open(in_bed) as f:
        for line in f:
            cols = line.strip().split("\t")
            chr_id = cols[0]
            start = int(cols[1])
            end = int(cols[2])
            reg_id = cols[3]
            pol = cols[5]
            
            if use_regions_as_ids:
                reg_id = chr_id + ":" + str(start) + "-" + str(end) + "," + pol

            if no_id_check:
                if reg_id in bed_regions_dic:
                    # Overwrite existing region.
                    print("WARNING: region ID \"%s\" already found in \"%s\". Overwriting existing region .." % (reg_id, in_bed))
            else:
                assert reg_id not in bed_regions_dic, "region ID \"%s\" already found in \"%s\". Please provide BED file with unique col4 IDs or set --use-regions" % (reg_id, in_bed)
            bed_regions_dic[reg_id] = [chr_id, start, end, pol]

    f.closed
    return bed_regions_dic


################################################################################

def get_region_dic_stats(region_dic):
    """
    region_dic format:
    {region_id: [chr_id, start, end, strand]} with 0-based start.

    """
    assert region_dic, "region_dic empty"
    # Get length stats.
    region_lengths =  []
    for reg_id, reg_info in region_dic.items():
        reg_len = reg_info[2] - reg_info[1]
        region_lengths.append(reg_len)
    # Get mean, median, min, max.
    mean_len = statistics.mean(region_lengths)
    median_len = statistics.median(region_lengths)
    min_len = min(region_lengths)
    max_len = max(region_lengths)
    # Get number of regions.
    num_regions = len(region_lengths)
    return num_regions, mean_len, median_len, min_len, max_len


################################################################################

def get_val_dic_stats(val_dic):
    """
    val_dic format:
    {data_id: score_value}

    """
    assert val_dic, "val_dic empty"
    val_list = []
    for data_id, val in val_dic.items():
        val_list.append(val)
    # Get mean, median, min, max.
    mean_val = statistics.mean(val_list)
    median_val = statistics.median(val_list)
    min_val = min(val_list)
    max_val = max(val_list)
    num_val = len(val_list)
    return num_val, mean_val, median_val, min_val, max_val


################################################################################

def create_con_sc_violin_plot(in_scores, control_scores, plot_out, 
                              pval=1.0,
                              con_type="phastCons",
                              in_id="Input sites",
                              control_id="Control sites",
                              add_pval=True,
                              plot_title=False,
                              plot_pdf=False):
    """
    Create conservation scores comparison plot.
    
    """

    plt.figure(figsize=(6, 4))
    plt.violinplot([in_scores, control_scores], showmeans=True, widths=0.8)
    plt.xticks([1, 2], [in_id, control_id])
    plt.ylabel("Mean site %s conservation score" %(con_type))
    if plot_title:
        plt.title("%s Conservation Score Distribution" %(con_type))

    # Optional: Add p-value to plot.
    if add_pval:
        plt.text(1.5, max(max(in_scores), max(control_scores)), f'p = {pval:.3g}', 
                ha='center', va='bottom', fontsize=10)

    plt.tight_layout()
    # plt.show()
    plt.savefig(plot_out, dpi=250)

    if plot_pdf and plot_out.endswith('.png'):
        pdf_out = plot_out[:-4] + '.pdf'
        # plt.savefig(pdf_out, dpi=110, bbox_inches='tight')
        plt.savefig(pdf_out, dpi=110)


################################################################################

def con_generate_html_report(args, stats_dic, benchlib_path,
                             html_report_out="report.rbpbench_con.html",
                             pc_plot_name="phastCons_scores.png",
                             pp_plot_name="phyloP_scores.png",
                             pc_pval=1.0,
                             pp_pval=1.0):

    """
    Create rbpbench con report.

    """

    # Check if plots are there.
    pc_plot = False
    pp_plot = False
    pc_plot_path_check = os.path.join(args.out_folder, pc_plot_name)
    pp_plot_path_check = os.path.join(args.out_folder, pp_plot_name)
    pc_plot_path = pc_plot_name
    pp_plot_path = pp_plot_name
    if os.path.exists(pc_plot_path_check):
        pc_plot = True
    if os.path.exists(pp_plot_path_check):
        pp_plot = True

    if not pp_plot and not pc_plot:
        assert False, "Neither %s nor %s plot found in output folder" %(pc_plot_name, pp_plot_name)

    assert "in_regions_stats" in stats_dic, "in_regions_stats not in stats_dic"
    assert "ctrl_regions_stats" in stats_dic, "ctrl_regions_stats not in stats_dic"

    out_folder = args.out_folder
    # Use absolute paths?
    if args.plot_abs_paths:
        out_folder = os.path.abspath(out_folder)
        pc_plot_path = os.path.join(out_folder, pc_plot_name)
        pp_plot_path = os.path.join(out_folder, pp_plot_name)

    # Version string.
    version_str = "v" + args.version

    html_out = os.path.join(out_folder, html_report_out)
    md_out = os.path.join(out_folder, "report.rbpbench_con.md")

    """
    Setup sorttable.js to make tables in HTML sortable.

    """
    sorttable_js_path = benchlib_path + "/content/sorttable.js"
    assert os.path.exists(sorttable_js_path), "sorttable.js not at %s" %(sorttable_js_path)
    sorttable_js_html = '<script src="' + sorttable_js_path + '" type="text/javascript"></script>'
    if args.sort_js_mode == 2:
        shutil.copy(sorttable_js_path, out_folder)
        sorttable_js_html = '<script src="sorttable.js" type="text/javascript"></script>'  # path in HTML code.
    elif args.sort_js_mode == 3:
        js_code = read_file_content_into_str_var(sorttable_js_path)
        sorttable_js_html = "<script>\n" + js_code + "\n</script>\n"

    # Logo path.
    logo_path_html = os.path.join(out_folder, "logo.png")
    if not os.path.exists(logo_path_html):
        logo_path = benchlib_path + "/content/logo.png"
        assert os.path.exists(logo_path), "logo.png not found in %s" %(logo_path)
        shutil.copy(logo_path, out_folder)

    # HTML head section.
    html_head = """<!DOCTYPE html>
<html>
<head>
<title>RBPBench - Conservation Scores Comparison Report</title>

<style>
    th, td {
        border: 1px solid black;
        padding: 8px;
        text-align: center;
    }
    th {
        background-color: #f2f2f2;
    }
    .page {
        page-break-after: always;
    }
    .title-container {
        display: flex;
        align-items: center;
    }
    .title-container img {
        margin-right: 10px; /* Adjust the spacing as needed */
    }
</style>

</head>

<div class="title-container">
    <img src="logo.png" alt="Logo" width="175">
    <h1>Conservation Scores Comparison Report</h1>
</div>

<body>
"""

    # HTML tail section.
    html_tail = """
%s
</body>
</html>
""" %(sorttable_js_html)

    # Markdown part.
    mdtext = """

List of available statistics and plots generated
by RBPBench (%s, rbpbench con):

- [Site length statistics](#site-len-stats)
""" %(version_str)

    if pc_plot:
        mdtext += "- [PhastCons conservation scores comparison](#pc-comp)\n"
    if pp_plot:
        mdtext += "- [PhyloP conservation scores comparison](#pp-comp)\n"

    mdtext += "\n"
    mdtext += "\n&nbsp;\n"


    """
    Site lengths statistics table.

    """

    mdtext += """
## Site length statistics ### {#site-len-stats}

**Table:** Site length statistics for provided BED region files (input and control sites).

"""

    mdtext += '<table style="max-width: 750px; width: 100%; border-collapse: collapse; line-height: 0.8;">' + "\n"
    mdtext += "<thead>\n"
    mdtext += "<tr>\n"
    mdtext += "<th>Dataset ID</th>\n"
    mdtext += "<th># sites</th>\n"
    mdtext += "<th>Mean length</th>\n"
    mdtext += "<th>Median length</th>\n"
    mdtext += "<th>Min length</th>\n"
    mdtext += "<th>Max length</th>\n"
    mdtext += "</tr>\n"
    mdtext += "</thead>\n"
    mdtext += "<tbody>\n"

    for dataset_id, dataset_stats in stats_dic.items():
        if dataset_id == "in_regions_stats" or dataset_id == "ctrl_regions_stats":
            data_id = "Input sites"
            if dataset_id == "ctrl_regions_stats":
                data_id = "Control sites"
            mdtext += "<tr>\n"
            mdtext += "<td>%s</td>\n" %(data_id)
            mdtext += "<td>%i</td>\n" %(dataset_stats[0])
            mdtext += "<td>%.2f</td>\n" %(dataset_stats[1])
            mdtext += "<td>%.2f</td>\n" %(dataset_stats[2])
            mdtext += "<td>%i</td>\n" %(dataset_stats[3])
            mdtext += "<td>%i</td>\n" %(dataset_stats[4])
            mdtext += "</tr>\n"

    mdtext += "</tbody>\n"
    mdtext += "</table>\n"

    mdtext += "\n&nbsp;\n&nbsp;\n"
    mdtext += "\nColumn IDs have the following meanings: "
    mdtext += "**Dataset ID** -> dataset ID (input or control), sites provided via --in (input) and --ctrl-in (control), "
    mdtext += '**# sites** -> number of genomic sites in dataset, '
    mdtext += "**Mean length** -> mean site length in dataset, "
    mdtext += "**Median length** -> median site length in dataset, "
    mdtext += "**Min length** -> minimum site length in dataset, "
    mdtext += "**Min length** -> minimum site length in dataset, "
    mdtext += "**Max length** -> maximum site length in dataset." + "\n"
    mdtext += "\n&nbsp;\n"


    """"
    PhastCons conservation scores comparison plot.

    """

    if pc_plot:

        assert "in_phastcons_stats" in stats_dic, "in_phastcons_stats not in stats_dic"
        assert "ctrl_phastcons_stats" in stats_dic, "ctrl_phastcons_stats not in stats_dic"

        mdtext += """
## PhastCons conservation scores comparison ### {#pc-comp}

Distribution of phastCons conservation scores in input and control sites, taking the average conservation
score for each site (i.e., averaged over all site positions).

"""
        # Stats.
        in_c_sites = stats_dic["in_phastcons_stats"][0]
        ctrl_c_sites = stats_dic["ctrl_phastcons_stats"][0]
        in_mean = stats_dic["in_phastcons_stats"][1]
        ctrl_mean = stats_dic["ctrl_phastcons_stats"][1]
        in_median = stats_dic["in_phastcons_stats"][2]
        ctrl_median = stats_dic["ctrl_phastcons_stats"][2]
        in_min = stats_dic["in_phastcons_stats"][3]
        ctrl_min = stats_dic["ctrl_phastcons_stats"][3]
        in_max = stats_dic["in_phastcons_stats"][4]
        ctrl_max = stats_dic["ctrl_phastcons_stats"][4]

        # pc_pval
        pval_info = "Small test p-values (< 0.05) indicate that input sites have significantly higher conservation scores than control sites."
        if args.wrs_mode == 2:
            pval_info = "Small test p-values (< 0.05) indicate that input sites have significantly lower conservation scores than control sites."

        mdtext += '<img src="' + pc_plot_path + '" alt="phastCons scores comparison plot"' + "\n"
        mdtext += 'title="phastCons scores comparison plot" width="900" />' + "\n"
        mdtext += """
**Figure:** Distribution of phastCons conservation scores in input and control sites.
For each site, the average phastCons score is used (i.e., average over all genomic site positions).
Wilcoxon rank-sum test p-value = %s.
Wilcoxon rank-sum test is applied to check for significant differences between input and control site scores.
%s
 # input sites = %i, # control sites = %i, mean input score = %.2f, mean control score = %.2f,
 # median input score = %.2f, median control score = %.2f,
 # min input score = %.2f, min control score = %.2f,
 # max input score = %.2f, max control score = %.2f.

&nbsp;

""" %(str(pc_pval), pval_info, in_c_sites, ctrl_c_sites, in_mean, ctrl_mean, in_median, ctrl_median, in_min, ctrl_min, in_max, ctrl_max)



    """"
    PhyloP conservation scores comparison plot.

    """

    if pp_plot:

        assert "in_phylop_stats" in stats_dic, "in_phylop_stats not in stats_dic"
        assert "ctrl_phylop_stats" in stats_dic, "ctrl_phylop_stats not in stats_dic"

        mdtext += """
## PhyloP conservation scores comparison ### {#pp-comp}

Distribution of phyloP conservation scores in input and control sites, taking the average conservation
score for each site (i.e., averaged over all site positions).

"""
        # Stats.
        in_c_sites = stats_dic["in_phylop_stats"][0]
        ctrl_c_sites = stats_dic["ctrl_phylop_stats"][0]
        in_mean = stats_dic["in_phylop_stats"][1]
        ctrl_mean = stats_dic["ctrl_phylop_stats"][1]
        in_median = stats_dic["in_phylop_stats"][2]
        ctrl_median = stats_dic["ctrl_phylop_stats"][2]
        in_min = stats_dic["in_phylop_stats"][3]
        ctrl_min = stats_dic["ctrl_phylop_stats"][3]
        in_max = stats_dic["in_phylop_stats"][4]
        ctrl_max = stats_dic["ctrl_phylop_stats"][4]

        # pc_pval
        pval_info = "Small test p-values (< 0.05) indicate that input sites have significantly higher conservation scores than control sites."
        if args.wrs_mode == 2:
            pval_info = "Small test p-values (< 0.05) indicate that input sites have significantly lower conservation scores than control sites."

        mdtext += '<img src="' + pp_plot_path + '" alt="phyloP scores comparison plot"' + "\n"
        mdtext += 'title="phyloP scores comparison plot" width="900" />' + "\n"
        mdtext += """
**Figure:** Distribution of phyloP conservation scores in input and control sites.
For each site, the average phyloP score is used (i.e., average over all genomic site positions).
Wilcoxon rank-sum test p-value = %s.
Wilcoxon rank-sum test is applied to check for significant differences between input and control site scores.
%s
 # input sites = %i, # control sites = %i, mean input score = %.2f, mean control score = %.2f,
 # median input score = %.2f, median control score = %.2f,
 # min input score = %.2f, min control score = %.2f,
 # max input score = %.2f, max control score = %.2f.

&nbsp;

""" %(str(pp_pval), pval_info, in_c_sites, ctrl_c_sites, in_mean, ctrl_mean, in_median, ctrl_median, in_min, ctrl_min, in_max, ctrl_max)


    # Convert mdtext to html.
    md2html = markdown(mdtext, extensions=['attr_list', 'tables'])

    # OUTMD = open(md_out,"w")
    # OUTMD.write("%s\n" %(mdtext))
    # OUTMD.close()

    html_content = html_head + md2html + html_tail
    OUTHTML = open(html_out,"w")
    OUTHTML.write("%s\n" %(html_content))
    OUTHTML.close()


################################################################################

def compare_conservation_scores(args,
                                in_regions_dic, control_regions_dic, 
                                benchlib_path,
                                pc_bw=False,
                                pp_bw=False,
                                html_report_out="report.rbpbench_con.html",
                                in_con_sc_name="in_regions.avg_con_sc.tsv",
                                ctrl_con_sc_name="control_regions.avg_con_sc.tsv",
                                pc_plot_name="phastCons_scores.png",
                                pp_plot_name="phyloP_scores.png"):
    """
    Extract and compare phylogenomic conservation scores for given genomic input 
    and control sites.

    in_regions_dic / control_regions_dic format:
    {region_id: [chr_id, start, end, strand]} with 0-based start.

    """

    assert in_regions_dic, "in_regions_dic empty"
    assert control_regions_dic, "control_regions_dic empty"
    assert pc_bw or pp_bw, "No conservation scores provided. Please provide --phastcons and/or --phylop files"
    assert os.path.exists(args.out_folder), "given --out folder \"%s\" not found" % (args.out_folder)
    if pc_bw:
        assert os.path.exists(pc_bw), "given --phastcons file \"%s\" not found" % (pc_bw)
    if pp_bw:
        assert os.path.exists(pp_bw), "given --phylop file \"%s\" not found" % (pp_bw)

    import pyBigWig

    # Get region length stats.
    in_regions_stats = get_region_dic_stats(in_regions_dic)
    ctrl_regions_stats = get_region_dic_stats(control_regions_dic)

    in_reg_avg_phastcons_dic = {}
    in_reg_avg_phylop_dic = {}
    ctrl_reg_avg_phastcons_dic = {}
    ctrl_reg_avg_phylop_dic = {}

    pc_plot_path = os.path.join(args.out_folder, pc_plot_name)
    pp_plot_path = os.path.join(args.out_folder, pp_plot_name)
    html_report_path = os.path.join(args.out_folder, html_report_out)

    # Remove old plots and HTML report.
    if os.path.exists(pc_plot_path):
        os.remove(pc_plot_path)
    if os.path.exists(pp_plot_path):
        os.remove(pp_plot_path)
    if os.path.exists(html_report_path):
        os.remove(html_report_path)

    # Wilcoxon rank-sum test / Mann Whitney U test mode.
    wrs_alt_hypo = "greater"
    if args.wrs_mode == 1:
        wrs_alt_hypo = "greater"
        print("Check if --in sites have significantly higher conservation scores ... ")
    elif args.wrs_mode == 2:
        wrs_alt_hypo = "less"
        print("Check if --in sites have significantly lower conservation scores ... ")
    else:
        assert False, "Invalid Wilcoxon rank-sum (Mann Whitney U) test mode: %i" %(args.wrs_mode)


    pc_pval = 1.0
    pp_pval = 1.0

    """
    PhastCons conservation scores.

    """

    if pc_bw:

        print("Read in phastCons conservation scores ... ")

        pc_bw_data = pyBigWig.open(pc_bw)

        for reg_id, reg_info in in_regions_dic.items():

            chr_id = reg_info[0]
            start = reg_info[1]
            end = reg_info[2]

            try:
                # Get conservation scores for the region.
                scores = pc_bw_data.values(chr_id, start, end, numpy=False)
                # Convert NaN values to 0.0.
                scores = [0.0 if np.isnan(s) else s for s in scores]
                avg_score = statistics.mean(scores) if scores else 0.0
                in_reg_avg_phastcons_dic[reg_id] = avg_score
            except RuntimeError:
                print(f"Skipping --in site {chr_id}:{start}-{end} (coordinates not in bigWig)")

        for reg_id, reg_info in control_regions_dic.items():

            chr_id = reg_info[0]
            start = reg_info[1]
            end = reg_info[2]

            try:
                # Get conservation scores for the region.
                scores = pc_bw_data.values(chr_id, start, end, numpy=False)
                # Convert NaN values to 0.0.
                scores = [0.0 if np.isnan(s) else s for s in scores]
                avg_score = statistics.mean(scores) if scores else 0.0
                ctrl_reg_avg_phastcons_dic[reg_id] = avg_score
            except RuntimeError:
                print(f"Skipping --ctrl-in site {chr_id}:{start}-{end} (coordinates not in bigWig)")

        pc_bw_data.close()

        # Make violin plot.
        print("Create plot ... ")

        in_scores = []
        for reg_id, reg_info in in_regions_dic.items():
            in_scores.append(in_reg_avg_phastcons_dic[reg_id])
        control_scores = []
        for reg_id, reg_info in control_regions_dic.items():
            control_scores.append(ctrl_reg_avg_phastcons_dic[reg_id])

        # Compare distributions.
        stat, pc_pval = mannwhitneyu(in_scores, control_scores, alternative=wrs_alt_hypo)
        # Round p-value to 4 significant digits.
        pc_pval = round_to_n_significant_digits_v2(pc_pval, 4,
                                                   min_val=0.0)

        print(f"Wilcoxon rank-sum test: U = {stat:.2f}, p = {pc_pval:.4g}")

        create_con_sc_violin_plot(in_scores, control_scores, pc_plot_path, 
                                  pval=pc_pval,
                                  con_type="phastCons",
                                  add_pval=True,
                                  plot_title=False,
                                  plot_pdf=args.plot_pdf)

    """
    PhyloP conservation scores.

    """

    if pp_bw is not None:

        print("Read in phyloP conservation scores ... ")

        pp_bw_data = pyBigWig.open(pp_bw)

        for reg_id, reg_info in in_regions_dic.items():

            chr_id = reg_info[0]
            start = reg_info[1]
            end = reg_info[2]

            try:
                # Get conservation scores for the region.
                scores = pp_bw_data.values(chr_id, start, end, numpy=False)
                # Convert NaN values to 0.0.
                scores = [0.0 if np.isnan(s) else s for s in scores]
                avg_score = statistics.mean(scores) if scores else 0.0
                in_reg_avg_phylop_dic[reg_id] = avg_score
            except RuntimeError:
                print(f"Skipping --in site {chr_id}:{start}-{end} (coordinates not in bigWig)")

        for reg_id, reg_info in control_regions_dic.items():

            chr_id = reg_info[0]
            start = reg_info[1]
            end = reg_info[2]

            try:
                # Get conservation scores for the region.
                scores = pp_bw_data.values(chr_id, start, end, numpy=False)
                # Convert NaN values to 0.0.
                scores = [0.0 if np.isnan(s) else s for s in scores]
                avg_score = statistics.mean(scores) if scores else 0.0
                ctrl_reg_avg_phylop_dic[reg_id] = avg_score
            except RuntimeError:
                print(f"Skipping --ctrl-in site {chr_id}:{start}-{end} (coordinates not in bigWig)")

        pp_bw_data.close()

        # Make violin plot.
        print("Create plot ... ")

        in_scores = []
        for reg_id, reg_info in in_regions_dic.items():
            in_scores.append(in_reg_avg_phylop_dic[reg_id])
        control_scores = []
        for reg_id, reg_info in control_regions_dic.items():
            control_scores.append(ctrl_reg_avg_phylop_dic[reg_id])

        # Compare distributions.
        stat, pp_pval = mannwhitneyu(in_scores, control_scores, alternative=wrs_alt_hypo)
        # Round p-value to 4 significant digits.
        pp_pval = round_to_n_significant_digits_v2(pp_pval, 4,
                                                   min_val=0.0)
        
        print(f"Wilcoxon rank-sum test: U = {stat:.2f}, p = {pp_pval:.4g}")

        create_con_sc_violin_plot(in_scores, control_scores, pp_plot_path, 
                                  pval=pp_pval,
                                  con_type="phyloP",
                                  add_pval=True,
                                  plot_title=False,
                                  plot_pdf=args.plot_pdf)

    stats_dic = {}
    stats_dic["in_regions_stats"] = in_regions_stats
    stats_dic["ctrl_regions_stats"] = ctrl_regions_stats

    if pc_bw:
        in_phastcons_stats = get_val_dic_stats(in_reg_avg_phastcons_dic)
        ctrl_phastcons_stats = get_val_dic_stats(ctrl_reg_avg_phastcons_dic)
        stats_dic["in_phastcons_stats"] = in_phastcons_stats
        stats_dic["ctrl_phastcons_stats"] = ctrl_phastcons_stats

    if pp_bw:
        in_phylop_stats = get_val_dic_stats(in_reg_avg_phylop_dic)
        ctrl_phylop_stats = get_val_dic_stats(ctrl_reg_avg_phylop_dic)
        stats_dic["in_phylop_stats"] = in_phylop_stats
        stats_dic["ctrl_phylop_stats"] = ctrl_phylop_stats


    """
    Create report.

    """

    print("Create HTML report ... ")

    con_generate_html_report(args, stats_dic, benchlib_path,
                             html_report_out=html_report_out,
                             pc_plot_name=pc_plot_name,
                             pp_plot_name=pp_plot_name,
                             pc_pval=pc_pval,
                             pp_pval=pp_pval)

    """
    Output site conservation scores.

    """

    print("Output site conservation scores ... ")

    in_reg_con_sc_out = os.path.join(args.out_folder, in_con_sc_name)
    ctrl_reg_con_sc_out = os.path.join(args.out_folder, ctrl_con_sc_name)

    SCOUT = open(in_reg_con_sc_out, "w")
    SCOUT.write("site_id\tchr_id\tsite_s\tsite_e\tphastcons_sc\tphylop_sc\n")

    for reg_id, reg_info in in_regions_dic.items():
        chr_id = reg_info[0]
        start = reg_info[1]
        end = reg_info[2]
        avg_phastcons = in_reg_avg_phastcons_dic.get(reg_id, "-")
        avg_phylop = in_reg_avg_phylop_dic.get(reg_id, "-")
        avg_phastcons = str(avg_phastcons)
        avg_phylop = str(avg_phylop)
        row_str = "%s\t%s\t%s\t%s\t%s\t%s\n" % (reg_id, chr_id, start, end, avg_phastcons, avg_phylop)
        SCOUT.write("%s" % (row_str))

    SCOUT.close()

    SCOUT = open(ctrl_reg_con_sc_out, "w")
    SCOUT.write("site_id\tchr_id\tsite_s\tsite_e\tphastcons_sc\tphylop_sc\n")

    for reg_id, reg_info in control_regions_dic.items():
        chr_id = reg_info[0]
        start = reg_info[1]
        end = reg_info[2]
        avg_phastcons = ctrl_reg_avg_phastcons_dic.get(reg_id, "-")
        avg_phylop = ctrl_reg_avg_phylop_dic.get(reg_id, "-")
        avg_phastcons = str(avg_phastcons)
        avg_phylop = str(avg_phylop)
        row_str = "%s\t%s\t%s\t%s\t%s\t%s\n" % (reg_id, chr_id, start, end, avg_phastcons, avg_phylop)
        SCOUT.write("%s" % (row_str))

    SCOUT.close()


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
    # Extract last output line, in case of warnings printed.
    version_string = output.split('\n')[-1]
    version_string = version_string.strip()
    tool_version = version.parse(version_string)
    min_version = version.parse(min_version)
    check = tool_version >= min_version
    return check, tool_version


################################################################################

def is_valid_regex(regex):
    """
    Check if regex string is valid regex.

    >>> is_valid_regex(".*")
    True
    >>> is_valid_regex(".*[")
    False
    >>> is_valid_regex("ACGT")
    True

    """

    try:
        re.compile(regex)
        return True
    except re.error:
        return False


################################################################################

def get_seqs_from_regex(regex):
    """
    Given a regex, e.g., A[CG]T[AG], return all possible sequences 
    that match the regex.

    >>> get_seqs_from_regex("A[CG]T[AG]")
    ['ACTA', 'ACTG', 'AGTA', 'AGTG']
    
    """

    assert regex, "given regex empty"
    assert is_valid_regex(regex), "Regex \"%s\" is not a valid regex" %(regex)

    # Extract fixed parts and variable parts.
    fixed_parts = re.split(r'\[.*?\]', regex)
    variable_parts = re.findall(r'\[(.*?)\]', regex)

    # Generate all combinations of the variable parts.
    combinations = list(product(*variable_parts))

    # Assemble all sequences matched by the regex.
    seqs_list = []
    for combo in combinations:
        seq = ''.join([fixed_part + combo_part for fixed_part, combo_part in zip(fixed_parts, combo + ('',))])
        seqs_list.append(seq)

    assert seqs_list, "no matching sequences extracted from regex \"%s\"" %(regex)

    return seqs_list


################################################################################

def get_seq_parts_from_regex(regex):
    """
    Split a regex into single sequence parts, e.g.
    'A[CG]AT[AG][CT]C' -> ['A', 'CG', 'AT', 'AG', 'CT', 'C']

    >>> get_seq_parts_from_regex('A[CG]AT[AG]AA')
    ['A', 'CG', 'A', 'T', 'AG', 'A', 'A']
    >>> get_seq_parts_from_regex('ACGT')
    ['A', 'C', 'G', 'T']
    >>> get_seq_parts_from_regex('[ACGT]')
    ['ACGT']
    
    """
    assert regex, "given regex empty"
    assert is_valid_regex(regex), "Regex \"%s\" is not a valid regex" %(regex)

    search_regex = r'[A-Z]+|\[[A-Z]+\]'

    parts = re.findall(search_regex, regex)

    split_parts = []
    for part in parts:
        if part.startswith('[') and part.endswith(']'):
            split_parts.append(part)
        else:
            split_parts.extend(list(part))

    # Remove brackets from variable parts.
    final_parts = [part.strip('[]') for part in split_parts]

    assert final_parts, "no sequence parts extracted from regex \"%s\"" %(regex)

    return final_parts


################################################################################

def seq_check_alphabet(seq, 
                       alphabet=["A", "C", "G", "T"]):
    """
    Check if sequence uses only characters from given alphabet.

    >>> seq_check_alphabet("ACGT", alphabet=["A", "C", "G", "T"])
    True
    >>> seq_check_alphabet("ACGT", alphabet=["A", "C", "G"])
    False

    """

    assert seq, "given seq empty"

    for nt in seq:
        if nt not in alphabet:
            return False
    
    return True


################################################################################

def seq_parts_to_motif_block(seq_parts_list,
                             alphabet=["A", "C", "G", "T"]):
    """
    Convert sequence parts lists with format:
    ['A', 'CG', 'A', 'T', 'AG', 'A', 'A']
    to motif block format.

    >>> seq_parts_to_motif_block(["A"], alphabet=["A", "C", "G", "T"])
    ['letter-probability matrix: alength= 4 w= 1 nsites= 20 E= 0', '1.00000  0.00000  0.00000  0.00000']
    >>> seq_parts_to_motif_block(["ACG"], alphabet=["A", "C", "G", "T"])
    ['letter-probability matrix: alength= 4 w= 1 nsites= 20 E= 0', '0.33333  0.33333  0.33333  0.00000']
    
    """

    assert seq_parts_list, "given seq_parts_list empty"

    l_ab = len(alphabet)
    matrix_head = "letter-probability matrix: alength= %i w= %i nsites= 20 E= 0" %(l_ab, len(seq_parts_list))
    seq_motif_block = [matrix_head]
    for part in seq_parts_list:
        # part can be single nt or multiple nts string.
        nts = set(list(part))
        single_prob = 1.0 / len(nts)
        line = ""
        for i in range(l_ab):
            if alphabet[i] in nts:
                # Make single_prob 5 decimal places format, so always print 5 decimal places.
                line += f"{single_prob:.5f}" + "  "
            else:
                line += "0.00000  "

        seq_motif_block.append(line.strip())

    return seq_motif_block


################################################################################

def seq_to_motif_block(seq,
                       alphabet=["A", "C", "G", "T"]):
    """
    Convert sequence to motif block format.
    
    >>> seq_to_motif_block("A", alphabet=["A", "C", "G", "T"])
    ['letter-probability matrix: alength= 4 w= 1 nsites= 20 E= 0', '1.00000  0.00000  0.00000  0.00000']

    """

    assert seq, "given seq empty"

    l_seq = len(seq)
    l_ab = len(alphabet)
    matrix_head = "letter-probability matrix: alength= %i w= %i nsites= 20 E= 0" %(l_ab, l_seq)
    seq_motif_block = [matrix_head]
    for nt in seq:
        line = ""
        for i in range(l_ab):
            if alphabet[i] == nt:
                line += "1.00000  "
            else:
                line += "0.00000  "

        seq_motif_block.append(line.strip())

    return seq_motif_block

    """
    seq_motif_block:
    ['letter-probability matrix: alength= 4 w= 6 nsites= 20 E= 0', ' 0.054545  0.636364  0.145455  0.163636 ', ' 0.814815  0.055555  0.000000  0.129630 ', ' 0.168381  0.000000  0.831619  0.000000 ', ' 0.163636  0.072727  0.181818  0.581819 ', ' 0.200000  0.127273  0.636363  0.036364 ', ' 0.290909  0.400000  0.090909  0.218182 ']

    MOTIF XRN2_1 
    letter-probability matrix: alength= 4 w= 9 nsites= 20 E= 0
    0.034700  0.046200  0.855500  0.063600 
    0.041100  0.057800  0.862900  0.038200 
    0.051200  0.048900  0.822700  0.077200 
    0.069107  0.053105  0.827683  0.050105 
    0.026900  0.130200  0.812900  0.030000 
    0.025900  0.527700  0.350300  0.096100 
    0.033500  0.823800  0.099300  0.043400 
    0.083192  0.575342  0.241976  0.099490 
    0.058200  0.577800  0.297300  0.066700 
    """


################################################################################

def get_godag_obo_file_online(godag_obo_file,
                              godag_obo_url="https://purl.obolibrary.org/obo/go/go-basic.obo"):
    """
    Try to download godag obo file via given godag_obo_url. Store file in path:
    godag_obo_file

    """

    import requests
    from requests.exceptions import RequestException

    # Check if file is available via given URL.
    try:
        response = requests.head(godag_obo_url)
        # if response.status_code == 200:
        #     print("go-basic.obo file is available.")
        # else:
        #     print("go-basic.obo File is not available.")
    except RequestException as e:
        print("Internet connection for downloading go-basic.obo file not available")
        print(f"Error details: {e}")

    # Download the file.
    try:
        response = requests.get(godag_obo_url)
        with open(godag_obo_file, 'wb') as f:
            f.write(response.content)
    except RequestException as e:
        print(f"Error downloading go-basic.obo file: {e}")


################################################################################

def get_godag_obo_file_local(godag_obo_gz_file, godag_obo_file):
    """
    Extract godag obo file from godag obo gz file.

    """

    with gzip.open(godag_obo_gz_file, 'rb') as f_in:
        with open(godag_obo_file, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)


################################################################################

def get_gid2go_mapping(gid2go_file,
                       go_dag=False,
                       remove_version_numbers=False,
                       list2set=True):
    """
    Given an (Ensembl) gene ID to GO ID mapping file, read the file and return a 
    dictionary with gene IDs as keys and GO ID list as values.
    Gene ID to GO mapping file format:
    gene_id<tab>go_id1,go_id2,go_id3,...
    gid2go_file format:
    {gene_id: [go_id1, go_id2, ...]}

    list2set:
        If True, convert GO ID lists to sets.
    go_dag:
        If go_dag (get via: go_dag = GODag(godag_obo_file)) provided, only 
        keep GO IDs that are present in the GO DAG.

    """

    # Gene ID -> GO IDs list mapping.
    gid2go_dic = {}
    with gzip.open(gid2go_file, 'rt') as f:
        for line in f:
            if line.startswith("gene_id"):
                continue
            columns = line.strip().split('\t')
            gene_id = columns[0]
            if remove_version_numbers:
                gene_id = re.sub(r"\.\d+$", "", gene_id)

            go_ids = columns[1].split(",")
            filtered_go_ids = []
            for go_id in go_ids:
                if go_dag:
                    if go_id in go_dag:
                        filtered_go_ids.append(go_id)
                else:
                    filtered_go_ids.append(go_id)
            if not filtered_go_ids:
                continue
            if list2set:
                filtered_go_ids = set(filtered_go_ids)
            gid2go_dic[gene_id] = filtered_go_ids

    assert gid2go_dic, "gid2go_dic empty. No GO terms stored for any gene ID. Make sure %s file has correct row format: gene_id<tab>go_id1,go_id2,..." %(gid2go_file)

    return gid2go_dic


################################################################################

class GeneDesc:
    """
    Stores gene description (compact gene infos objects).

    """

    def __init__(self,
                 gene_id: str,
                 gene_name: str,
                 gene_synonyms: str,
                 gene_type: str,
                 gene_desc: str) -> None:
        self.gene_id = gene_id
        self.gene_name = gene_name
        self.gene_synonyms = gene_synonyms
        self.gene_type = gene_type
        self.gene_desc = gene_desc


################################################################################

def get_gene_descriptions(target_reg_annot_file):
    """
    Read in gene descriptions.

    E.g. get from ensembl_gene_infos.biomart.GRCh38.112.tsv.gz:
    gene_id	gene_name	gene_synonyms	gene_type	gene_description
    ENSG00000286112	KYAT1	CCBL1;GTK;KAT1;KATI	protein_coding	kynurenine aminotransferase 1 [Source:NCBI gene (formerly Entrezgene);Acc:883]
    ENSG00000171097	KYAT1	CCBL1;GTK;KAT1;KATI;CCBL1;GTK;KATI	protein_coding	kynurenine aminotransferase 1 [Source:HGNC Symbol;Acc:HGNC:1564]
    ENSG00000291087	CRYBB2P1	CRYB2B	lncRNA	crystallin beta B2 pseudogene 1 [Source:NCBI gene (formerly Entrezgene);Acc:1416]
    ...
    
    """

    gene_desc_dic = {}
    with gzip.open(target_reg_annot_file, 'rt') as f:
        for line in f:
            if line.startswith("gene_id"):
                continue
            cols = line.strip().split("\t")
            gene_id = cols[0]
            gene_name = cols[1]
            gene_synonyms = cols[2]
            gene_type = cols[3]
            gene_desc = cols[4]
            gene_desc_obj = GeneDesc(gene_id, gene_name, gene_synonyms, gene_type, gene_desc)
            gene_desc_dic[gene_id] = gene_desc_obj

    assert gene_desc_dic, "no gene descriptions read in from %s" %(target_reg_annot_file)
    return gene_desc_dic


################################################################################

def output_target_reg_annot(target_genes_dic, gene_infos_file, target_reg_annot_file,
                            remove_version_numbers=True,
                            gid2tid_dic=None,
                            tid2tio_dic=None):

    """
    Output target region annotations (target regions used in GOA).

    Note that if internal gene_infos_file is outdated compared to given GOA/GTF files, 
    this could result in missing gene descriptions for some genes.
                              
    """

    assert target_genes_dic, "No target genes provided"
    # check if gene_infos_file file exists.
    assert os.path.exists(gene_infos_file), "gene infos file %s not found" %(gene_infos_file)

    gene_desc_dic = get_gene_descriptions(gene_infos_file)

    OUTANNOT = open(target_reg_annot_file, "w")

    OUTANNOT.write("gene_id\tgene_name\tgene_synonyms\tgene_type\tgene_description\ttr_id\ttr_type\n")

    for gene_id in target_genes_dic:

        gene_id_full = gene_id

        if remove_version_numbers:
            gene_id = re.sub(r"\.\d+$", "", gene_id)

        gene_name = "-"
        gene_synonyms = "-"
        gene_type = "-"
        gene_desc = "-"
        if gene_id in gene_desc_dic:
            gene_name = gene_desc_dic[gene_id].gene_name
            gene_synonyms = gene_desc_dic[gene_id].gene_synonyms
            gene_type = gene_desc_dic[gene_id].gene_type
            gene_desc = gene_desc_dic[gene_id].gene_desc

        tr_id = "-"
        tr_type = "-"
        if gid2tid_dic is not None:
            if gene_id in gid2tid_dic or gene_id_full in gid2tid_dic:
                tr_id = gid2tid_dic[gene_id]
                if tid2tio_dic is not None:
                    if tr_id in tid2tio_dic:
                        tr_type = tid2tio_dic[tr_id].tr_biotype

        OUTANNOT.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %(gene_id, gene_name, gene_synonyms, gene_type, gene_desc, tr_id, tr_type))

    OUTANNOT.close()


################################################################################

def run_go_analysis(target_genes_dic, background_genes_dic, 
                    gid2go_file, out_folder,
                    pval_thr=0.05,
                    excluded_terms = ["GO:0005575", "GO:0008150", "GO:0003674"],
                    goa_obo_mode=1,
                    stats_dic=None,
                    propagate_counts=True,
                    store_gene_names=True,
                    goa_obo_file=False):

    """
    Run GO enrichment analysis for target genes, using goatools:
    conda install goatools
    
    pval_thr:
        BH corrected p-value threshold.
    goa_obo_mode:
        1: Get GO DAG file online.
        2: Get GO DAG file locally.
        3: Use provided GO DAG file.
    store_gene_names:
        If True, use values of background_genes_dic (which are their gene names),
        instead of gene IDs in dataframe.
    
    """

    assert target_genes_dic, "No target genes provided"
    assert background_genes_dic, "No background genes provided"

    godag_obo_file = out_folder + "/go-basic.obo"
    if stats_dic is None:
        stats_dic = {}

    print("Run GO enrichment analysis ...")

    # Get GO DAG file.
    if goa_obo_mode == 1:

        print("Get GO DAG obo file online ...")
        get_godag_obo_file_online(godag_obo_file)
        assert os.path.exists(godag_obo_file), "Downloaded GO DAG file not found in output folder. Please check internet connection or use local file via --goa_obo_mode 2"

    elif goa_obo_mode == 2:

        print("Use local GO DAG obo file ...")
        get_godag_obo_file_local(goa_obo_file, godag_obo_file)
        assert os.path.exists(godag_obo_file), "Local GO DAG file not found in output folder. Please contact developers"

    elif goa_obo_mode == 3:

        assert goa_obo_file, "--goa-obo-mode 3, but no GO DAG obo file provided. Please provide --goa-obo-file or change --goa-obo-mode"
        print("Use GO DAG obo file provided via --goa-obo-file ...")
        # Copy goa_obo_file to godag_obo_file.
        shutil.copyfile(goa_obo_file, godag_obo_file)

    # Read in GO DAG.
    print("Load GO DAG file ...")
    go_dag = GODag(godag_obo_file)

    # Check if gene IDs have version numbers.
    id_has_version = False
    version_pattern = re.compile(r"\.\d+$")
    for gene_id in target_genes_dic:
        # if re.search("\.\d+$", gene_id):
        if version_pattern.search(gene_id):
            id_has_version = True
        break

    # Read in gene ID 2 GO ID(s) mapping.
    print("Load gene ID to GO ID mapping ...")
    gid2go_dic = get_gid2go_mapping(gid2go_file,
                                    remove_version_numbers=id_has_version,
                                    go_dag=go_dag,  # Filter out GO terms not present in GO DAG.
                                    list2set=True)  # Convert GO ID lists to sets.

    print("# genes in gene ID to GO ID mapping:", len(gid2go_dic))

    stats_dic["c_genes_from_gid2go"] = len(gid2go_dic)

    print("# of target genes (before GO ID filtering):    ", len(target_genes_dic))
    print("# of background genes (before GO ID filtering):", len(background_genes_dic))

    # Remove version numbers for GOA.
    print("Filter to keep only genes with associated GO IDs ...")
    target_genes = []
    for gene_id in target_genes_dic:
        if id_has_version:
            gene_id = re.sub(r"\.\d+$", "", gene_id)
        if gene_id in gid2go_dic:  # Can only work with genes that have associated GO terms.
            target_genes.append(gene_id)
    background_genes = []
    gid2gn_dic = {}
    for gene_id in background_genes_dic:
        new_gene_id = gene_id
        if id_has_version:
            new_gene_id = re.sub(r"\.\d+$", "", gene_id)
        gid2gn_dic[new_gene_id] = background_genes_dic[gene_id]
        if new_gene_id in gid2go_dic:
            background_genes.append(new_gene_id)

    assert target_genes, "No target genes left after filtering for GO terms. Make sure to provide --gtf with compatible gene IDs (internal IDs are ENSEMBL style gene IDs). Alternatively, provide gene ID to GO ID mapping file via --goa-gene2go-file"
    assert background_genes, "No background genes left after filtering for GO terms. Make sure to provide --gtf with compatible gene IDs (internal IDs are ENSEMBL style gene IDs). Alternatively, provide gene ID to GO ID mapping file via --goa-gene2go-file"

    print("# of target genes with GO IDs (used for GOA):    ", len(target_genes))
    print("# of background genes with GO IDS (used for GOA):", len(background_genes))

    stats_dic["c_target_genes_goa"] = len(target_genes)
    stats_dic["c_background_genes_goa"] = len(background_genes)

    # Namespace short to long name mapping.
    ns2name_dic = {
        'BP': 'biological_process',
        'MF': 'molecular_function',
        'CC': 'cellular_component'
    }

    # Create GOEnrichmentStudy object.
    goeaobj = GOEnrichmentStudy(
        background_genes,   # List of background genes.
        gid2go_dic,         # Gene ID -> GO IDs mapping.
        go_dag,             # GO directed acyclic graph.
        propagate_counts=propagate_counts,
        alpha=0.05,
        methods=['fdr_bh']  # Multiple testing correction method.
    )

    # Run GO enrichment analysis.
    goea_results_all = goeaobj.run_study(target_genes)

    results_dic = {
        "GO": [],
        "term": [],
        "class": [],
        "p": [],
        "p_corr": [],
        "enrichment": [],
        "depth": [],
        "n_children": [],
        "n_genes": [],
        "n_study": [],
        "perc_genes": [],
        "study_genes": []
    }

    stats_dic["c_sig_go_terms_e"] = 0
    stats_dic["c_sig_go_terms_p"] = 0

    # Loop over significant results.
    for res in goea_results_all:

        go_id = res.GO

        if res.p_fdr_bh > pval_thr:
            continue

        if go_id in excluded_terms:
            continue

        go_n_children = len(res.goterm.get_all_children()) 
        go_name = res.name
        namespace = res.NS
        class_name = ns2name_dic[namespace]
        p_uncorrected = res.p_uncorrected
        p_corrected = res.p_fdr_bh
        study_count = res.study_count  # number of genes for this result.
        study_items = res.study_items  # Gene names for this result.
        study_n = res.study_n  # number target genes.
        enrichment = res.enrichment  # e for enriched, p for purified.
        go_depth = res.depth

        # if filter_purified and enrichment == "p":
        #     continue

        perc_in_study = 0.0
        if res.ratio_in_study[1] > 0:
            # Rounded to 2 decimal places percentage.
            perc_in_study = round((res.ratio_in_study[0] / res.ratio_in_study[1]) * 100, 2)

        p_uncorrected = round_to_n_significant_digits_v2(p_uncorrected, 4)
        p_corrected = round_to_n_significant_digits_v2(p_corrected, 4)

        gene_list = []
        for gene_id in study_items:
            gene_name = gene_id
            if store_gene_names:
                if gene_id in gid2gn_dic:
                    gene_name = gid2gn_dic[gene_id]
            gene_list.append(gene_name)

        gene_list_str = ""
        if len(study_items) > 0:
            gene_list_str = ",".join(gene_list)
        
        results_dic["GO"].append(go_id)
        results_dic["term"].append(go_name)
        results_dic["class"].append(class_name)
        results_dic["p"].append(p_uncorrected)
        results_dic["p_corr"].append(p_corrected)
        results_dic["enrichment"].append(enrichment)
        results_dic["depth"].append(go_depth)
        results_dic["n_children"].append(go_n_children)
        results_dic["n_genes"].append(study_count)
        results_dic["n_study"].append(study_n)
        results_dic["perc_genes"].append(perc_in_study)
        results_dic["study_genes"].append(gene_list_str)
    
        if enrichment == "e":
            stats_dic["c_sig_go_terms_e"] += 1
        elif enrichment == "p":
            stats_dic["c_sig_go_terms_p"] += 1
        
    # Make dataframe out of results dictionary.
    goa_results_df = pd.DataFrame(results_dic)

    c_goa_results = len(goa_results_df)

    print("# of significant GO terms:", c_goa_results)
    stats_dic["c_sig_go_terms"] = c_goa_results
    if excluded_terms:
        excluded_terms_str = ",".join(excluded_terms)
        stats_dic["excluded_terms"] = excluded_terms_str

    return goa_results_df


################################################################################

def round_to_n_significant_digits_v2(num, n, zero_check_val=1e-304,
                                     min_val=0):
    """
    Round float / scientific notation number to n significant digits.
    This function only works for positive numbers.

    >>> round_to_n_significant_digits_v2(4.0980000000000007e-38, 4)
    4.098e-38
    >>> round_to_n_significant_digits_v2(4.410999999999999e-81, 4)
    4.411e-81
    >>> round_to_n_significant_digits_v2(0.0000934234823499234, 4)
    9.342e-05
    >>> round_to_n_significant_digits_v2(0.0112123123123123, 4)
    0.01121
    >>> round_to_n_significant_digits_v2(0.0, 2)
    0
    >>> round_to_n_significant_digits_v2(1e-300, 3)
    1e-300

    """

    getcontext().prec = n  # Set precision to n.

    if num < zero_check_val:
        return min_val
    else:
        d_num = Decimal(num)  # Convert float to decimal.
        return float(round(d_num, -int(floor(log10(abs(d_num))) - (n - 1))))


################################################################################

def round_to_n_significant_digits(num, n,
                                  zero_check_val=1e-300):
    """
    Round float / scientific notation number to n significant digits.
    This function only works for positive numbers.
    
    >>> round_to_n_significant_digits(0.296538210, 3)
    0.297
    >>> round_to_n_significant_digits(0.0002964982, 3)
    0.000296
    >>> round_to_n_significant_digits(0.0000296498, 3)
    2.96e-05
    >>> round_to_n_significant_digits(0.0000296498, 2)
    3e-05
    >>> round_to_n_significant_digits(0.0, 2)
    0
    >>> round_to_n_significant_digits(1e-300, 3)
    1e-300

    
    """

    if num < zero_check_val:
        return 0
    else:
        return round(num, -int(floor(log10(abs(num))) - (n - 1)))


################################################################################

def read_in_tomtom_tsv(tomtom_tsv):
    """
    Read in tomtom.tsv file (produced by tomtom call).

    """

    pair2sim_dic = {}
    motif_ids_dic = {}

    line_c = 0
    with open(tomtom_tsv) as f:
        for line in f:
            line_c += 1
            if line_c == 1:
                continue
            cols = line.strip().split("\t")
            if len(cols) != 10:
                continue
            q_id = cols[0]
            t_id = cols[1]
            opt_offset = int(cols[2])
            pval = float(cols[3])
            e_val = float(cols[4])
            q_val = float(cols[5])
            overlap = int(cols[6])
            q_cons = cols[7]
            t_cons = cols[8]
            orient = cols[9]

            pval = round_to_n_significant_digits_v2(pval, 4,
                                                    min_val=1e-304)

            log_pval = log_tf_pval(pval)
            log_pval = round_to_n_significant_digits_v2(log_pval, 4,
                                                        min_val=0)

            motif_ids_dic[q_id] = 1

            pair_str1 = q_id + "," + t_id
            pair_str2 = t_id + "," + q_id
            pair2sim_dic[pair_str1] = log_pval
            pair2sim_dic[pair_str2] = log_pval

    f.closed

    return pair2sim_dic, motif_ids_dic


################################################################################

def output_tomtom_sim_results(motif_ids_dic, pair2sim_dic, sim_out_tsv,
                              header=False):
    """
    Output tomtom similarity results.

    Similarity == -log10(pval_tomtom)
    pval_tomtom measuring how similar two motifs are (the lower the more similar).
    
    """

    motif_ids_list = [x for x in motif_ids_dic.keys()]
    motif_ids_list.sort()

    motif_pairs = list(combinations(motif_ids_list, 2))

    OUTSIM = open(sim_out_tsv, "w")
    if header:
        OUTSIM.write("motif1\tmotif2\tsimilarity\n")

    for pair in motif_pairs:
        pair = list(pair)
        pair.sort()
        pair_str = pair[0] + "," + pair[1]

        assert pair_str in pair2sim_dic, "motif pair string \"%s\" not in pair2sim_dic" %(pair_str)

        sim = pair2sim_dic[pair_str]

        OUTSIM.write("%s\t%s\t%s\n" %(pair[0], pair[1], str(sim)))

    # Also get motif ID with itself pairs.
    for motif_id in motif_ids_list:
        
        pair_str = motif_id + "," + motif_id

        assert pair_str in pair2sim_dic, "motif pair string \"%s\" not in pair2sim_dic" %(pair_str)

        sim = pair2sim_dic[pair_str]

        OUTSIM.write("%s\t%s\t%s\n" %(motif_id, motif_id, str(sim)))

    OUTSIM.close()


################################################################################

def read_in_tomtom_sim_results(sim_out_tsv,
                               motif_sim_cap=50):
    """
    Read in tomtom similarity results.

    """

    assert os.path.exists(sim_out_tsv), "sim_out_tsv %s not found" %(sim_out_tsv)

    pair2sim_dic = {}

    if re.search(r".+\.gz$", sim_out_tsv):
        f = gzip.open(sim_out_tsv, 'rt')
    else:
        f = open(sim_out_tsv, "r")
    
    for line in f:
        cols = line.strip().split("\t")
        if cols[2] == "similarity":
            continue
        motif1 = cols[0]
        motif2 = cols[1]
        sim = float(cols[2])

        if motif_sim_cap:
            if sim > motif_sim_cap:
                sim = motif_sim_cap

        pair_str1 = motif1 + "," + motif2
        pair_str2 = motif2 + "," + motif1
        pair2sim_dic[pair_str1] = sim
        pair2sim_dic[pair_str2] = sim

    f.closed

    assert pair2sim_dic, "pair2sim_dic empty (no similarity values read in from %s)" %(sim_out_tsv)

    return pair2sim_dic


################################################################################

def calc_tomtom_sim(seq_motifs_db_file, out_folder):
    """
    Based on given seq_motifs_db_file (MEME motif format file), 
    calculate similarities between all motifs in file.
    Return motif pair to similarity score dictionary, format:
    "motif_id1,motif_id2" -> similarity_score 


    tomtom -norc -oc mdb3_sim_out -dist ed -thresh 1.0 
    """
    assert os.path.exists(seq_motifs_db_file), "seq_motifs_db_file not found"

    # Bang the tomtom slowly, dumbass.
    run_tomtom(seq_motifs_db_file, seq_motifs_db_file, out_folder,
                tomtom_bfile=False,
                tomtom_thresh=1.0,
                tomtom_evalue=False,
                tomtom_m=False,
                tomtom_min_overlap=1,
                params="-norc -dist ed",
                call_dic=None,
                print_output=False,
                error_check=False)

    tomtom_tsv = out_folder + "/tomtom.tsv"
    pair2sim_dic, motif_ids_dic = read_in_tomtom_tsv(tomtom_tsv)

    tomtom_sim_out = out_folder + "/tomtom_sim.tsv"

    output_tomtom_sim_results(motif_ids_dic, pair2sim_dic, tomtom_sim_out, header=True)

    motif_pair2sim_dic = read_in_tomtom_sim_results(tomtom_sim_out)

    return motif_pair2sim_dic


################################################################################

def search_regex_in_seqs_dic(regex, seqs_dic,
                             step_size_one=False,
                             case_sensitive=True):
    """
    Get regex matches in sequences stored in seqs_dic.
    match.start : 0-based index of the start of the match
    match.end : 1-based index of the end of the match
    match.group : the matched string
    
    >>> seqs_dic = {'seq1': 'XXXXXXX'}
    >>> regex = "AA"
    >>> search_regex_in_seqs_dic(regex, seqs_dic)
    {}
    >>> seqs_dic = {'seq1': 'ATXXAGX'}
    >>> regex = "A[GT]"
    >>> search_regex_in_seqs_dic(regex, seqs_dic)
    {'seq1': [[0, 2, 'AT'], [4, 6, 'AG']]}
    >>> regex = "a[GT]"
    >>> seqs_dic = {'seq1': 'ATXXAGX'}
    >>> search_regex_in_seqs_dic(regex, seqs_dic, case_sensitive=False)
    {'seq1': [[0, 2, 'AT'], [4, 6, 'AG']]}
    >>> regex = "AAAC"
    >>> seqs_dic = {'seq1': 'AAA'}
    >>> search_regex_in_seqs_dic(regex, seqs_dic)
    {}
    >>> regex = "AA"
    >>> seqs_dic = {'seq1': 'AAAA', 'seq2': 'TTTT'}
    >>> search_regex_in_seqs_dic(regex, seqs_dic, step_size_one=True)
    {'seq1': [[0, 2, 'AA'], [1, 3, 'AA'], [2, 4, 'AA']]}

    """

    if not is_valid_regex(regex):
        return "Invalid regex"
    
    hits_dic = {}
    seq_c = 0

    # Pre-compile the regular expression.
    flags = 0
    if not case_sensitive:
        flags |= re.IGNORECASE
    if all(ord(char) < 128 for char in regex):
        flags |= re.ASCII
    compiled_regex = re.compile(regex, flags)

    if step_size_one:

        for seq_name, seq in seqs_dic.items():
            seq_c += 1
            if seq_c % 1000 == 0:
                print(f"{seq_c} sequences scanned ... ")
            
            seq_length = len(seq)
            for i in range(seq_length):
                for match in compiled_regex.finditer(seq[i:]):
                    if match.start() != 0:
                        break
                    if seq_name not in hits_dic:
                        hits_dic[seq_name] = [[i + match.start(), i + match.end(), match.group()]]
                    else:
                        hits_dic[seq_name].append([i + match.start(), i + match.end(), match.group()])

    else:

        for seq_name, seq in seqs_dic.items():
            seq_c += 1
            if seq_c % 1000 == 0:
                print(f"{seq_c} sequences scanned ... ")
            
            for match in compiled_regex.finditer(seq):
                if seq_name not in hits_dic:
                    hits_dic[seq_name] = [[match.start(), match.end(), match.group()]]
                else:
                    hits_dic[seq_name].append([match.start(), match.end(), match.group()])

    return hits_dic


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

def check_table_file(in_file):
    """
    Check if file is a 4 column tab-separated table file, with last column containing 
    valid file paths.

    Table file containing dataset infos, with columns:
    rbp_id method_id data_id path_to_BED_file
    PUM1 	clipper_rep1 	k562_eclip 	path/to/PUM1.k562_eclip.clipper_rep1.bed
    PUM1 	clipper_rep2 	k562_eclip 	path/to/PUM1.k562_eclip.clipper_rep2.bed
    PUM1 	clipper_idr 	k562_eclip 	path/to/PUM1.k562_eclip.clipper_idr.bed
    PUM1 	dewseq_w100_s5 	k562_eclip 	path/to/PUM1.k562_eclip.dewseq_w100_s5.bed
    RBFOX2 	clipper_idr 	hepg2_eclip 	path/to/RBFOX2.hepg2_eclip.clipper_idr.bed
    RBFOX2 	clipper_idr 	k562_eclip 	path/to/RBFOX2.k562_eclip.clipper_idr.bed 

    """

    is_table = True

    with open(in_file) as f:
        for line in f:
            cols = line.strip().split("\t")
            if len(cols) != 4:
                is_table = False
            if not os.path.exists(cols[3]):
                is_table = False
            break
    f.closed

    return is_table


################################################################################

def read_in_table_file(in_file):
    """
    Read in table file, containing dataset infos, with tab-separated columns:
    rbp_id method_id data_id path_to_BED_file

    """

    dataset_list = []
    with open(in_file) as f:
        for line in f:
            cols = line.strip().split("\t")
            rbp_id = cols[0]
            method_id = cols[1]
            data_id = cols[2]
            path = cols[3]
            assert os.path.exists(path), "BED file path %s in --bed table file not found. Please provide valid file paths inside table file" %(path)
            dataset_list.append([rbp_id, method_id, data_id, path])
    f.closed
    assert dataset_list, "no dataset infos read in from --bed table file %s" %(in_file)

    return dataset_list


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
            if re.search(r".+\.%s" %(file_ending), df):
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

def closest_dist_motif_pos_lists(l1, l2):
    """
    Given two lists of motif center positions, calculate closest distance 
    between any two positions in the two lists.

    >>> l1 = [30, 90]
    >>> l2 = [50, 100]
    >>> closest_dist_motif_pos_lists(l1, l2)
    10
    >>> l1 = []
    >>> print(closest_dist_motif_pos_lists(l1, l2))
    None

    """

    min_dist = None
    for i in l1:
        for j in l2:
            dist = abs(i-j)
            if min_dist is None:
                min_dist = dist
            else:
                if dist < min_dist:
                    min_dist = dist
    return min_dist


################################################################################

def run_k_nt_shuffling(seqs_fa, out_fa,
                       kmer_size=2,
                       params="-dna",
                       error_check=True,
                       tag="_shuf",
                       seed=None):
    """
    
    Do k-nucleotide (kmer_size-nt) shuffling of FASTA sequences in seqs_fa. 
    Output shuffled sequences to out_fa.

    https://meme-suite.org/meme/doc/fasta-shuffle-letters.html

    fasta-shuffle-letters [options] <sequence file> [<output file>]

    -kmer 2
    -dna
    -tag	The name of the sequence will have text appended to it.	
            default: The name of the sequence will have "_shuf" appended to it.
    -seed
            
    """

    assert is_tool("fasta-shuffle-letters"), "fasta-shuffle-letters not in PATH. MEME suite (v5) not installed?"

    params += " -kmer %s " %(str(kmer_size))

    if seed is not None:
        params += " -seed %s " %(str(seed))

    params += ' -tag "%s" ' %(tag)

    check_cmd = "fasta-shuffle-letters " + params + " " + seqs_fa + " " + out_fa
    output = subprocess.getoutput(check_cmd)

    if error_check:
        error = False
        if output:
            error = True
        assert error == False, "fasta-shuffle-letters is complaining:\n%s\n%s" %(check_cmd, output)


################################################################################

def make_contingency_table_2x2_v2(region_labels_dic, idx1, idx2,
                                  rid2rbpidx2hcp_dic,
                                  max_motif_dist=50):
    """
    Make a contingency table 2x2, using region_labels_dic.
    region_labels_dic format:
    region_id -> [False, True, False ... ]
    with list number of RBP IDs (len_rbp_list), alphabetically sorted.
    True: region covered by RBP at idx (i.e. motif hits)
    False: region not covered by RBP at idx (i.e. no motif hits)
    Also use hit center position(s) list from rid2rbpidx2hcp_dic to calculate
    average distance of closest pairs of motif hits for both RBPs.
    rid2rbpidx2hcp_dic format:
    region_id -> rbp_idx -> hit center position(s) list

    Return table format:
    table = [[A, B], [C, D]]

    https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.fisher_exact.html
                   List 1              Not in List 1
    List 2         A                   B
    Not in List 2  C                   D

    table = [[A, B], [C, D]]

    >>> region_labels_dic = {'r1': [True, False], 'r2': [False, True], 'r3': [True, True], 'r4': [False, False], 'r5': [True, True]}
    >>> rid2rbpidx2hcp_dic = {'r1': {0: [30]}, 'r2': {1: [100]}, 'r3': {0: [30], 1: [50, 70]}, 'r4': {}, 'r5': {0: [30, 40], 1: [20, 100]}}
    >>> make_contingency_table_2x2_v2(region_labels_dic, 0, 1, rid2rbpidx2hcp_dic, max_motif_dist=50)
    ([[2, 1], [1, 1]], '15.0', '100.0')
    >>> make_contingency_table_2x2_v2(region_labels_dic, 0, 1, rid2rbpidx2hcp_dic, max_motif_dist=10)
    ([[2, 1], [1, 1]], '15.0', '50.0')
    >>> make_contingency_table_2x2_v2(region_labels_dic, 0, 1, rid2rbpidx2hcp_dic, max_motif_dist=5)
    ([[2, 1], [1, 1]], '15.0', '0.0')
    
    """

    cont_a = 0
    cont_b = 0
    cont_c = 0
    cont_d = 0

    sum_min_dist = 0
    c_cont_a_max_motif_dist = 0

    for reg_id in region_labels_dic:
        in1 = region_labels_dic[reg_id][idx1]
        in2 = region_labels_dic[reg_id][idx2]

        if in1 and in2:
            cont_a += 1
            min_dist = closest_dist_motif_pos_lists(rid2rbpidx2hcp_dic[reg_id][idx1], rid2rbpidx2hcp_dic[reg_id][idx2])
            if min_dist is not None:
                # print(reg_id, min_dist)
                sum_min_dist += min_dist
                if min_dist <= max_motif_dist:
                    c_cont_a_max_motif_dist += 1

        elif not in1 and in2:
            cont_b += 1
        elif in1 and not in2:
            cont_c += 1
        else: # not in1 and not in2
            cont_d += 1

    avg_min_dist = "-"
    if cont_a:
        avg_min_dist = sum_min_dist / cont_a
        # Round avg_min_dist to 1 decimal.
        avg_min_dist = str(round(avg_min_dist, 1))

    # Percentage of regions with motif hits for both RBPs, where closest motif hits are <= max_motif_dist away.
    perc_close_hits = "-"
    if cont_a:
        perc_close_hits = 100 * c_cont_a_max_motif_dist / cont_a
        # Round perc_close_hits to 2 decimals.
        perc_close_hits = str(round(perc_close_hits, 2))

    table = [[cont_a, cont_b], [cont_c, cont_d]]
    return table, avg_min_dist, perc_close_hits


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
            if re.search(r"^INFERNAL", line):
                blocks_list.append(line)
                idx += 1
            elif re.search(r"^ACC\s+\w+", line):
                m = re.search(r"^ACC\s+(\w+)", line)
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
            # if re.search("^#", line):
            if line.startswith("#"):
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
            if re.search(r"^ACC\s+\w+", line):
                m = re.search(r"^ACC\s+(\w+)", line)
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
        if re.search(r"^>", line):
            if full_header:
                m = re.search(r"^>(.+)", line)
                seq_id = m.group(1)
                seq_ids_dic[seq_id] = 1
            else:
                m = re.search(r"^>(\S+)", line)
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

def run_streme(in_fa, out_folder,
               neg_fa=False,
               streme_bfile=False,
               streme_evalue=False,
               streme_thresh=0.05,
               streme_minw=6,
               streme_maxw=15,
               streme_seed=0,
               streme_order=2,
               params="--dna",
               call_dic=None,
               print_output=True,
               error_check=False):
    """
    Run STREME on input FASTA file in_fa.
    
    """
    assert is_tool("streme"), "streme not in PATH"
    assert os.path.exists(in_fa), "in_fa \"%s\" does not exist" %(in_fa)

    if neg_fa:
        params += " --n " + neg_fa
    if streme_bfile:
        params += " --bfile " + streme_bfile
    if streme_evalue:
        params += " --evalue"
    params += " --thresh " + str(streme_thresh)
    params += " --minw " + str(streme_minw)
    params += " --maxw " + str(streme_maxw)
    params += " --seed " + str(streme_seed)
    params += " --order " + str(streme_order)
    params += " --oc " + out_folder

    check_cmd = "streme " + params + " -p " + in_fa

    output = subprocess.getoutput(check_cmd)

    if call_dic is not None:
        call_dic["streme_call"] = check_cmd

    if error_check:
        error = False
        if output:
            error = True
        assert error == False, "streme is complaining:\n%s\n%s" %(check_cmd, output)
    if print_output:
        if output:
            print("")
            print("STREME COMMAND:\n%s" %(check_cmd))
            print("STREME OUTPUT:\n%s" %(output))
            print("")


################################################################################

def run_tomtom(query_meme, db_meme, out_folder,
               tomtom_bfile=False,
               tomtom_thresh=0.5,
               tomtom_evalue=False,
               tomtom_m=False,
               tomtom_min_overlap=1,
               params="-norc",
               call_dic=None,
               print_output=True,
               error_check=False):

    """
    Run TOMTOM using query motifs file query_meme and database motifs file db_meme.

    Output folder content:
    tomtom.tsv
    tomtom.xml
    tomtom.html

    -dist ed (default) not compatible with a set -bfile. Using 'allr' instead.

    """
    assert is_tool("tomtom"), "tomtom not in PATH"
    assert os.path.exists(query_meme), "query_meme \"%s\" does not exist" %(query_meme)
    assert os.path.exists(db_meme), "db_meme \"%s\" does not exist" %(db_meme)
    
    if tomtom_bfile:
        params += " -bfile %s" %(tomtom_bfile)
        params += " -dist allr"
    params += " -thresh %s" %(str(tomtom_thresh))
    if tomtom_evalue:
        params += " -evalue"
    if tomtom_m:
        params += " -m"
    params += " -min-overlap %s" %(str(tomtom_min_overlap))
    params += " -oc %s" %(out_folder)

    check_cmd = "tomtom " + params + " " + query_meme + " " + db_meme
    output = subprocess.getoutput(check_cmd)

    if call_dic is not None:
        call_dic["tomtom_call"] = check_cmd

    if error_check:
        error = False
        if output:
            error = True
        assert error == False, "tomtom is complaining:\n%s\n%s" %(check_cmd, output)
    if print_output:
        if output:
            print("")
            print("TOMTOM COMMAND:\n%s" %(check_cmd))
            print("TOMTOM OUTPUT:\n%s" %(output))
            print("")


################################################################################

def bed_extract_sequences_from_fasta(in_bed, in_fa, out_fa,
                                     add_param="",
                                     print_warnings=False,
                                     ignore_errors=False):
    """
    Extract sequences from genome (provide genome .fasta file).
    Gzipped FASTA not supported by bedtools getfasta so far ...

    print_warnings:
        Instead of raising an error when output is encountered,
        print output (i.e. bedtools getfasta warnings)

    Q: what happens if BED file contains regions not fully included in FASTA?
    This should in region sequence not being extracted at all (even if partially overlapping 
    with FASTA), plus throwing error:
    Feature (t1:19-25) beyond the length of t1 size (20 bp).  Skipping.

        
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

    check_cmd = "bedtools getfasta -fi " + in_fa +  " -bed " + in_bed + " " + add_param + " -s  -fo " +  out_fa
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

def plot_seq_len_distr(seq_len_list1, seq_len_list2, plot_out,
                       label1='Set1',
                       label2='Set2',
                       density=True):
    """
    Plot histograms of sequence length distributions.

    """

    plt.figure(figsize=(10, 6))
    plt.hist(seq_len_list1, bins=50, alpha=0.5, label=label1, density=density)
    plt.hist(seq_len_list2, bins=50, alpha=0.5, label=label2, density=density)

    plt.xlabel('Sequence Length')
    plt.ylabel('Density')
    plt.title('Length Distribution of Sequence sets')
    plt.legend(loc='upper right')

    plt.savefig(plot_out)


################################################################################

def read_fasta_into_dic(fasta_file,
                        seqs_dic=False,
                        ids_dic=False,
                        dna=False,
                        full_header=False,
                        report=1,
                        all_uc=False,
                        name_bed=False,
                        empty_check=True,
                        id_check=True,
                        skip_data_id="set",
                        new_header_id="site",
                        remove_regex=False,  # e.g. r"[ :\(\)]"
                        make_uniq_headers=False,
                        skip_n_seqs=True):
    """
    Read in FASTA sequences, store in dictionary and return dictionary.
    FASTA file can be plain text or gzipped (watch out for .gz ending).

    remove_regex:
        If regex given, use this to remove special characters from the header ID.
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
    header_idx = 0

    # Compile regex patterns.
    header_pattern = re.compile(r"^>(.+)" if full_header else r"^>(\S+)")  # \S any non-whitespaces match.
    bed_pattern = re.compile(r"^(.+)::")
    seq_pattern = re.compile(r"[ACGTUN]+", re.I)
    n_pattern = re.compile(r"N", re.I)

    # Open FASTA either as .gz or as text file.
    if re.search(r".+\.gz$", fasta_file):
        f = gzip.open(fasta_file, 'rt')
    else:
        f = open(fasta_file, "r")
    
    for line in f:
        # line = line.strip()
        if line.startswith(">"):

            m = header_pattern.search(line)
            assert m, f'header ID extraction failed for FASTA header line "{line}"'
            seq_id = m.group(1)

            if remove_regex:
                seq_id = re.sub(remove_regex, "", seq_id)
                assert seq_id, "filtering FASTA sequence header ID \"%s\" by regex \"%s\" resulted in string" %(seq_id, remove_regex)

            # If name_bed, get first part of ID (before "::").
            if name_bed:
                m = bed_pattern.search(seq_id)
                assert m, f'BED column 4 ID extraction failed for FASTA header "{seq_id}"'
                seq_id = m.group(1)

            if id_check:
                assert seq_id not in seqs_dic, "non-unique FASTA header \"%s\" in \"%s\"" % (seq_id, fasta_file)

            if make_uniq_headers:
                header_idx += 1
                seq_id = new_header_id + "_" + str(header_idx)

            if ids_dic:
                if seq_id in ids_dic:
                    seqs_dic[seq_id] = ""
            else:
                seqs_dic[seq_id] = ""
            # elif re.search("[ACGTUN]+", line, re.I):
            #     m = re.search("([ACGTUN]+)", line, re.I)
            #     seq = m.group(1)

        elif seq_pattern.search(line):
            m = seq_pattern.search(line)
            seq = m.group(0)
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
            # if re.search("N", seq, re.I):
            if n_pattern.search(seq):
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
                     tr2gid_dic=False,
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
    tr2gid_dic:
        If set add gene ID to header.

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
        if tr2gid_dic:
            if seq_id in tr2gid_dic:
                out_id = out_id + "," + str(tr2gid_dic[seq_id])
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

    motif_blocks_dic entry:
    ['letter-probability matrix: alength= 4 w= 10 nsites= 20 E= 0', 
    ' 0.050000  0.400000  0.050000  0.500000 ', 
    ' 0.285714  0.000000  0.000000  0.714286 ', 
    ' 0.252525  0.141414  0.606061  0.000000 ', 
    ' 0.505051  0.090909  0.404040  0.000000 ', 
    ' 0.714286  0.000000  0.000000  0.285714 ', 
    ' 0.081633  0.000000  0.918367  0.000000 ', 
    ' 0.505051  0.090909  0.404040  0.000000 ', 
    ' 0.505051  0.090909  0.404040  0.000000 ', 
    ' 0.353535  0.191919  0.454545  0.000000 ', 
    ' 0.404040  0.242424  0.353535  0.000000 ']

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

def get_consensus_motif_from_seq_block(seq_block):
    """
    seq_block is list with following format:
    ['letter-probability matrix: alength= 4 w= 10 nsites= 20 E= 0', 
    ' 0.050000  0.400000  0.050000  0.500000 ', 
    ' 0.285714  0.000000  0.000000  0.714286 ', 
    ' 0.252525  0.141414  0.606061  0.000000 ', 
    ' 0.505051  0.090909  0.404040  0.000000 ', 
    ' 0.714286  0.000000  0.000000  0.285714 ', 
    ' 0.081633  0.000000  0.918367  0.000000 ', 
    ' 0.505051  0.090909  0.404040  0.000000 ', 
    ' 0.505051  0.090909  0.404040  0.000000 ', 
    ' 0.353535  0.191919  0.454545  0.000000 ', 
    ' 0.404040  0.242424  0.353535  0.000000 ']
    
    Return consensus motif sequence.

    >>> seq_block = ['letter-probability matrix: alength= 4 w= 10 nsites= 20 E= 0', ' 0.250000  0.250000  0.250000  0.250000 ', ' 0.285714  0.000000  0.000000  0.714286 ', ' 0.081633  0.000000  0.918367  0.000000 ']
    >>> get_consensus_motif_from_seq_block(seq_block)
    'ATG'
    
    """

    idx2nt_dic = {0: "A", 1: "C", 2: "G", 3: "T"}

    assert seq_block, "seq_block empty"
    assert len(seq_block) > 1, "seq_block empty (only header present?)"

    # Get consensus motif sequence.
    consensus_seq = ""
    for i in range(1, len(seq_block)):
        line = seq_block[i]
        line = line.strip()
        if not line:
            continue
        # Split line.
        ll = line.split()
        assert len(ll) == 4, "invalid line format in seq_block (expected 4 columns)"
        # Get max value index.
        max_val = 0
        max_idx = 0
        for j in range(4):
            val = float(ll[j])
            if val > max_val:
                max_val = val
                max_idx = j
        # Get consensus base.
        max_nt = idx2nt_dic[max_idx]
        consensus_seq += max_nt

    return consensus_seq


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
        if re.search(r"^MOTIF\s\w+", l):
            m = re.search(r"MOTIF (\w+)", l)
            motif_id = m.group(1)
            new_motif_id = remove_special_chars_from_str(motif_id)
            assert new_motif_id, "no characters left after removal of special characters from motif ID \"%s\". Please use valid MEME XML motif IDs (i.e., modify MOTIF column strings in motifs xml file)" %(motif_id)
            motif_id = new_motif_id
        else:
            if motif_id and l:
                if re.search(r"^URL", l):  # Skip URL rows e.g. from Ray2013 meme file. format: URL http:// ...
                    continue
                # Also remove <tab> characters from lines.
                l = l.replace("\t", "")
                if motif_id in motif_blocks_dic:
                    motif_blocks_dic[motif_id].append(l)
                else:
                    motif_blocks_dic[motif_id] = [l]

    return motif_blocks_dic


################################################################################

def blocks_to_xml_string(motif_blocks_dic, motif_ids_dic,
                         mid2rid_dic=None,
                         out_xml=False):
    """
    Return MEME XML string based on given motif IDs dictionary and available
    motif blocks dictionary.

    mid2rid_dic:
        motif ID -> RBP ID mapping dictionary.
        Adds RBP ID to motif ID in MEME XML output.

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
            motif_str = motif_id
            if mid2rid_dic is not None:  # if motif ID to RBP ID mapping given.
                if motif_id in mid2rid_dic:
                    motif_str = motif_id + " " + mid2rid_dic[motif_id]
            block_str_final = "MOTIF " + motif_str + " \n" + block_str + "\n\n"
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
        if re.search(r"^chr", chr_id):
            if chr_id in add_chr_names_dic or re.search(r"^chr[\dMXY]+$", chr_id):
                return chr_id
            else:
                return False
        else:
            # Convert to "chr" IDs.
            if chr_id == "MT": # special case MT -> chrM.
                return "chrM"
            if chr_id in add_chr_names_dic or re.search(r"^[\dXY]+$", chr_id):
                return "chr" + chr_id
            else:
                return False

    elif id_style == 2:

        if re.search(r"^chr", chr_id):
            if chr_id == "chrM": # special case chrM -> MT.
                return "MT"
            if chr_id in add_chr_names_dic or re.search(r"^chr[\dXY]+$", chr_id):
                # Cut out chr suffix.
                m = re.search(r"chr(.+)", chr_id)
                assert m, "no match for regex search"
                chr_suffix = m.group(1)
                return chr_suffix
            else:
                return False

        else:
            if chr_id == "MT": # special case MT.
                return chr_id
            if chr_id in add_chr_names_dic or re.search(r"^[\dXY]+$", chr_id):
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

def gtf_tr_feat_to_bed(in_gtf, out_bed, tr_feat,
                       chr_style=0,
                       uniq_reg=False,
                       codon_len_check=True,
                       tr_ids_dic=False):
    """
    Extract transcript feature regions to BED file.

    tr_ids_dic:
        If transcript IDs dictionary given, only extract regions for these.
    chr_style:
        0: do not change
        1: change to chr1, chr2 ...
        2: change to 1, 2, 3, ...
    uniq_reg:
        If True, do not output same genomic region twice.
    codon_len_check:
        If True, length check for stop_codon and start_codon features is
        performed and filtered.

    """

    OUTBED = open(out_bed, "w")

    seen_reg_dic = {}

    if re.search(r".+\.gz$", in_gtf):
        f = gzip.open(in_gtf, 'rt')
    else: 
        f = open(in_gtf, "r")
    
    for line in f:

        if line.startswith("#"):
            continue

        cols = line.strip().split("\t")
        feature = cols[2]
        if feature != tr_feat:
            continue

        chr_id = cols[0]
        feat_s = int(cols[3]) - 1
        feat_e = int(cols[4])
        feat_pol = cols[6]
        infos = cols[8]

        reg_id = chr_id + ":" + str(feat_s) + "-" + str(feat_e) + "(" + feat_pol + ")"
        reg_len = feat_e - feat_s

        if codon_len_check:
            if feature == "start_codon":
                if reg_len != 3:
                    print("WARNING: start_codon region length != 3 in GTF file \"%s\", line \"%s\". Skipping region ..." %(in_gtf, line))
                    continue
            elif feature == "stop_codon":
                if reg_len != 3:
                    print("WARNING: stop_codon region length != 3 in GTF file \"%s\", line \"%s\". Skipping region ..." %(in_gtf, line))
                    continue

        if uniq_reg:
            if reg_id in seen_reg_dic:
                continue
            else:
                seen_reg_dic[reg_id] = 1
    
        chr_id = check_convert_chr_id(chr_id, id_style=chr_style)
        # If not one of standard chromosomes, continue.
        if not chr_id:
            continue

        assert feat_e > feat_s, "feature end <= feature start in GTF file \"%s\", line \"%s\". This should not happen" %(in_gtf, line)

        m = re.search('transcript_id "(.+?)"', infos)
        assert m, "transcript_id entry missing in GTF file \"%s\", line \"%s\"" %(in_gtf, line)
        tr_id = m.group(1)
        if tr_ids_dic:
            if tr_id not in tr_ids_dic:
                continue

        out_id = tr_id + ";" + tr_feat

        OUTBED.write("%s\t%i\t%i\t%s\t0\t%s\n" %(chr_id, feat_s, feat_e, out_id, feat_pol))

    f.close()
    OUTBED.close()


################################################################################

def gtf_read_in_gene_infos(in_gtf,
                           tr2gid_dic=None,
                           tr_types_dic=None,
                           check_chr_ids_dic=None,
                           chr_style=0,
                           skip_gene_biotype_dic=None,
                           gene_ids_dic=False,
                           empty_check=False):
    """
    Read in gene infos into GeneInfo objects, including information on 
    transcript isoforms for the gene. Note that only features on standard 
    chromosomes (1,2,...,X Y MT) are currently used.

    Assuming gtf file with order: gene,transcript(s),exon(s) ...

    tr_types_dic:
        Store transcript biotype IDs and number of appearances.
        transcript biotype ID -> # appearances
    
    gene_ids_dic:
        If gene IDs dictionary given, only extract these genes.
        
    chr_style:
        0: do not change
        1: change to chr1, chr2 ...
        2: change to 1, 2, 3, ...

    """

    # Transcript ID to exonic length dictionary.
    tr2len_dic = {}
    # Gene info objects dictionary (gene_id -> gene info object).
    gid2gio_dic = {}

    if skip_gene_biotype_dic is None:
        skip_gene_biotype_dic = {"TEC" : 1}

    if re.search(r".+\.gz$", in_gtf):
        f = gzip.open(in_gtf, 'rt')
    else: 
        f = open(in_gtf, "r")
    for line in f:
        # Skip header.
        # if re.search("^#", line):
        if line.startswith("#"):
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
            assert chr_id in check_chr_ids_dic, "chromosome ID \"%s\" from GTF file \"%s\" not found in --genome FASTA file. Make sure input GTF + FASTA + BED files use compatible chromosome IDs" %(chr_id, in_gtf)

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

            if gene_ids_dic:
                if gene_id not in gene_ids_dic:
                    continue

            if gene_biotype in skip_gene_biotype_dic:
                continue

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

            # assert gene_id in gid2gio_dic, "gene_id %s belonging to transcript ID %s not (yet) encountered. Gene feature expected to come before transcript and exon features in GTF file \"%s\"" %(gene_id, tr_id, in_gtf)
            if gene_id not in gid2gio_dic:
                continue

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
            # MANE select.
            mane_select = 0
            m = re.search('tag "MANE_Select"', infos)
            if m:
                mane_select = 1
            # Transcript support level (TSL).
            # transcript_support_level "NA (assigned to previous version 1)"
            m = re.search('transcript_support_level "(.+?)"', infos)
            tsl_id = "NA"
            if m:
                tsl_id = m.group(1)
                if re.search("assigned to previous", tsl_id):
                    m = re.search(r"(.+?) \(", tsl_id)
                    tsl_id = m.group(1)
            # Dummy length for now.
            tr_length = 0

            # Update gene infos.
            gid2gio_dic[gene_id].tr_ids.append(tr_id)
            gid2gio_dic[gene_id].tr_biotypes.append(tr_biotype)
            gid2gio_dic[gene_id].tr_basic_tags.append(basic_tag)
            gid2gio_dic[gene_id].tr_ensembl_canonical_tags.append(ensembl_canonical)
            gid2gio_dic[gene_id].tr_mane_select_tags.append(mane_select)
            gid2gio_dic[gene_id].tr_tsls.append(tsl_id)
            gid2gio_dic[gene_id].tr_lengths.append(tr_length)

        elif feature == "exon":
            m = re.search('gene_id "(.+?)"', infos)
            assert m, "gene_id entry missing in GTF file \"%s\", line \"%s\"" %(in_gtf, line)
            gene_id = m.group(1)
            if gene_id not in gid2gio_dic:
                continue
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

    if re.search(r".+\.gz$", in_gtf):
        f = gzip.open(in_gtf, 'rt')
    else: 
        f = open(in_gtf, "r")
    for line in f:
        # Skip header.
        # if re.search("^#", line):
        if line.startswith("#"):
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
            # MANE select.
            mane_select = 0
            m = re.search('tag "MANE_Select"', infos)
            if m:
                mane_select = 1
            # Transcript support level (TSL).
            # transcript_support_level "NA (assigned to previous version 1)"
            m = re.search('transcript_support_level "(.+?)"', infos)
            tsl_id = "NA"
            if m:
                tsl_id = m.group(1)
                if re.search("assigned to previous", tsl_id):
                    m = re.search(r"(.+?) \(", tsl_id)
                    tsl_id = m.group(1)

            if tr_types_dic is not None:
                if tr_biotype not in tr_types_dic:
                    tr_types_dic[tr_biotype] = 1
                else:
                    tr_types_dic[tr_biotype] += 1

            tr_infos = TranscriptInfo(tr_id, tr_biotype, chr_id, feat_s, feat_e, feat_pol, gene_id,
                                      tr_length=0,
                                      basic_tag=basic_tag,  # int
                                      ensembl_canonical=ensembl_canonical,  # int
                                      mane_select=mane_select,  # int
                                      tsl_id=tsl_id,  # int
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

            m = re.search(r'exon_number "(\d+?)"', infos)
            if not m:
                m = re.search(r'exon_number (\d+?);', infos)  # GENCODE encoding.
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

    """
    Add intron coordinates and total intron length to TranscriptInfo objects.

    """

    if tid2tio_dic:
        for tid in tid2tio_dic:        
            tr_pol = tid2tio_dic[tid].tr_pol
            exon_coords = tid2tio_dic[tid].exon_coords
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

            total_intron_len = 0
            for intron in intron_coords:
                intron_len = intron[1] - intron[0] + 1
                total_intron_len += intron_len
                tid2tio_dic[tid].intron_coords.append(intron)

            tid2tio_dic[tid].total_intron_len = total_intron_len

    return tid2tio_dic


################################################################################

def get_exon_pos_count_list_dic(tid2tio_dic,
                                tr_ids_dic=None):
    """
    Get exon ID -> exon position count list dictionary.

    Expects exon_coords in correct order (i.e. exon 1 on - most downstream, 
    exon 1 on + most upstream).
    
    tr_ids_dic:
        If set, only consider transcript IDs in tr_ids_dic.
    
    """

    exon2pcl_dic = {}

    for tr_id in tid2tio_dic:
        
        if tr_ids_dic is not None:
            if tr_id not in tr_ids_dic:
                continue

        exon_coords = tid2tio_dic[tr_id].exon_coords
        assert exon_coords is not None, "exon coordinates list not set for transcript ID %s" %(tr_id)
    
        for idx, exon in enumerate(exon_coords):

            exon_len = exon[1] - exon[0] + 1
            exon_id = "exon;%s;%i" %(tr_id, idx+1)
            exon2pcl_dic[exon_id] = [0] * exon_len

    return exon2pcl_dic


################################################################################

class ExonIntronOverlap:
    """
    Exon intron overlap stats class.

    """
    def __init__(self,
                 dataset_id: str,
                 c_input_sites: int,
                 c_exon_sites = 0,
                 c_intron_sites = 0,
                 c_intergenic_sites = 0,
                 c_us_ib_sites = 0,
                 c_ds_ib_sites = 0,
                 c_us_ib_dist_sites = 0,  # upstream intron border distant (>intron_border_len) sites.
                 c_ds_ib_dist_sites = 0,  # downstream intron border distant (>intron_border_len) sites.
                 c_eib_sites = 0,
                 c_first_exon_sites = 0,
                 c_last_exon_sites = 0,
                 c_single_exon_sites = 0,
                 min_overlap = 0.9,
                 intron_border_len = 250,
                 ei_border_len = 50,
                 c_tr_ids = 0,
                 c_tr_ids_with_sites = 0) -> None:

        self.dataset_id = dataset_id
        self.c_input_sites = c_input_sites
        self.c_exon_sites = c_exon_sites
        self.c_intron_sites = c_intron_sites
        self.c_intergenic_sites = c_intergenic_sites
        self.c_us_ib_sites = c_us_ib_sites
        self.c_ds_ib_sites = c_ds_ib_sites
        self.c_us_ib_dist_sites = c_us_ib_dist_sites
        self.c_ds_ib_dist_sites = c_ds_ib_dist_sites
        self.c_eib_sites = c_eib_sites
        self.c_first_exon_sites = c_first_exon_sites
        self.c_last_exon_sites = c_last_exon_sites
        self.c_single_exon_sites = c_single_exon_sites
        self.min_overlap = min_overlap
        self.intron_border_len = intron_border_len
        self.ei_border_len = ei_border_len
        self.c_tr_ids = c_tr_ids
        self.c_tr_ids_with_sites = c_tr_ids_with_sites


################################################################################

def get_intron_exon_ol_counts(overlap_ei_regions_bed):
    """
    Get exon + intron overlap counts

    overlap_ei_regions_bed format (old):
    chr1	3385286	3396490	intron;ENST00000270722	0	+
    chr1	3396593	3402790	intron;ENST00000270722	0	+

    overlap_ei_regions_bed format:
    $ intersectBed -a intron_exon_regions.tmp.bed -b in_sites.filtered.bed -s -f 0.5 -F 0.5 -e  -wb
    chr1	21727871	21728037	intron;ENST00000308271	0	-	chr1	21727871	21728037	chr1:21727871-21728037(-)	6.18043979583829	-
    chr5	138561019	138561076	intron;ENST00000297185	0	-	chr5	138561019	138561076	chr5:138561019-138561076(-)	3.36528644625518	-
    chr6	16265910	16265991	intron;ENST00000259727	0	+	chr6	16265910	16265991	chr6:16265910-16265991(+)	3.80557689731445	+
    chr6	16265981	16266114	intron;ENST00000259727	0	+	chr6	16265981	16266114	chr6:16265981-16266114(+)	3.39282576293052	+

    """

    assert os.path.exists(overlap_ei_regions_bed), "file \"%s\" does not exist" %(overlap_ei_regions_bed)

    c_exon_ol = 0
    c_intron_ol = 0
    reg2exc_dic = {}  # only count a region once as exon region.
    reg2inc_dic = {}  # only count a region once as intron region.

    with open(overlap_ei_regions_bed, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            cols = line.split("\t")
            reg_tr_id = cols[3]
            reg_id = cols[9]
            reg_annot = reg_tr_id.split(";")[0]
            if reg_annot == "intron":
                if reg_id not in reg2inc_dic:
                    c_intron_ol += 1
                    reg2inc_dic[reg_id] = 1
            else:
                if reg_id not in reg2exc_dic:
                    c_exon_ol += 1
                    reg2exc_dic[reg_id] = 1

    return c_exon_ol, c_intron_ol


################################################################################

def get_eib_ol_counts(overlap_eib_regions_bed):
    """
    Get exon-intron border overlap counts.

    overlap_eib_regions_bed old format:
    chr1	231373978	231374002	eib	0	-
    chr1	244408941	244408985	ds_ib	0	-
    chr1	151266003	151266053	eib	0	+
    chr1	161157671	161157694	eib	0	+
    chr1	110952350	110952391	us_ib	0	-
    chr1	66958933	66958983	eib	0	+

    New format:
    chr1	231373978	231374002	eib	0	-	chr1	231373978	231374002	chr1:231373978-231374002(-)	6.6772932452923	-
    chr1	244408941	244408985	ds_ib	0	-	chr1	244408941	244408985	chr1:244408941-244408985(-)	3.84108057430655	-
    chr1	151266003	151266053	eib	0	+	chr1	151266003	151266058	chr1:151266003-151266058(+)	3.62307325778397	+
    chr1	161157671	161157694	eib	0	+	chr1	161157671	161157694	chr1:161157671-161157694(+)	5.94342578410804	+
    chr1	27907232	27907277	eib	0	-	chr1	27907225	27907277	chr1:27907225-27907277(-)	4.17126552599614	-
    chr1	110407601	110407635	eib	0	-	chr1	110407601	110407660	chr1:110407601-110407660(-)	3.33819938009558	-
    chr1	110952350	110952391	us_ib	0	-	chr1	110952350	110952391	chr1:110952350-110952391(-)	3.3312393913757	-

    """

    assert os.path.exists(overlap_eib_regions_bed), "file \"%s\" does not exist" %(overlap_eib_regions_bed)

    c_eib_ol = 0
    c_us_ib_ol = 0
    c_ds_ib_ol = 0
    reg2eib_dic = {}  # only count a region once as eib region.
    reg2usib_dic = {}  # only count a region once as us_ib region.
    reg2dsib_dic = {}  # only count a region once as ds_ib region.

    with open(overlap_eib_regions_bed, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            cols = line.split("\t")
            reg_annot = cols[3]
            reg_id = cols[9]
            if reg_annot == "eib":
                if reg_id not in reg2eib_dic:
                    c_eib_ol += 1
                    reg2eib_dic[reg_id] = 1
            elif reg_annot == "us_ib":
                if reg_id not in reg2usib_dic:
                    c_us_ib_ol += 1
                    reg2usib_dic[reg_id] = 1
            elif reg_annot == "ds_ib":
                if reg_id not in reg2dsib_dic:
                    c_ds_ib_ol += 1
                    reg2dsib_dic[reg_id] = 1
            else:
                assert False, "invalid region ID \"%s\" found in line \"%s\"" %(reg_annot, line)

    return c_eib_ol, c_us_ib_ol, c_ds_ib_ol


################################################################################

def exon_intron_border_regions_to_bed(tid2tio_dic, out_bed,
                                      tr_ids_dic=None,
                                      intron_border_len=250,
                                      ei_border_len=50):
    """
    Output border regions (intron border regions, and exon-intron border regions)
    to out_bed. 
    and exon border regions to BED.

    These can then be used to determin what regions are covered by
    input sites.

    Correct exon_coords order expected (i.e. exon 1 on - most downstream, and most upstream on +).

    Based on:
    chr1	3000	3100	exon	0	-
    chr1	2000	2100	exon	0	-
    chr1	1000	1100	exon	0	-
    chr1	2100	3000	intron	0	-
    chr1	1100	2000	intron	0	-
    chr1	1000	1100	exon	0	+
    chr1	2000	2100	exon	0	+
    chr1	3000	3100	exon	0	+
    chr1	1100	2000	intron	0	+
    chr1	2100	3000	intron	0	+

    Generate this:
    chr1	1100	1350	us_ib	0	+
    chr1	1750	2000	ds_ib	0	+
    chr1	2100	2350	us_ib	0	+
    chr1	2750	3000	ds_ib	0	+
    chr1	1050	1150	eib	0	+
    chr1	1950	2050	eib	0	+
    chr1	2050	2150	eib	0	+
    chr1	2950	3050	eib	0	+
    chr1	2750	3000	us_ib	0	-
    chr1	2100	2350	ds_ib	0	-
    chr1	1750	2000	us_ib	0	-
    chr1	1100	1350	ds_ib	0	-
    chr1	2050	2150	W	0	-
    chr1	2950	3050	eib	0	-
    chr1	1050	1150	eib	0	-
    chr1	1950	2050	eib	0	-
    
    >>> out_bed = "test_data/test.exon_intron_border_reg.tmp.bed"
    >>> exp_bed = "test_data/test.exon_intron_border_reg.exp.bed"
    >>> tid2tio_dic = {}
    >>> tid2tio_dic["ENST1"] = TranscriptInfoExonTest("ENST1", "chr1", "+", [[1001, 1100], [2001, 2100], [3001, 3100]])
    >>> tid2tio_dic["ENST2"] = TranscriptInfoExonTest("ENST2", "chr1", "-", [[3001, 3100], [2001, 2100], [1001, 1100]])
    >>> exon_intron_border_regions_to_bed(tid2tio_dic, out_bed)
    >>> diff_two_files_identical(out_bed, exp_bed)
    True


    """

    assert tid2tio_dic, "given tid2tio_dic empty"

    OUTBED = open(out_bed, "w")

    min_intron_len = intron_border_len * 2

    for tr_id in tid2tio_dic:

        if tr_ids_dic is not None:
            if tr_id not in tr_ids_dic:
                continue

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
                intron_coords.append([exon_coords[i+1][1]+1, exon_coords[i][0]-1])
        else:
            assert False, "invalid strand given (%s) for transcript ID %s" %(tr_pol, tr_id)

        # # Output exon regions.
        # for idx, exon in enumerate(exon_coords):
        #     c_exon_out += 1
        #     exon_id = "exon"
        #     OUTBED.write("%s\t%i\t%i\t%s\t0\t%s\n" % (chr_id, exon[0]-1, exon[1], exon_id, tr_pol))

        for idx, intron in enumerate(intron_coords):
            # intron_id = "intron"
            intron_s = intron[0]-1
            intron_e = intron[1]
            # # Output intron regions.
            # OUTBED.write("%s\t%i\t%i\t%s\t0\t%s\n" % (chr_id, intron_s, intron_e, intron_id, tr_pol))
            # Intron border regions.
            intron_len = intron_e - intron_s
            if intron_len >= min_intron_len:
                # + case.
                us_intron_s = intron_s
                us_intron_e = intron_s + intron_border_len
                ds_intron_s = intron_e - intron_border_len
                ds_intron_e = intron_e
                us_intron_id = "us_ib"  # upstream intron border region.
                ds_intron_id = "ds_ib"  # downstream intron border region.
                if tr_pol == "-":
                    us_intron_s = intron_e - intron_border_len
                    us_intron_e = intron_e
                    ds_intron_s = intron_s
                    ds_intron_e = intron_s + intron_border_len
                OUTBED.write("%s\t%i\t%i\t%s\t0\t%s\n" % (chr_id, us_intron_s, us_intron_e, us_intron_id, tr_pol))
                OUTBED.write("%s\t%i\t%i\t%s\t0\t%s\n" % (chr_id, ds_intron_s, ds_intron_e, ds_intron_id, tr_pol))

            # Exon-intron-border regions.
            us_exon_s = exon_coords[idx][0] - 1
            if tr_pol == "-":
                us_exon_s = exon_coords[idx+1][0] - 1            
            us_exon_e = intron_s
            us_exon_len = us_exon_e - us_exon_s
            ds_exon_s = intron_e
            ds_exon_e = exon_coords[idx+1][1]
            if tr_pol == "-":
                ds_exon_e = exon_coords[idx][1]
            ds_exon_len = ds_exon_e - ds_exon_s
            us_eib_s = intron_s - ei_border_len
            if us_exon_len < ei_border_len:
                us_eib_s = us_exon_s
            us_eib_e = intron_s + ei_border_len
            if intron_len < ei_border_len:
                us_eib_e = intron_e
            ds_eib_s = intron_e - ei_border_len
            if intron_len < ei_border_len:
                ds_eib_s = intron_s
            ds_eib_e = intron_e + ei_border_len
            if ds_exon_len < ei_border_len:
                ds_eib_e = ds_exon_e
            eib_id = "eib"
            OUTBED.write("%s\t%i\t%i\t%s\t0\t%s\n" % (chr_id, us_eib_s, us_eib_e, eib_id, tr_pol))
            OUTBED.write("%s\t%i\t%i\t%s\t0\t%s\n" % (chr_id, ds_eib_s, ds_eib_e, eib_id, tr_pol))

    OUTBED.close()

"""
	-f	Minimum overlap required as a fraction of A.
		- Default is 1E-9 (i.e., 1bp).
		- FLOAT (e.g. 0.50)

	-F	Minimum overlap required as a fraction of B.
		- Default is 1E-9 (i.e., 1bp).
		- FLOAT (e.g. 0.50)

	-r	Require that the fraction overlap be reciprocal for A AND B.
		- In other words, if -f is 0.90 and -r is used, this requires
		  that B overlap 90% of A and A _also_ overlaps 90% of B.

	-e	Require that the minimum fraction be satisfied for A OR B.
		- In other words, if -e is used with -f 0.90 and -F 0.10 this requires
		  that either 90% of A is covered OR 10% of  B is covered.
		  Without -e, both fractions would have to be satisfied.

"""

################################################################################

def output_eib_pos_to_bed(tid2tio_dic, out_bed,
                          min_intron_len=False,
                          min_exon_len=False,
                          border_mode=1):

    """
    Output exon-intron border positions (first or last intron positions) 
    to out_bed.

    border_mode:
        1: upstream borders
        2: downstream borders

    Correct exon_coords order expected (i.e. exon 1 on - most downstream, and most upstream on +).

    >>> out_bed = "test_data/test.exon_intron_border_pos.tmp.bed"
    >>> exp_bed = "test_data/test.exon_intron_border_pos.exp.bed"
    >>> tid2tio_dic = {}
    >>> tid2tio_dic["ENST1"] = TranscriptInfoExonTest("ENST1", "chr1", "+", [[1001, 1100], [2001, 2100], [3001, 3100]])
    >>> tid2tio_dic["ENST2"] = TranscriptInfoExonTest("ENST2", "chr1", "-", [[3001, 3100], [2001, 2100], [1001, 1100]])
    >>> output_eib_pos_to_bed(tid2tio_dic, out_bed)
    >>> diff_two_files_identical(out_bed, exp_bed)
    True
    >>> exp_bed = "test_data/test.exon_intron_border_pos2.exp.bed"
    >>> output_eib_pos_to_bed(tid2tio_dic, out_bed, border_mode=2)
    >>> diff_two_files_identical(out_bed, exp_bed)
    True

    
    """

    assert tid2tio_dic, "given tid2tio_dic empty"

    OUTBED = open(out_bed, "w")

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
                intron_coords.append([exon_coords[i+1][1]+1, exon_coords[i][0]-1])
        else:
            assert False, "invalid strand given (%s) for transcript ID %s" %(tr_pol, tr_id)

        for idx, intron in enumerate(intron_coords):

            intron_s = intron[0]-1
            intron_e = intron[1]
            eib_id = "eib"

            intron_len = intron_e - intron_s
            if min_intron_len:
                if intron_len < min_intron_len:
                    continue

            # Exon lengths.
            us_exon_s = exon_coords[idx][0] - 1
            if tr_pol == "-":
                us_exon_s = exon_coords[idx+1][0] - 1            
            us_exon_e = intron_s
            us_exon_len = us_exon_e - us_exon_s
            ds_exon_s = intron_e
            ds_exon_e = exon_coords[idx+1][1]
            if tr_pol == "-":
                ds_exon_e = exon_coords[idx][1]
            ds_exon_len = ds_exon_e - ds_exon_s

            us_eib_s = intron_s
            us_eib_e = intron_s + 1
            ds_eib_s = intron_e - 1
            ds_eib_e = intron_e

            if tr_pol == "-":
                us_eib_s = intron_e - 1
                us_eib_e = intron_e
                ds_eib_s = intron_s
                ds_eib_e = intron_s + 1

            if border_mode == 1:
                if min_exon_len:
                    if us_exon_len < min_exon_len:
                        continue
                eib_id += ";us;%i;%s" %(idx+1, tr_id)
                OUTBED.write("%s\t%i\t%i\t%s\t0\t%s\n" % (chr_id, us_eib_s, us_eib_e, eib_id, tr_pol))
            elif border_mode == 2:
                if min_exon_len:
                    if ds_exon_len < min_exon_len:
                        continue
                eib_id += ";ds;%i;%s" %(idx+1, tr_id)
                OUTBED.write("%s\t%i\t%i\t%s\t0\t%s\n" % (chr_id, ds_eib_s, ds_eib_e, eib_id, tr_pol))

    OUTBED.close()


################################################################################

def output_transcript_info_intron_exon_to_bed(tid2tio_dic, out_bed,
                                              output_mode=1,
                                              tr_ids_dic=None,
                                              report_counts=True,
                                              add_tr_id=True,
                                              add_numbers=False,
                                              number_format=1,
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
    add_numbers:
        If True, add exon + intron numbers, with format: current_number-total_number
    number_format:
        1: exon;ENST1;1-3
        2: exon;ENST1;1
    
    """
    assert tid2tio_dic, "given tid2tio_dic empty"
    assert number_format in [1,2], "invalid number_format given"

    OUTBED = open(out_bed, "w")

    c_exon_out = 0
    c_intron_out = 0

    for tr_id in tid2tio_dic:

        if tr_ids_dic is not None:
            if tr_id not in tr_ids_dic:
                continue

        chr_id = tid2tio_dic[tr_id].chr_id
        tr_pol = tid2tio_dic[tr_id].tr_pol
        exon_c = tid2tio_dic[tr_id].exon_c
        exon_coords = tid2tio_dic[tr_id].exon_coords
        assert exon_coords is not None, "exon coordinates list not set for transcript ID %s" %(tr_id)
    
        # Get intron coordinates.
        intron_coords = []
        intron_c = 0
        if tr_pol == "+":
            for i in range(len(exon_coords) - 1):
                intron_coords.append([exon_coords[i][1]+1, exon_coords[i+1][0]-1])
        elif tr_pol == "-":
            for i in range(len(exon_coords) - 1):
                intron_coords.append([exon_coords[i+1][1]+1,exon_coords[i][0]-1])
        else:
            assert False, "invalid strand given (%s) for transcript ID %s" %(tr_pol, tr_id)

        if intron_coords:
            intron_c = len(intron_coords)
        assert exon_c == intron_c + 1, "exon count (%i) does not match intron count (%i) + 1 for transcript ID %s" %(exon_c, intron_c, tr_id)

        if output_mode == 1:
            for idx, exon in enumerate(exon_coords):
                c_exon_out += 1
                exon_id = "exon"
                if add_tr_id:
                    exon_id += ";" + tr_id
                if add_numbers:
                    if number_format == 1:
                        exon_id += ";%i-%i" %(idx+1, exon_c)
                    elif number_format == 2:
                        exon_id += ";%i" %(idx+1)
                OUTBED.write("%s\t%i\t%i\t%s\t0\t%s\n" % (chr_id, exon[0]-1, exon[1], exon_id, tr_pol))

            for idx, intron in enumerate(intron_coords):
                c_intron_out += 1
                intron_id = "intron"
                if add_tr_id:
                    intron_id += ";" + tr_id
                if add_numbers:
                    if number_format == 1:
                        intron_id += ";%i-%i" %(idx+1, intron_c)
                    elif number_format == 2:
                        intron_id += ";%i" %(idx+1)
                OUTBED.write("%s\t%i\t%i\t%s\t0\t%s\n" % (chr_id, intron[0]-1, intron[1], intron_id, tr_pol))

        elif output_mode == 2:
            for idx, exon in enumerate(exon_coords):
                c_exon_out += 1
                exon_id = "exon"
                if add_tr_id:
                    exon_id += ";" + tr_id
                if add_numbers:
                    if number_format == 1:
                        exon_id += ";%i-%i" %(idx+1, exon_c)
                    elif number_format == 2:
                        exon_id += ";%i" %(idx+1)
                OUTBED.write("%s\t%i\t%i\t%s\t0\t%s\n" % (chr_id, exon[0]-1, exon[1], exon_id, tr_pol))

        elif output_mode == 3:
            for idx, intron in enumerate(intron_coords):
                c_intron_out += 1
                intron_id = "intron"
                if add_tr_id:
                    intron_id += ";" + tr_id
                if add_numbers:
                    if number_format == 1:
                        intron_id += ";%i-%i" %(idx+1, intron_c)
                    elif number_format == 2:
                        intron_id += ";%i" %(idx+1)
                OUTBED.write("%s\t%i\t%i\t%s\t0\t%s\n" % (chr_id, intron[0]-1, intron[1], intron_id, tr_pol))

    OUTBED.close()

    if report_counts:
        print("# of exon features output to BED:   ", c_exon_out)
        print("# of intron features output to BED: ", c_intron_out)
    if empty_check:
        assert (c_exon_out+c_intron_out) > 0, "no exon/intron features output to BED"


################################################################################

def get_mrna_tids_and_sites(in_bed):
    """
    Given BED file with format:
    chr1	1080	1085	exon;ENST01;1	0	-	chr1	1070	1085	s1	0	-
    chr2	1010	1020	exon;ENST02;1	0	+	chr2	1010	1030	s2	0	+
    chr2	1015	1020	exon;ENST02;2	0	+	chr2	1015	1030	s3	0	+

    Extract transcript IDs (ENST01,..) and count site IDs (s1,...).
    
    """

    c_ol_mrna_sites = 0
    ol_mrna_tids_dic = {}
    seen_site_ids_dic = {}

    with open(in_bed, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            cols = line.split("\t")
            exon_id = cols[3]
            tr_id = exon_id.split(";")[1]
            site_id = cols[9]
            ol_mrna_tids_dic[tr_id] = 1
            seen_site_ids_dic[site_id] = 1

    c_ol_mrna_sites = len(seen_site_ids_dic)

    return c_ol_mrna_sites, ol_mrna_tids_dic


################################################################################

def fill_exon_pos_count_lists(exon_cov_bed, tid2tio_dic, exon2pcl_dic):
    """
    Fill position count lists for exons in exon2pcl_dic.

    exon_cov_bed has following format:
    chr1	1080	1085	exon;ENST01;1	0	-
    chr2	1010	1020	exon;ENST02;1	0	+
    chr2	1010	1020	exon;ENST02;1	0	+

    tid2tio_dic["ENST01"].exon_coords format (e1, e2, e3 order):
    ENST01 exon_coords: [[1081, 1100], [1041, 1050], [1001, 1010]]
    exon2pcl_dic format:
    {..., 'exon;2;ENST01': [0, 0, 0, 0, 0, 0, 0, 0, 0, 0], ...}

    test.exon_cov.bed:
    chr1	1050	1055	exon;2;ENST01	0	-
    chr1	1050	1054	exon;2;ENST01	0	-
    chr1	1050	1053	exon;2;ENST01	0	-
    chr1	1050	1052	exon;2;ENST01	0	-
    chr1	1050	1051	exon;2;ENST01	0	-
    chr2	1005	1010	exon;1;ENST02	0	+
    chr2	1007	1010	exon;1;ENST02	0	+
    chr2	1050	1053	exon;2;ENST02	0	+

    >>> exon_cov_bed = "test_data/test.exon_cov.bed"
    >>> tid2tio_dic = {}
    >>> tid2tio_dic["ENST02"] = TranscriptInfoExonTest("ENST02", "chr2", "+", [[1001, 1010], [1051, 1055]])
    >>> exon2pcl_dic = {'exon;ENST02;1': [0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 'exon;ENST02;2': [0, 0, 0, 0, 0], 'exon;666;1': [0]}
    >>> fill_exon_pos_count_lists(exon_cov_bed, tid2tio_dic, exon2pcl_dic)
    >>> exon2pcl_dic
    {'exon;ENST02;1': [0, 0, 0, 0, 0, 1, 1, 2, 2, 2], 'exon;ENST02;2': [1, 1, 1, 0, 0], 'exon;666;1': [0]}
    >>> tid2tio_dic = {}
    >>> tid2tio_dic["ENST01"] = TranscriptInfoExonTest("ENST01", "chr1", "-", [[1096, 1100], [1051, 1060]])
    >>> exon2pcl_dic = {'exon;ENST01;1': [0, 0, 0, 0, 0], 'exon;ENST01;2': [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]}
    >>> fill_exon_pos_count_lists(exon_cov_bed, tid2tio_dic, exon2pcl_dic)
    >>> exon2pcl_dic
    {'exon;ENST01;1': [0, 0, 0, 0, 0], 'exon;ENST01;2': [0, 0, 0, 0, 0, 1, 2, 3, 4, 5]}

    """

    assert exon2pcl_dic, "exon2pcl_dic is empty"
    assert tid2tio_dic, "tid2tio_dic is empty"
    assert os.path.exists(exon_cov_bed), "exon_cov_bed does not exist"

    with open(exon_cov_bed, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            cols = line.split("\t")
            chr_id = cols[0]
            hit_s = int(cols[1])  # 0-based.
            hit_e = int(cols[2])  # 1-based.
            ex_id = cols[3]
            ex_strand = cols[5]

            tr_id = ex_id.split(";")[1]
            
            if tr_id not in tid2tio_dic:
                continue

            ex_nr = int(ex_id.split(";")[2])
            ex_s = tid2tio_dic[tr_id].exon_coords[ex_nr-1][0]
            ex_e = tid2tio_dic[tr_id].exon_coords[ex_nr-1][1]

            if ex_id not in exon2pcl_dic:
                continue

            for i in range(hit_s, hit_e):
                if ex_strand == "+":
                    exon2pcl_dic[ex_id][i - ex_s + 1] += 1
                else:
                    exon2pcl_dic[ex_id][ex_e - i - 1] += 1


################################################################################

def count_lines_in_file(file_path):
    """
    Count number of lines in file.

    >>> file_path = "test_data/file1"
    >>> count_lines_in_file(file_path)
    9
    >>> file_path = "test_data/empty_file"
    >>> count_lines_in_file(file_path)
    0

    """
    with open(file_path, 'r') as file:
        line_count = sum(1 for line in file)
    return line_count


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

def bed_intersect_files_count_lines(a_bed, b_bed,
                                    params="-s -u"):
    """
    Intersect two files using bedtools intersect.
    Count number of lines in output.

    """
    assert os.path.exists(a_bed), "a_bed does not exist"
    assert os.path.exists(b_bed), "b_bed does not exist"

    check_cmd = "intersectBed -a " + a_bed + " -b " + b_bed + " " + params + " | wc -l"
    output = subprocess.getoutput(check_cmd)

    count = 0
    for line in output.split('\n'):
        count = int(line.strip())
        break

    assert count >= 0, "intersectBed + line count has problems with your input:\n%s\n%s" %(check_cmd, output)

    return count


################################################################################

class TranscriptInfoExonTest:
    """
    Transcript infos exon coordinates test class.

    """
    def __init__(self,
                 tr_id: str,
                 chr_id: str,
                 tr_pol: str,
                 exon_coords=None) -> None:

        self.tr_id = tr_id
        self.chr_id = chr_id
        self.tr_pol = tr_pol
        if exon_coords is None:
            self.exon_coords = []
        else:
            self.exon_coords = exon_coords



################################################################################

class MrnaRegionProfile:
    """
    mRNA input region coverage profile class. 

    """
    def __init__(self,
                 dataset_id: str,
                 utr5_len_norm: float,
                 cds_len_norm: float,
                 utr3_len_norm: float,
                 norm_mode: str,
                 c_ol_sites: Optional[int] = None,
                 c_all_sites: Optional[int] = None,
                 c_ol_mrnas: Optional[int] = None,
                 utr5_pc_list=None,
                 cds_pc_list=None,
                 utr3_pc_list=None) -> None:

        self.dataset_id = dataset_id
        self.utr5_len_norm = utr5_len_norm
        self.cds_len_norm = cds_len_norm
        self.utr3_len_norm = utr3_len_norm
        self.norm_mode = norm_mode
        self.c_ol_sites = c_ol_sites
        self.c_all_sites = c_all_sites
        self.c_ol_mrnas = c_ol_mrnas
        if utr5_pc_list is None:
            self.utr5_pc_list = []
        else:
            self.utr5_pc_list = utr5_pc_list
        if cds_pc_list is None:
            self.cds_pc_list = []
        else:
            self.cds_pc_list = cds_pc_list
        if utr3_pc_list is None:
            self.utr3_pc_list = []
        else:
            self.utr3_pc_list = utr3_pc_list


################################################################################

def get_mrna_reg_norm_len(tid2regl_dic,
                          mrna_norm_mode=1):
    """
    Get normalized mRNA region lengths.

    mrna_norm_mode:
        1: median
        2: mean
    
    """

    assert tid2regl_dic, "given tid2regl_dic empty"

    # Get mean / median UTR CDS lengths.
    utr5_len_list = []
    cds_len_list = []
    utr3_len_list = []
    for tr_id in tid2regl_dic:
        utr5_len_list.append(tid2regl_dic[tr_id][0])
        cds_len_list.append(tid2regl_dic[tr_id][1])
        utr3_len_list.append(tid2regl_dic[tr_id][2])
    
    utr5_len_norm = 100
    cds_len_norm = 100
    utr3_len_norm = 100
    norm_mode = "uniform"

    if mrna_norm_mode == 1:
        # Median.
        utr5_len_norm = statistics.median(utr5_len_list)
        cds_len_norm = statistics.median(cds_len_list)
        utr3_len_norm = statistics.median(utr3_len_list)
        norm_mode = "median"
        # print("Median lengths of mRNA regions:")

    elif mrna_norm_mode == 2:
        # Mean.
        utr5_len_norm = statistics.mean(utr5_len_list)
        cds_len_norm = statistics.mean(cds_len_list)
        utr3_len_norm = statistics.mean(utr3_len_list)
        norm_mode = "mean"
        # print("Mean lengths of mRNA regions:")

    else:
        assert False, "invalid --mrna-norm-mode %i set" %(mrna_norm_mode)

    # print("5'UTR = ", utr5_len_norm)
    # print("CDS   = ", cds_len_norm)
    # print("3'UTR = ", utr3_len_norm)

    return utr5_len_norm, cds_len_norm, utr3_len_norm, norm_mode


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
                 total_intron_len = 0,
                 tr_length: Optional[int] = None,  # spliced transcript length.
                 exon_c: Optional[int] = None,
                 basic_tag: Optional[int] = None,
                 ensembl_canonical: Optional[int] = None,
                 mane_select: Optional[int] = None,
                 tsl_id: Optional[str] = None,
                 cds_s: Optional[int] = None,
                 cds_e: Optional[int] = None,
                 intron_coords=None,  # intron_coords + exon_coords with 1-based starts and ends.
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
        self.total_intron_len = total_intron_len
        self.basic_tag = basic_tag
        self.ensembl_canonical = ensembl_canonical
        self.mane_select = mane_select
        self.tsl_id = tsl_id

        if intron_coords is None:
            self.intron_coords = []
        else:
            self.intron_coords = intron_coords
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
                 tr_mane_select_tags=None,
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
        if tr_mane_select_tags is None:
            self.tr_mane_select_tags = []
        else:
            self.tr_mane_select_tags = tr_mane_select_tags
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

def get_transcript_sequences_from_gtf(tid2tio_dic, in_genome_fasta,
                                      tr_ids_dic=False,
                                      dna=True,
                                      all_uc=True,
                                      tmp_out_folder=False):
    """
    Given tid2tio_dic, extract transcript sequences from genome FASTA file.
    
    Return dictionary with transcript ID -> sequence mapping.

    tr_ids_dic:
        Defines transcript IDs for which sequences should be extracted.
    
    """

    assert tid2tio_dic, "given tid2tio_dic empty"

    tr_seqs_dic = {}

    tids_to_extract = []
    if tr_ids_dic:
        tids_to_extract = list(tr_ids_dic.keys())
    else:
        for tid in tid2tio_dic:
            tids_to_extract.append(tr_id)
    

    # Generate .tmp files.
    tmp_bed = "exon_regions.tmp.bed"
    tmp_fa = "exon_regions.tmp.fa"
    if tmp_out_folder:
        tmp_bed = tmp_out_folder + "/" + tmp_bed
        tmp_fa = tmp_out_folder + "/" + tmp_fa


    """
    Output exon regions to BED.

    """

    OUTEXONBED = open(tmp_bed, "w")

    extracted_exon_ids_dic = {}

    for tr_id in tids_to_extract:

        assert tr_id in tid2tio_dic, "transcript ID \"%s\" not in tid2tio_dic" %(tr_id)
        tio = tid2tio_dic[tr_id]

        chr_id = tio.chr_id
        tr_pol = tio.tr_pol

        for idx, exon in enumerate(tio.exon_coords):
            exon_s = exon[0] - 1
            exon_e = exon[1]
            idx += 1
            exon_id = tr_id + "_e" + str(idx)
            extracted_exon_id = exon_id + "::" + chr_id + ":" + str(exon_s) + "-" + str(exon_e) + "(" + tr_pol + ")"
            extracted_exon_ids_dic[exon_id] = extracted_exon_id

            OUTEXONBED.write("%s\t%i\t%i\t%s\t0\t%s\n" % (chr_id, exon_s, exon_e, exon_id, tr_pol))

    OUTEXONBED.close()

    """
    Get exon region sequences from genome FASTA.

    chr6	113868012	113873351	mrocki_201	0	-
    -name option allows to keep the BED name field in the FASTA header.
    Unfortunately additional crap is added as well to the header (no matter what option,
    -name, nameOnly ...)
    >mrocki_201::chr6:113868012-113873351(-)
    
    """

    bed_extract_sequences_from_fasta(tmp_bed, in_genome_fasta, tmp_fa,
                                     add_param="-name",
                                     print_warnings=False,
                                     ignore_errors=False)

    exon_seqs_dic = read_fasta_into_dic(tmp_fa,
                                        dna=dna,
                                        all_uc=all_uc,
                                        skip_n_seqs=False)

    """
    Concatenate exon region sequences to transcript sequences.

    """

    for tr_id in tids_to_extract:
        ex_c = tid2tio_dic[tr_id].exon_c

        for i in range(ex_c):
            i += 1
            ex_id = tr_id + "_e" + str(i)
            extr_ex_id = extracted_exon_ids_dic[ex_id]

            if extr_ex_id in exon_seqs_dic:
                ex_seq = exon_seqs_dic[extr_ex_id]
                if tr_id not in tr_seqs_dic:
                    tr_seqs_dic[tr_id] = ex_seq
                else:
                    tr_seqs_dic[tr_id] += ex_seq
            else:
                print("WARNING: no sequence extracted for exon ID \"%s\". Skipping \"%s\" .. " %(ex_id, tr_id))
                if tr_id in tr_seqs_dic:
                    del tr_seqs_dic[tr_id]
                break

    assert tr_seqs_dic, "tr_seqs_dic empty (no FASTA sequences extracted?)"
    for tr_id in tr_seqs_dic:
        tr_len = len(tr_seqs_dic[tr_id])
        exp_len = tid2tio_dic[tr_id].tr_length
        assert tr_len == exp_len, "BED transcript length != FASTA transcript length for \"%s\"" %(tr_id)

    # Take out the trash.
    if os.path.exists(tmp_bed):
        os.remove(tmp_bed)
    if os.path.exists(tmp_fa):
        os.remove(tmp_fa)

    return tr_seqs_dic


################################################################################

def output_mrna_regions_to_bed(tid2regl_dic, mrna_regions_bed,
                               empty_check=True):
    """
    Given dictionary with transcript ID -> list of 5'UTR, CDS, 3'UTR lengths,
    output the UTR and CDS regions to BED file.
    
    """

    if empty_check:
        assert tid2regl_dic, "given tid2regl_dic empty"

    OUTREGBED = open(mrna_regions_bed, "w")

    for tid in tid2regl_dic:
        trl = tid2regl_dic[tid]
        utr5_s = 0
        utr5_e = trl[0]
        cds_s = utr5_e
        cds_e = cds_s + trl[1]
        utr3_s = cds_e
        utr3_e = utr3_s + trl[2]

        utr5l = utr5_e - utr5_s
        cds_l = cds_e - cds_s
        utr3_l = utr3_e - utr3_s

        if utr5l > 0:
            OUTREGBED.write("%s\t%i\t%i\t%s;5'UTR\t0\t+\n" % (tid, utr5_s, utr5_e, tid))
        if cds_l > 0:
            OUTREGBED.write("%s\t%i\t%i\t%s;CDS\t0\t+\n" % (tid, cds_s, cds_e, tid))
        if utr3_l > 0:
            OUTREGBED.write("%s\t%i\t%i\t%s;3'UTR\t0\t+\n" % (tid, utr3_s, utr3_e, tid))

    OUTREGBED.close()


################################################################################

def get_mrna_region_lengths(tid2tio_dic):
    """
    Given a dictionary of TranscriptInfo objects (tid2tio_dic), calculate
    5'UTR, CDS, and 3'UTR lengths for each transcript.
    Return dictionary with transcript ID as key and list of 5'UTR, CDS, 3'UTR
    lengths as value.

    """
    tid2regl_dic = {}

    for tid in tid2tio_dic:
        tio = tid2tio_dic[tid]

        if tio.cds_s is not None:
            
            tid2regl_dic[tid] = [0, 0, 0]

            for exon in tio.exon_coords:
                cds_se, utr5_se, utr3_se, cds, utr5, utr3 = get_cds_exon_overlap(tio.cds_s, tio.cds_e, exon[0], exon[1], strand=tio.tr_pol)
                if cds:
                    cds_len = cds_se[1] - cds_se[0] + 1
                    tid2regl_dic[tid][1] += cds_len
                if utr5:
                    utr5_len = utr5_se[1] - utr5_se[0] + 1
                    tid2regl_dic[tid][0] += utr5_len
                if utr3:
                    utr3_len = utr3_se[1] - utr3_se[0] + 1
                    tid2regl_dic[tid][2] += utr3_len

    return tid2regl_dic


################################################################################

def bed_filter_by_seqs_dic(seqs_dic, in_bed, out_bed,
                           use_col4_id=False):
    """
    Keep only in_bed entries with column 4 ID present in seqs_dic.

    """

    assert os.path.isfile(in_bed), "in_bed file does not exist"
    assert seqs_dic, "seqs_dic empty"
    c_out = 0

    OUTBED = open(out_bed, "w")

    with open(in_bed) as f:
        for line in f:
            cols = line.strip().split("\t")
            chr_id = cols[0]
            reg_s = cols[1]
            reg_e = cols[2]
            reg_id = cols[3]
            sc = cols[4]
            strand = cols[5]

            exp_reg_id = chr_id + ":" + reg_s + "-" + reg_e + "(" + strand + ")"
            if use_col4_id:
                exp_reg_id = reg_id

            if exp_reg_id in seqs_dic:
                OUTBED.write("%s\t%s\t%s\t%s\t0\t%s\n" % (chr_id, reg_s, reg_e, exp_reg_id, strand))
                c_out += 1

    f.closed
    OUTBED.close()
    assert c_out, "no entries written to out_bed"


################################################################################

def output_exon_annotations(tid2tio_dic, out_bed,
                            custom_annot_dic=None,
                            add_numbers=False,
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

    for tr_id in tid2tio_dic:
        tio = tid2tio_dic[tr_id]
        exon_c = tio.exon_c
        cds_label = "CDS;%s" %(tr_id)
        utr5_label = "5'UTR;%s" %(tr_id)
        utr3_label = "3'UTR;%s" %(tr_id)

        if tio.cds_s is not None:
            for idx, exon in enumerate(tio.exon_coords):

                exon_nr_label = ""
                if add_numbers:
                    exon_nr = idx + 1
                    exon_nr_label =  ";%i-%i" %(exon_nr, exon_c)

                cds_se, utr5_se, utr3_se, cds, utr5, utr3 = get_cds_exon_overlap(tio.cds_s, tio.cds_e, exon[0], exon[1], strand=tio.tr_pol)
                if cds:
                    OUTEXAN.write("%s\t%i\t%i\t%s%s\t0\t%s\n" %(tio.chr_id, cds_se[0]-1, cds_se[1], cds_label, exon_nr_label, tio.tr_pol))
                if utr5:
                    OUTEXAN.write("%s\t%i\t%i\t%s%s\t0\t%s\n" %(tio.chr_id, utr5_se[0]-1, utr5_se[1], utr5_label, exon_nr_label, tio.tr_pol))
                if utr3:
                    OUTEXAN.write("%s\t%i\t%i\t%s%s\t0\t%s\n" %(tio.chr_id, utr3_se[0]-1, utr3_se[1], utr3_label, exon_nr_label, tio.tr_pol))
        else:

            label = other_annot
            if tio.tr_biotype in valid_annot_dic:
                label = valid_annot_dic[tio.tr_biotype]
            label += ";%s" %(tr_id)

            for idx, exon in enumerate(tio.exon_coords):

                exon_nr_label = ""
                if add_numbers:
                    exon_nr = idx + 1
                    exon_nr_label =  ";%i-%i" %(exon_nr, exon_c)

                OUTEXAN.write("%s\t%i\t%i\t%s%s\t0\t%s\n" %(tio.chr_id, exon[0]-1, exon[1], label, exon_nr_label, tio.tr_pol))

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

def get_regex_hit_region_annotations(overlap_annotations_bed,
                                     tid2tio_dic=None):
    """
    Given a BED file (overlap_annotations_bed) with format:
    chr12	57586611	57586617	AATAAA	0	+	chr12	57584217	57586633	3'UTR;ENST00000455537	0	+	6
    chr4	73075045	73075051	AATAAA	0	-	.	-1	-1	.	-1	.	0
    ...
    extract the region annotations and return a dictionary with region ID as key.
    Note that regex regions without overlapping annotation are expected to have
    "." in the annotation field, and will get the annotation "intergenic".
    
    Format reg2annot_dic:
    chr1:101-120(+) -> ["exon", "ENST00000455537"]

    tid2tio_dic:
        If provided, if two features have same overlap amount with region,
        choose the more prominent transcript ID (i.e., better annotation).

    """

    assert os.path.exists(overlap_annotations_bed), "file %s does not exist" %(overlap_annotations_bed)

    reg2maxol_dic = {}
    reg2annot_dic = {}
    annot_col = 9
    c_ol_nt_col = 12

    with open(overlap_annotations_bed) as f:
        for line in f:
            cols = line.strip().split("\t")

            chr_id = cols[0]
            reg_s = str(int(cols[1]) + 1)
            reg_e = cols[2]
            regex = cols[3]
            reg_strand = cols[5]

            reg_id = chr_id + ":" + reg_s + "-" + reg_e + "(" + reg_strand + ")"
            annot_id = "intergenic"
            tr_id = False

            if cols[annot_col] != ".":
                annot_ids = cols[annot_col].split(";")
                assert len(annot_ids) == 3, "len(annot_ids) != 3 (expected ; separated string, but got: \"%s\")" %(cols[9])
                annot_id = annot_ids[0]
                tr_id = annot_ids[1]

            c_overlap_nts = int(cols[c_ol_nt_col])

            if reg_id not in reg2maxol_dic:
                reg2maxol_dic[reg_id] = c_overlap_nts
                reg2annot_dic[reg_id] = [annot_id, tr_id]
            else:
                stored_c_overlap_nts = reg2maxol_dic[reg_id]
                if c_overlap_nts > stored_c_overlap_nts:
                    reg2maxol_dic[reg_id] = c_overlap_nts
                    reg2annot_dic[reg_id][0] = annot_id
                    reg2annot_dic[reg_id][1] = tr_id
                elif c_overlap_nts == stored_c_overlap_nts:  # if region overlaps with > 1 feature by same amount.
                    if tid2tio_dic is not None:
                        best_tid = reg2annot_dic[reg_id][1]
                        new_best_tid = select_more_prominent_tid(tr_id, best_tid, tid2tio_dic)
                        # If current tr_id has better annotation, update region annotation.
                        if new_best_tid == tr_id:
                            reg2maxol_dic[reg_id] = c_overlap_nts
                            reg2annot_dic[reg_id][0] = annot_id
                            reg2annot_dic[reg_id][1] = new_best_tid

    f.closed

    return reg2annot_dic


################################################################################

def get_normnalized_annot_counts(filtered_sites_bed, intron_exon_out_bed,
                                 rbp2motif2annot2c_dic, reg2pol_dic, out_folder):

    """
    Get normalized annotation counts, i.e., counts normalized by the 
    total length of input regions with the given annotations.

    filtered_sites_bed:
    chr8	126557892	127292862	chr8:126557892-127292862(+)	0.0	+
    chr2	20248736	20350844	chr2:20248736-20350844(-)	0.0	-

    intron_exon_out_bed:
    chr8	126713674	126877671	intron;ENST00000645463	0	+
    chr8	126877790	127006554	intron;ENST00000645463	0	+
    chr8	127006618	127292811	intron;ENST00000645463	0	+
    chr2	20350596	20350844	5'UTR;ENST00000361078	0	-
    chr2	20327309	20327360	CDS;ENST00000361078	0	-
    chr2	20327360	20327378	5'UTR;ENST00000361078	0	-
    chr2	20318536	20318645	CDS;ENST00000361078	0	-
    ...
    chr2	20256032	20256170	CDS;ENST00000361078	0	-
    chr2	20255215	20255341	CDS;ENST00000361078	0	-
    chr2	20254862	20254984	CDS;ENST00000361078	0	-
    chr2	20253821	20254014	CDS;ENST00000361078	0	-
    chr2	20251587	20251716	CDS;ENST00000361078	0	-
    chr2	20248736	20251587	3'UTR;ENST00000361078	0	-
    chr8	126557892	126557931	lncRNA;ENST00000645463	0	+
    chr8	126589191	126589290	lncRNA;ENST00000645463	0	+

    rbp2motif2annot2c_dic:
        rbp_id -> motif_id -> annot -> annot_c.
    
    $ intersectBed -a test_reg.bed -b test_annot.bed -s -wb
    chr1	1000	1100	s1	0	+	chr1	900	1100	5'UTR	0	+
    chr1	1100	1300	s1	0	+	chr1	1100	1300	intron	0	+
    chr1	1300	1500	s1	0	+	chr1	1300	1500	CDS	0	+
    chr1	1500	1700	s1	0	+	chr1	1500	1700	intron	0	+
    chr1	1700	1900	s1	0	+	chr1	1700	1900	3'UTR	0	+
    chr1	1600	1650	s1	0	+	chr1	1600	1650	lncRNA	0	+

    """

    assert os.path.exists(filtered_sites_bed), "filtered_sites_bed does not exist"
    assert os.path.exists(intron_exon_out_bed), "intron_exon_out_bed does not exist"
    assert rbp2motif2annot2c_dic, "rbp2motif2annot2c_dic empty"
    
    rbp2motif2annot2normc_dic = {}

    in_bed = out_folder + "/in_sites.filtered.eff_regs.tmp.bed"
    reg_len_dic = {}
    bed_get_effective_reg_bed(filtered_sites_bed, in_bed, reg2pol_dic,
                              reg_len_dic=reg_len_dic)
                              
    eff_reg_size = 0
    for reg_id in reg_len_dic:
        reg_len = reg_len_dic[reg_id]
        eff_reg_size += reg_len

    # # Get effective total input regions length.
    # eff_reg_size = get_uniq_gen_size(filtered_sites_bed)

    out_bed = out_folder + "/in_sites.filtered.annot_overlap.tmp.bed"

    bed_intersect_files(in_bed, intron_exon_out_bed, out_bed, params="-s -wb")

    annot_len_dic = {}

    with open(out_bed) as f:
        for line in f:
            cols = line.strip().split("\t")
            reg_s = int(cols[1])
            reg_e = int(cols[2])
            reg_id = cols[3]
            reg_len = reg_e - reg_s
            annot_info = cols[9]
            annot_id = annot_info.split(";")[0]
            tr_id = annot_info.split(";")[1]
            if annot_id not in annot_len_dic:
                annot_len_dic[annot_id] = reg_len
            else:
                annot_len_dic[annot_id] += reg_len

    f.closed

    # Get effective region length of of regions overlapping with annotations.
    eff_ol_reg_size = get_uniq_gen_size(out_bed)
    assert eff_ol_reg_size <= eff_reg_size, "eff_ol_reg_size > eff_reg_size"
    len_intergenic = eff_reg_size - eff_ol_reg_size
    annot_len_dic["intergenic"] = len_intergenic

    for rbp_id in rbp2motif2annot2c_dic:
        rbp2motif2annot2normc_dic[rbp_id] = {}
        for motif_id in rbp2motif2annot2c_dic[rbp_id]:
            rbp2motif2annot2normc_dic[rbp_id][motif_id] = {}
            for annot_id in rbp2motif2annot2c_dic[rbp_id][motif_id]:
                assert annot_id in annot_len_dic, "annot_id \"%s\" not in annot_len_dic" %(annot_id)
                annot_c = rbp2motif2annot2c_dic[rbp_id][motif_id][annot_id]
                annot_len_1000 = annot_len_dic[annot_id] / 1000
                norm_c = annot_c / annot_len_1000
                rbp2motif2annot2normc_dic[rbp_id][motif_id][annot_id] = norm_c
    
    return rbp2motif2annot2normc_dic, eff_reg_size


################################################################################

def bed_get_effective_reg_bed(in_bed, out_bed, reg2pol_dic,
                              reg2sc_dic=False,
                              reg_len_dic=None):
    """
    Convert regions BED file into effective regions BED file.
    Effective regions are unique regions, i.e., overlapping regions are merged.

    reg2pol_dic:
        region ID -> region polarity / strand. Needed since merge operation does not 
        store region polarity/ strand.

    >>> in_bed = "test_data/test.eff_reg.in.bed"
    >>> out_bed = "test_data/test.eff_reg.tmp.bed"
    >>> exp_bed = "test_data/test.eff_reg.exp.bed"
    >>> reg2pol_dic = {"s1":"+", "s2":"+", "s3":"+", "s4":"-", "s5":"-"}
    >>> bed_get_effective_reg_bed(in_bed, out_bed, reg2pol_dic)
    >>> diff_two_files_identical(out_bed, exp_bed)
    True
         
    """

    params_str = '-s -c 4 -o distinct -delim ";"'
    check_cmd = "sort -k1,1 -k2,2n " + in_bed + " | mergeBed -i stdin " + params_str
    output = subprocess.getoutput(check_cmd)

    BEDOUT = open(out_bed, "w")
    c_out = 0

    for line in output.split('\n'):
        cols = line.strip().split("\t")
        chr_id = cols[0]
        reg_s = int(cols[1])
        reg_e = int(cols[2])
        reg_ids_str = cols[3]
        reg_ids_list = reg_ids_str.split(";")
        reg_id_0 = reg_ids_list[0]
        assert reg_id_0 in reg2pol_dic, "region ID \"%s\" not in reg2pol_dic" %(reg_id_0)
        reg_strand = reg2pol_dic[reg_id_0]
        reg_len = reg_e - reg_s
        if reg_len_dic is not None:
            reg_len_dic[reg_ids_str] = reg_len
        c_out += 1
        comb_reg_sc = 0  # Combined region score.
        reg_sc_list = []
        if reg2sc_dic:
            for reg_id in reg_ids_list:
                reg_sc_list.append(float(reg2sc_dic[reg_id]))
            comb_reg_sc = sum(reg_sc_list) / len(reg_sc_list)

        BEDOUT.write("%s\t%i\t%i\t%s\t%s\t%s\n" % (chr_id, reg_s, reg_e, reg_ids_str, str(comb_reg_sc), reg_strand))

    BEDOUT.close()

    assert c_out, "mergeBed on in_bed %s to produce out_bed %s failed (in_bed empty?)" %(in_bed, out_bed)


################################################################################

def output_promoter_regions_to_bed(tid2tio_dic, out_bed,
                                   prom_min_tr_len=False,
                                   prom_mrna_only=False,
                                   mrna_biotype_label="protein_coding",
                                   prom_ext_up=1000,
                                   prom_ext_down=100,
                                   add_annot_stats_dic=None):
    """
    Output promoter regions to BED file, based on transcript start sites.

    """
    assert tid2tio_dic, "tid2tio_dic empty"

    OUTBED = open(out_bed, "w")

    c_filt_min_tr_len = 0
    c_filt_mrna_only = 0
    c_out = 0

    for tid in tid2tio_dic:
        tio = tid2tio_dic[tid]
        chr_id = tio.chr_id
        tr_pol = tio.tr_pol

        # Optionally filter by minimum transcript length and biotype.
        tr_length = tio.tr_length
        if prom_min_tr_len:
            if tr_length < prom_min_tr_len:
                c_filt_min_tr_len += 1
                continue

        tr_biotype = tio.tr_biotype
        if prom_mrna_only:
            if tr_biotype != mrna_biotype_label:
                c_filt_mrna_only += 1
                continue

        # Define putative promoter region.
        tts_s = tio.tr_s - 1  # make 0-based.
        tts_e = tts_s + 1
        prom_s = tts_s - prom_ext_up
        prom_e = tts_e + prom_ext_down
        if tr_pol == "-":  # minus strand case.
            tts_s = tio.tr_e - 1
            tts_e = tts_s + 1
            prom_s = tts_s - prom_ext_down
            prom_e = tts_e + prom_ext_up

        if prom_s < 0:
            prom_s = 0

        c_out += 1

        OUTBED.write("%s\t%i\t%i\t%s\t0\t%s\n" % (chr_id, prom_s, prom_e, tid, tr_pol))


    OUTBED.close()
    assert c_out, "no promoter regions written to out_bed %s" %(out_bed)

    if add_annot_stats_dic is not None:
        add_annot_stats_dic["c_promoters"] = c_out
        add_annot_stats_dic["c_filt_min_tr_len"] = c_filt_min_tr_len
        add_annot_stats_dic["c_filt_mrna_only"] = c_filt_mrna_only


################################################################################

def output_gene_regions_to_bed(gid2gio_dic, out_bed,
                               add_annot_stats_dic=None):
    """
    Output gene regions stored in gid2gio_dic to BED file.

    """
    assert gid2gio_dic, "gid2gio_dic empty"

    OUTBED = open(out_bed, "w")
    c_out = 0

    for gene_id in gid2gio_dic:
        gene_info = gid2gio_dic[gene_id]
        chr_id = gene_info.chr_id
        gene_s = gene_info.gene_s - 1  # make 0-based.
        gene_e = gene_info.gene_e
        gene_pol = gene_info.gene_pol

        c_out += 1

        OUTBED.write("%s\t%i\t%i\t%s\t0\t%s\n" % (chr_id, gene_s, gene_e, gene_id, gene_pol))

    OUTBED.close()
    assert c_out, "no gene regions written to out_bed %s" %(out_bed)
    if add_annot_stats_dic is not None:
        add_annot_stats_dic["c_genes"] = c_out


################################################################################

def bed_sort_file(in_bed, out_bed, 
                  params="-k1,1 -k2,2n"):
    """
    Sort BED file (necessary for some bedtools operations).

    """

    assert os.path.exists(in_bed), "in_bed does not exist"

    check_cmd = "sort " + params + " " + in_bed + " > " + out_bed
    output = subprocess.getoutput(check_cmd)
    error = False
    if output:
        error = True
    assert error == False, "sort has problems with your input:\n%s\n%s" %(check_cmd, output)


# ################################################################################

# def bed_merge_entries(in_bed, out_bed, 
#                       params="-k1,1 -k2,2n"):

#     """
#     Merge entries in in_bed file. Output is 3-col BED:
#     chr1	1000	1900

#     """

#     assert os.path.exists(in_bed), "in_bed does not exist"

#     check_cmd = "sort " + params + " " + in_bed + " > " + out_bed
#     output = subprocess.getoutput(check_cmd)
#     error = False
#     if output:
#         error = True
#     assert error == False, "sort has problems with your input:\n%s\n%s" %(check_cmd, output)



################################################################################

def get_motif_hit_region_annotations(overlap_annotations_bed,
                                     tid2tio_dic=None):
    """
    Get motif hit region annotations from overlap_annotations_bed. This file 
    was produced using intersectBed -wao, so also motif hit regions that do not
    overlap with genomic regions are present in the file.
    motif hit region id format: HNRNPL:HNRNPL_1;1;method_id:data_id
    genomic annotation region id format: intron;ENST00000434296

    tid2tio_dic:
        If provided, if two features have same overlap amount with region,
        choose the more prominent transcript ID (i.e., better annotation).

    """

    assert os.path.exists(overlap_annotations_bed), "file %s does not exist" %(overlap_annotations_bed)

    reg2maxol_dic = {}
    reg2annot_dic = {}
    # These shift since motif hits BED contains additional (4) p-value and score columns.
    annot_col = 13
    c_ol_nt_col = 16

    motif_id_pattern = re.compile(r"^.+?:(.+?);")

    with open(overlap_annotations_bed) as f:
        for line in f:
            cols = line.strip().split("\t")

            chr_id = cols[0]
            reg_s = str(int(cols[1]) + 1)
            reg_e = cols[2]
            reg_id = cols[3]
            reg_strand = cols[5]

            # col[3] has format: "rbp_id:motif_id;1;method_id:data_id". Extract motif_id from this string.
            # m = re.search("^.+?:(.+?);", reg_id)
            # assert m is not None, "Motif ID extraction failed for region ID \"%s\"" %(reg_id)
            # motif_id = m.group(1)
            m = motif_id_pattern.search(reg_id)
            assert m is not None, "Motif ID extraction failed for region ID \"%s\"" %(reg_id)
            motif_id = m.group(1)

            # motif_id = reg_id.split(":")[1].split(";")[0]
            reg_id = chr_id + ":" + reg_s + "-" + reg_e + "(" + reg_strand + ")" + motif_id

            annot_id = "intergenic"
            tr_id = False

            if cols[annot_col] != ".":
                annot_ids = cols[annot_col].split(";")
                assert len(annot_ids) == 2, "len(annot_ids) != 2 (expected ; separated string, but got: \"%s\")" %(cols[9])
                annot_id = annot_ids[0]
                tr_id = annot_ids[1]

            c_overlap_nts = int(cols[c_ol_nt_col])

            if reg_id not in reg2maxol_dic:
                reg2maxol_dic[reg_id] = c_overlap_nts
                reg2annot_dic[reg_id] = [annot_id, tr_id]
            else:
                stored_c_overlap_nts = reg2maxol_dic[reg_id]
                if c_overlap_nts > stored_c_overlap_nts:
                    reg2maxol_dic[reg_id] = c_overlap_nts
                    reg2annot_dic[reg_id][0] = annot_id
                    reg2annot_dic[reg_id][1] = tr_id

                elif c_overlap_nts == stored_c_overlap_nts:  # if region overlaps with > 1 feature by same amount.
                    if tid2tio_dic is not None:
                        best_tid = reg2annot_dic[reg_id][1]
                        new_best_tid = select_more_prominent_tid(tr_id, best_tid, tid2tio_dic)
                        # If current tr_id has better annotation, update region annotation.
                        if new_best_tid == tr_id:
                            reg2maxol_dic[reg_id] = c_overlap_nts
                            reg2annot_dic[reg_id][0] = annot_id
                            reg2annot_dic[reg_id][1] = new_best_tid

    f.closed

    return reg2annot_dic


################################################################################

def reformat_to_bed10(in_bed, out_bed):
    """
    Reformat to BED10 format.
    
    """

    assert os.path.exists(in_bed), "in_bed file does not exist"
    OUTBED = open(out_bed, "w")

    with open(in_bed) as f:
        for line in f:
            cols = line.strip().split("\t")
            OUTBED.write("%s\t%s\t%s\t%s\t%s\t%s\t-1.0\t-1.0\t-1.0\t-1.0\n" %(cols[0], cols[1], cols[2], cols[3], cols[4], cols[5]))
    f.closed
    OUTBED.close()


################################################################################

def output_bed6plus_to_bed6_format(in_bed, out_bed):
    """
    Output BED file with 6+ columns to BED file with 6 columns.

    Can be necessary for intersectBed to work.

    """

    assert os.path.exists(in_bed), "in_bed file does not exist"
    OUTBED = open(out_bed, "w")

    with open(in_bed) as f:
        for line in f:
            cols = line.strip().split("\t")
            assert len(cols) >= 6, "less than 6 columns in line \"%s\"" %(line)
            OUTBED.write("%s\t%s\t%s\t%s\t%s\t%s\n" %(cols[0], cols[1], cols[2], cols[3], cols[4], cols[5]))
    f.closed
    OUTBED.close()


################################################################################

def get_mrna_region_annotations_v2(overlap_annotations_bed,
                                   reg_ids_dic=None):
    """
    Get mRNA region annotations, i.e. sites are on transcripts and were overlapped
    with 5'UTR CDS 3'UTR regions.

    This function (v2) is designed for rbpbench enmo.

    """
    assert os.path.exists(overlap_annotations_bed), "file %s does not exist" %(overlap_annotations_bed)

    reg2maxol_dic = {}
    reg2annot_dic = {}
    annot_col = 9
    c_ol_nt_col = 12

    with open(overlap_annotations_bed) as f:
        for line in f:
            cols = line.strip().split("\t")
            reg_id = cols[3]  # reg_id format: ENST00000663363:36-136(+)

            annot_ids = cols[annot_col].split(";")  # annot_col format: ENST00000663363;5'UTR
            assert len(annot_ids) == 2, "len(annot_ids) != 2 (expected ; separated string, but got: \"%s\")" %(cols[9])
            tr_id = annot_ids[0]
            annot_id = annot_ids[1]

            c_overlap_nts = int(cols[c_ol_nt_col])

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
                tr_id = reg_id.split(":")[0]
                reg2annot_dic[reg_id] = ["ncRNA", tr_id]

    return reg2annot_dic


################################################################################

def get_dist_to_next_border(site_s, site_e, reg_s, reg_e, reg_strand,
                            skip_sites_not_within=True):
    """
    Get distance of site that is within a region to closest border of the region. 
    Take the center of the site to calculate the distance. Also return whether
    the closest border is upstream (i.e. 5' end) or downstream (i.e. 3' end) of
    the site.
    If site center position is not within region, return -1, "-"

    site_s + reg_s are 0-based coordinates.
    site_e + reg_e are 1-based coordinates.

    skip_sites_not_within:
        If True, return -1, "-" if site center position not within region.

    >>> get_dist_to_next_border(11, 20, 5, 20, "+")
    (3, 'down')
    >>> get_dist_to_next_border(11, 20, 5, 20, "-")
    (3, 'up')
    >>> get_dist_to_next_border(10, 11, 5, 20, "+")
    (5, 'up')
    >>> get_dist_to_next_border(10, 11, 5, 20, "-")
    (5, 'down')
    >>> get_dist_to_next_border(9, 10, 4, 15, "+")
    (5, 'down')
    >>> get_dist_to_next_border(9, 10, 4, 15, "-")
    (5, 'down')
    
    """
    # Get site center position.
    site_cp = get_center_position(site_s, site_e)  # 1-based center.

    reg_s += 1  # make 1-based.

    if site_cp < reg_s or site_cp > reg_e:
        if skip_sites_not_within:
            return -1, "-"
        else:
            assert False, "site center position not within region (site_cp: %i, region: %i-%i,%s)" %(site_cp, reg_s, reg_e, reg_strand)

    # Get distance of site_cp to closest region border.
    dist_up = site_cp - reg_s
    dist_down = reg_e - site_cp
    if reg_strand == "-":
        dist_up = reg_e - site_cp
        dist_down = site_cp - reg_s

    assert dist_up >= 0, "dist_up < 0 (dist_up = %i, site_cp: %i, reg_s: %i, reg_e: %i, reg_strand: %s)" %(dist_up, site_cp, reg_s, reg_e, reg_strand)
    assert dist_down >= 0, "dist_down < 0 (dist_down = %i, site_cp: %i, reg_s: %i, reg_e: %i, reg_strand: %s)" %(dist_down, site_cp, reg_s, reg_e, reg_strand)

    if dist_up < dist_down:
        return dist_up, "up"
    else:
        return dist_down, "down"


################################################################################

def get_eib_annot_c_strict(annot_list, eib_annot_c_dic,
                           ib_len=250,
                           eib_len=50):
    """
    Get exon+intron+eib overlap counts for example annot_list (with ib_strict in dictionary).
    Add counts to eib_annot_c_dic.
    Format:
    eib_annot_c_dic = {
        "exonic" : 0,
        "intronic" : 0,
        "intergenic" : 0,
        "eib" : 0,
        "us_ib_strict" : 0,
        "ds_ib_strict" : 0,
        "us_ib" : 0,
        "ds_ib" : 0,
        "first_exon" : 0,
        "last_exon" : 0,
        "single_exon" : 0
    }

    To add, counts in first + last exon (for > 1 exon transcripts).

    reg2annot_dic format: reg_id -> 
    [annot_id, tr_id, border_dist, us_ds_label, annot_reg_len, exon_intron_nr, gene_id, gene_name, tr_biotype]
    For intergenic regions:
    ["intergenic", False, -1, "-", -1, "-", "-", "-", "-"]
    For intronic/exonic regions with center outside:
    [annot_id, tr_id, -1, "-", annot_reg_len, exon_intron_nr, gene_id, gene_name, tr_biotype]

    >>> eib_annot_c_dic = {"exonic" : 0, "intronic" : 0, "intergenic" : 0, "eib" : 0, "us_ib_strict" : 0, "ds_ib_strict" : 0, "us_ib" : 0, "ds_ib" : 0, "first_exon" : 0, "last_exon" : 0, "single_exon" : 0}
    >>> annot_list = ["intron", "ENST1", 10, "down", 100, "1-1", "ENSG1", "GENE1", "protein_coding"]
    >>> get_eib_annot_c_strict(annot_list, eib_annot_c_dic)
    >>> eib_annot_c_dic
    {'exonic': 0, 'intronic': 1, 'intergenic': 0, 'eib': 1, 'us_ib_strict': 0, 'ds_ib_strict': 0, 'us_ib': 0, 'ds_ib': 1, 'first_exon': 0, 'last_exon': 0, 'single_exon': 0}
    >>> annot_list = ["intron", "ENST1", 10, "down", 1000, "1-1", "ENSG1", "GENE1", "protein_coding"]
    >>> get_eib_annot_c_strict(annot_list, eib_annot_c_dic)
    >>> eib_annot_c_dic
    {'exonic': 0, 'intronic': 2, 'intergenic': 0, 'eib': 2, 'us_ib_strict': 0, 'ds_ib_strict': 1, 'us_ib': 0, 'ds_ib': 2, 'first_exon': 0, 'last_exon': 0, 'single_exon': 0}
    >>> annot_list = ["exon", "ENST1", 10, "down", 1000, "1-1", "ENSG1", "GENE1", "protein_coding"]
    >>> get_eib_annot_c_strict(annot_list, eib_annot_c_dic)
    >>> eib_annot_c_dic
    {'exonic': 1, 'intronic': 2, 'intergenic': 0, 'eib': 2, 'us_ib_strict': 0, 'ds_ib_strict': 1, 'us_ib': 0, 'ds_ib': 2, 'first_exon': 0, 'last_exon': 0, 'single_exon': 1}
    >>> annot_list = ["exon", "ENST1", 10, "down", 1000, "1-2", "ENSG1", "GENE1", "protein_coding"]
    >>> get_eib_annot_c_strict(annot_list, eib_annot_c_dic)
    >>> eib_annot_c_dic
    {'exonic': 2, 'intronic': 2, 'intergenic': 0, 'eib': 3, 'us_ib_strict': 0, 'ds_ib_strict': 1, 'us_ib': 0, 'ds_ib': 2, 'first_exon': 1, 'last_exon': 0, 'single_exon': 1}
    >>> annot_list = ["exon", "ENST1", 10, "down", 1000, "2-2", "ENSG1", "GENE1", "protein_coding"]
    >>> get_eib_annot_c_strict(annot_list, eib_annot_c_dic)
    >>> eib_annot_c_dic
    {'exonic': 3, 'intronic': 2, 'intergenic': 0, 'eib': 3, 'us_ib_strict': 0, 'ds_ib_strict': 1, 'us_ib': 0, 'ds_ib': 2, 'first_exon': 1, 'last_exon': 1, 'single_exon': 1}
    >>> annot_list = ["intergenic", False, -1, "-", -1, "0-0", "-", "-", "-"]
    >>> get_eib_annot_c_strict(annot_list, eib_annot_c_dic)
    >>> eib_annot_c_dic
    {'exonic': 3, 'intronic': 2, 'intergenic': 1, 'eib': 3, 'us_ib_strict': 0, 'ds_ib_strict': 1, 'us_ib': 0, 'ds_ib': 2, 'first_exon': 1, 'last_exon': 1, 'single_exon': 1}
    >>> annot_list = ["exon", "ENST1", 10, "up", 1000, "2-3", "ENSG1", "GENE1", "protein_coding"]
    >>> get_eib_annot_c_strict(annot_list, eib_annot_c_dic)
    >>> eib_annot_c_dic
    {'exonic': 4, 'intronic': 2, 'intergenic': 1, 'eib': 4, 'us_ib_strict': 0, 'ds_ib_strict': 1, 'us_ib': 0, 'ds_ib': 2, 'first_exon': 1, 'last_exon': 1, 'single_exon': 1}
    >>> annot_list = ["exon", "ENST1", 60, "up", 1000, "2-3", "ENSG1", "GENE1", "protein_coding"]
    >>> get_eib_annot_c_strict(annot_list, eib_annot_c_dic)
    >>> eib_annot_c_dic
    {'exonic': 5, 'intronic': 2, 'intergenic': 1, 'eib': 4, 'us_ib_strict': 0, 'ds_ib_strict': 1, 'us_ib': 0, 'ds_ib': 2, 'first_exon': 1, 'last_exon': 1, 'single_exon': 1}
    >>> annot_list = ["intron", "ENST1", 100, "up", 1000, "1-1", "ENSG1", "GENE1", "protein_coding"]
    >>> get_eib_annot_c_strict(annot_list, eib_annot_c_dic)
    >>> eib_annot_c_dic
    {'exonic': 5, 'intronic': 3, 'intergenic': 1, 'eib': 4, 'us_ib_strict': 1, 'ds_ib_strict': 1, 'us_ib': 1, 'ds_ib': 2, 'first_exon': 1, 'last_exon': 1, 'single_exon': 1}
    >>> annot_list = ["intron", "ENST1", 300, "up", 1000, "1-1", "ENSG1", "GENE1", "protein_coding"]
    >>> get_eib_annot_c_strict(annot_list, eib_annot_c_dic)
    >>> eib_annot_c_dic
    {'exonic': 5, 'intronic': 4, 'intergenic': 1, 'eib': 4, 'us_ib_strict': 1, 'ds_ib_strict': 1, 'us_ib': 1, 'ds_ib': 2, 'first_exon': 1, 'last_exon': 1, 'single_exon': 1}

    """

    assert eib_annot_c_dic, "eib_annot_c_dic empty"
    assert len(annot_list) == 9, "len(annot_list) != 9 for annot_list \"%s\"" %(annot_list)
    exp_intron_len = ib_len * 2

    annot_id = annot_list[0]
    border_dist = annot_list[2]
    us_ds_label = annot_list[3]
    annot_reg_len = annot_list[4]
    exon_intron_nr_c = annot_list[5]
    exon_intron_nr = int(exon_intron_nr_c.split("-")[0])
    exon_intron_c = int(exon_intron_nr_c.split("-")[1])

    if annot_id == "intron":
        assert exon_intron_nr > 0, "exon_intron_nr <= 0 for annot_list \"%s\"" %(annot_list)
        eib_annot_c_dic["intronic"] += 1
        if border_dist >= 0:
            if border_dist <= ib_len:
                if us_ds_label == "up":
                    eib_annot_c_dic["us_ib"] += 1
                elif us_ds_label == "down":
                    eib_annot_c_dic["ds_ib"] += 1
                if annot_reg_len >= exp_intron_len:
                    if us_ds_label == "up":
                        eib_annot_c_dic["us_ib_strict"] += 1
                    elif us_ds_label == "down":
                        eib_annot_c_dic["ds_ib_strict"] += 1
            if border_dist <= eib_len:
                eib_annot_c_dic["eib"] += 1

    elif annot_id == "intergenic":
        eib_annot_c_dic["intergenic"] += 1
    else:  # exonic sites.
        assert exon_intron_nr > 0, "exon_intron_nr <= 0 for annot_list \"%s\"" %(annot_list)
        eib_annot_c_dic["exonic"] += 1
        # Only intron-exon borders check ...
        if border_dist >= 0:
            # If in first exon.
            if exon_intron_nr == 1 and exon_intron_c > 1:
                eib_annot_c_dic["first_exon"] += 1
            # If in last exon.
            if exon_intron_nr == exon_intron_c and exon_intron_c > 1:
                eib_annot_c_dic["last_exon"] += 1
            # If on single exon.
            if exon_intron_c == 1:
                eib_annot_c_dic["single_exon"] += 1
            else:
                # If first exon, only look at distance to 3' end.
                if exon_intron_nr == 1:
                    if us_ds_label == "down":
                        if border_dist <= eib_len:
                            eib_annot_c_dic["eib"] += 1
                elif exon_intron_nr == exon_intron_c:  # last exon, only look at distance to 5' end.
                    if us_ds_label == "up":
                        if border_dist <= eib_len:
                            eib_annot_c_dic["eib"] += 1
                else:  # exon in the middle, look at both borders.
                    if border_dist <= eib_len:
                        eib_annot_c_dic["eib"] += 1


################################################################################

def get_eib_annot_c(annot_list, eib_annot_c_dic,
                    ib_len=250,
                    eib_len=50):
    """
    Get exon+intron+eib overlap counts for example annot_list (with ib_dist in dictionary).
    Add counts to eib_annot_c_dic.
    Format:
    eib_annot_c_dic = {
        "exonic" : 0,
        "intronic" : 0,
        "intergenic" : 0,
        "eib" : 0,
        "us_ib_dist" : 0,
        "ds_ib_dist" : 0,
        "us_ib" : 0,
        "ds_ib" : 0,
        "first_exon" : 0,
        "last_exon" : 0,
        "single_exon" : 0
    }

    To add, counts in first + last exon (for > 1 exon transcripts).

    reg2annot_dic format: reg_id -> 
    [annot_id, tr_id, border_dist, us_ds_label, annot_reg_len, exon_intron_nr, gene_id, gene_name, tr_biotype]
    For intergenic regions:
    ["intergenic", False, -1, "-", -1, "-", "-", "-", "-"]
    For intronic/exonic regions with center outside:
    [annot_id, tr_id, -1, "-", annot_reg_len, exon_intron_nr, gene_id, gene_name, tr_biotype]

    >>> eib_annot_c_dic = {"exonic" : 0, "intronic" : 0, "intergenic" : 0, "eib" : 0, "us_ib_dist" : 0, "ds_ib_dist" : 0, "us_ib" : 0, "ds_ib" : 0, "first_exon" : 0, "last_exon" : 0, "single_exon" : 0}
    >>> annot_list = ["intron", "ENST1", 10, "down", 100, "1-1", "ENSG1", "GENE1", "protein_coding"]
    >>> get_eib_annot_c(annot_list, eib_annot_c_dic)
    >>> eib_annot_c_dic
    {'exonic': 0, 'intronic': 1, 'intergenic': 0, 'eib': 1, 'us_ib_dist': 0, 'ds_ib_dist': 0, 'us_ib': 0, 'ds_ib': 1, 'first_exon': 0, 'last_exon': 0, 'single_exon': 0}
    >>> annot_list = ["intron", "ENST1", 10, "down", 1000, "1-1", "ENSG1", "GENE1", "protein_coding"]
    >>> get_eib_annot_c(annot_list, eib_annot_c_dic)
    >>> eib_annot_c_dic
    {'exonic': 0, 'intronic': 2, 'intergenic': 0, 'eib': 2, 'us_ib_dist': 0, 'ds_ib_dist': 0, 'us_ib': 0, 'ds_ib': 2, 'first_exon': 0, 'last_exon': 0, 'single_exon': 0}
    >>> annot_list = ["exon", "ENST1", 10, "down", 1000, "1-1", "ENSG1", "GENE1", "protein_coding"]
    >>> get_eib_annot_c(annot_list, eib_annot_c_dic)
    >>> eib_annot_c_dic
    {'exonic': 1, 'intronic': 2, 'intergenic': 0, 'eib': 2, 'us_ib_dist': 0, 'ds_ib_dist': 0, 'us_ib': 0, 'ds_ib': 2, 'first_exon': 0, 'last_exon': 0, 'single_exon': 1}
    >>> annot_list = ["exon", "ENST1", 10, "down", 1000, "1-2", "ENSG1", "GENE1", "protein_coding"]
    >>> get_eib_annot_c(annot_list, eib_annot_c_dic)
    >>> eib_annot_c_dic
    {'exonic': 2, 'intronic': 2, 'intergenic': 0, 'eib': 3, 'us_ib_dist': 0, 'ds_ib_dist': 0, 'us_ib': 0, 'ds_ib': 2, 'first_exon': 1, 'last_exon': 0, 'single_exon': 1}
    >>> annot_list = ["exon", "ENST1", 10, "down", 1000, "2-2", "ENSG1", "GENE1", "protein_coding"]
    >>> get_eib_annot_c(annot_list, eib_annot_c_dic)
    >>> eib_annot_c_dic
    {'exonic': 3, 'intronic': 2, 'intergenic': 0, 'eib': 3, 'us_ib_dist': 0, 'ds_ib_dist': 0, 'us_ib': 0, 'ds_ib': 2, 'first_exon': 1, 'last_exon': 1, 'single_exon': 1}
    >>> annot_list = ["intergenic", False, -1, "-", -1, "0-0", "-", "-", "-"]
    >>> get_eib_annot_c(annot_list, eib_annot_c_dic)
    >>> eib_annot_c_dic
    {'exonic': 3, 'intronic': 2, 'intergenic': 1, 'eib': 3, 'us_ib_dist': 0, 'ds_ib_dist': 0, 'us_ib': 0, 'ds_ib': 2, 'first_exon': 1, 'last_exon': 1, 'single_exon': 1}
    >>> annot_list = ["exon", "ENST1", 10, "up", 1000, "2-3", "ENSG1", "GENE1", "protein_coding"]
    >>> get_eib_annot_c(annot_list, eib_annot_c_dic)
    >>> eib_annot_c_dic
    {'exonic': 4, 'intronic': 2, 'intergenic': 1, 'eib': 4, 'us_ib_dist': 0, 'ds_ib_dist': 0, 'us_ib': 0, 'ds_ib': 2, 'first_exon': 1, 'last_exon': 1, 'single_exon': 1}
    >>> annot_list = ["exon", "ENST1", 60, "up", 1000, "2-3", "ENSG1", "GENE1", "protein_coding"]
    >>> get_eib_annot_c(annot_list, eib_annot_c_dic)
    >>> eib_annot_c_dic
    {'exonic': 5, 'intronic': 2, 'intergenic': 1, 'eib': 4, 'us_ib_dist': 0, 'ds_ib_dist': 0, 'us_ib': 0, 'ds_ib': 2, 'first_exon': 1, 'last_exon': 1, 'single_exon': 1}
    >>> annot_list = ["intron", "ENST1", 100, "up", 1000, "1-1", "ENSG1", "GENE1", "protein_coding"]
    >>> get_eib_annot_c(annot_list, eib_annot_c_dic)
    >>> eib_annot_c_dic
    {'exonic': 5, 'intronic': 3, 'intergenic': 1, 'eib': 4, 'us_ib_dist': 0, 'ds_ib_dist': 0, 'us_ib': 1, 'ds_ib': 2, 'first_exon': 1, 'last_exon': 1, 'single_exon': 1}
    >>> annot_list = ["intron", "ENST1", 300, "up", 1000, "1-1", "ENSG1", "GENE1", "protein_coding"]
    >>> get_eib_annot_c(annot_list, eib_annot_c_dic)
    >>> eib_annot_c_dic
    {'exonic': 5, 'intronic': 4, 'intergenic': 1, 'eib': 4, 'us_ib_dist': 1, 'ds_ib_dist': 0, 'us_ib': 1, 'ds_ib': 2, 'first_exon': 1, 'last_exon': 1, 'single_exon': 1}
    >>> annot_list = ["intron", "ENST1", 300, "down", 1000, "1-1", "ENSG1", "GENE1", "protein_coding"]
    >>> get_eib_annot_c(annot_list, eib_annot_c_dic)
    >>> eib_annot_c_dic
    {'exonic': 5, 'intronic': 5, 'intergenic': 1, 'eib': 4, 'us_ib_dist': 1, 'ds_ib_dist': 1, 'us_ib': 1, 'ds_ib': 2, 'first_exon': 1, 'last_exon': 1, 'single_exon': 1}
    >>> annot_list = ["exon", "ENST1", 10, "down", 50, "3-3", "ENSG1", "GENE1", "protein_coding"]
    >>> get_eib_annot_c(annot_list, eib_annot_c_dic)
    >>> eib_annot_c_dic
    {'exonic': 6, 'intronic': 5, 'intergenic': 1, 'eib': 5, 'us_ib_dist': 1, 'ds_ib_dist': 1, 'us_ib': 1, 'ds_ib': 2, 'first_exon': 1, 'last_exon': 2, 'single_exon': 1}
    >>> annot_list = ["exon", "ENST1", 10, "up", 50, "1-3", "ENSG1", "GENE1", "protein_coding"]
    >>> get_eib_annot_c(annot_list, eib_annot_c_dic)
    >>> eib_annot_c_dic
    {'exonic': 7, 'intronic': 5, 'intergenic': 1, 'eib': 6, 'us_ib_dist': 1, 'ds_ib_dist': 1, 'us_ib': 1, 'ds_ib': 2, 'first_exon': 2, 'last_exon': 2, 'single_exon': 1}

    """

    assert eib_annot_c_dic, "eib_annot_c_dic empty"
    assert len(annot_list) == 9, "len(annot_list) != 9 for annot_list \"%s\"" %(annot_list)
    # exp_intron_len = ib_len * 2

    annot_id = annot_list[0]
    border_dist = annot_list[2]
    us_ds_label = annot_list[3]
    annot_reg_len = annot_list[4]
    exon_intron_nr_c = annot_list[5]
    exon_intron_nr = int(exon_intron_nr_c.split("-")[0])
    exon_intron_c = int(exon_intron_nr_c.split("-")[1])

    if annot_id == "intron":
        assert exon_intron_nr > 0, "exon_intron_nr <= 0 for annot_list \"%s\"" %(annot_list)
        eib_annot_c_dic["intronic"] += 1
        if border_dist >= 0:
            if border_dist <= ib_len:
                if us_ds_label == "up":
                    eib_annot_c_dic["us_ib"] += 1
                elif us_ds_label == "down":
                    eib_annot_c_dic["ds_ib"] += 1
            else:
                if us_ds_label == "up":
                    eib_annot_c_dic["us_ib_dist"] += 1
                elif us_ds_label == "down":
                    eib_annot_c_dic["ds_ib_dist"] += 1
            if border_dist <= eib_len:
                eib_annot_c_dic["eib"] += 1

    elif annot_id == "intergenic":
        eib_annot_c_dic["intergenic"] += 1
    else:  # exonic sites.
        assert exon_intron_nr > 0, "exon_intron_nr <= 0 for annot_list \"%s\"" %(annot_list)
        eib_annot_c_dic["exonic"] += 1
        if border_dist >= 0:
            # If in first exon.
            if exon_intron_nr == 1 and exon_intron_c > 1:
                eib_annot_c_dic["first_exon"] += 1
            # If in last exon.
            if exon_intron_nr == exon_intron_c and exon_intron_c > 1:
                eib_annot_c_dic["last_exon"] += 1
            # If on single exon.
            if exon_intron_c == 1:
                eib_annot_c_dic["single_exon"] += 1
            else:
                if exon_intron_nr == 1:  # First exon.
                    if us_ds_label == "down":
                        if border_dist <= eib_len:
                            eib_annot_c_dic["eib"] += 1
                    if us_ds_label == "up":
                        if annot_reg_len - border_dist <= eib_len:
                            eib_annot_c_dic["eib"] += 1
                elif exon_intron_nr == exon_intron_c:  # Last exon.
                    if us_ds_label == "up":
                        if border_dist <= eib_len:
                            eib_annot_c_dic["eib"] += 1
                    if us_ds_label == "down":
                        if annot_reg_len - border_dist <= eib_len:
                            eib_annot_c_dic["eib"] += 1
                else:  # Exon in the middle.
                    if border_dist <= eib_len:
                        eib_annot_c_dic["eib"] += 1


################################################################################

def get_region_annotations(overlap_annotations_bed,
                           tid2tio_dic,
                           motif_hits=False,
                           reg_ids_dic=None):
    """
    Get region annotations from overlapping genomic regions with exon / intron 
    / transript biotype annotations.

    tid2tio_dic:
        To check that if two features have same overlap amount with region,
        choose the more prominent transcript ID (i.e., better annotation).
        Also, for exonic sites to get distances from exon borders.
        Exonic sites == everything not "intron" or "intergenic".
    motif_hits:
        If True, the -a file is the motif hits BED file. In this case, the region ID
        has to be reconstructed from the BED region info.
        Format of reg_ids_dic key if motif_hits=True:
        "chr1:10-15(+)motif_id"
    reg_ids_dic:
        If set, compare genomic region IDs with IDs in dictionary. If region ID 
        from dictionary not in overlap_annotations_bed, set to "intergenic".

    """

    assert os.path.exists(overlap_annotations_bed), "file %s does not exist" %(overlap_annotations_bed)

    reg2maxol_dic = {}
    reg2annot_dic = {}
    annot_col = 9
    c_ol_nt_col = 12

    motif_id_pattern = re.compile(r"^.+?:(.+?);")

    with open(overlap_annotations_bed) as f:
        for line in f:
            cols = line.strip().split("\t")
            reg_id = cols[3]  # format: chr8:90314134-90314381(+) (if motif_hits=False)
            reg_strand = cols[5]
            reg_s = int(cols[1])
            reg_e = int(cols[2])

            # If motif_ids, construct new unique region_id from BED region info.
            if motif_hits:
                chr_id = cols[0]

                # col[3] has format: "rbp_id:motif_id;1;method_id:data_id". Extract motif_id from this string.
                # m = re.search("^.+?:(.+?);", reg_id)
                # assert m is not None, "Motif ID extraction failed for region ID \"%s\"" %(reg_id)
                # motif_id = m.group(1)
                m = motif_id_pattern.search(reg_id)
                assert m is not None, "Motif ID extraction failed for region ID \"%s\"" %(reg_id)
                motif_id = m.group(1)

                # motif_id = reg_id.split(":")[1].split(";")[0]
                reg_id = chr_id + ":" + str(reg_s+1) + "-" + str(reg_e) + "(" + reg_strand + ")" + motif_id
                annot_col = 13  # These shift since motif hits BED contains additional (4) p-value and score columns.
                c_ol_nt_col = 16

            annot_reg_s = int(cols[annot_col-2])
            annot_reg_e = int(cols[annot_col-1])

            annot_ids = cols[annot_col].split(";")
            assert len(annot_ids) == 3, "len(annot_ids) != 3 (expected ; separated string, but got: \"%s\")" %(cols[9])
            annot_id = annot_ids[0]
            tr_id = annot_ids[1]
            c_overlap_nts = int(cols[c_ol_nt_col])

            # Get distanced to intron borders (only for intron annotations).
            border_dist = -1
            us_ds_label = "-"
            # Get exon/intron number.
            exon_intron_nr_c = annot_ids[2]  # format: "exon_nr-exon_c" or "intron_nr-intron_c"
            # Get exon/intron number + total count for transcript.
            exon_intron_nr = int(exon_intron_nr_c.split("-")[0])
            exon_intron_c = int(exon_intron_nr_c.split("-")[1])

            # If not intronic site, we need to get exon start and end positions to calculate distance to border.
            if annot_id != "intron":
                c_exon_coords = len(tid2tio_dic[tr_id].exon_coords)
                assert exon_intron_c == c_exon_coords, "exon_intron_c != c_exon_coords for annotation \"%s\" (%i != %i)" %(annot_ids, exon_intron_c, c_exon_coords)
                annot_reg_s, annot_reg_e = tid2tio_dic[tr_id].exon_coords[exon_intron_nr-1]

            annot_reg_len = annot_reg_e - annot_reg_s

            # Get distance to closest annotation region border + if it is upstream (us) or downstream (ds).
            # if annot_id == "intron":
            border_dist, us_ds_label = get_dist_to_next_border(reg_s, reg_e, annot_reg_s, annot_reg_e, reg_strand,
                                                               skip_sites_not_within=True)  # Return False if site not within.

            # if border_dist == -1:
            #     print(reg_id, "reg:", reg_s, reg_e, "annot_reg:", annot_reg_s, annot_reg_e, "annot_id:", annot_id, "tr_id:", tr_id)

            if reg_id not in reg2maxol_dic:
                reg2maxol_dic[reg_id] = c_overlap_nts
                reg2annot_dic[reg_id] = [annot_id, tr_id, border_dist, us_ds_label, annot_reg_len, exon_intron_nr_c]
            else:
                stored_c_overlap_nts = reg2maxol_dic[reg_id]
                if c_overlap_nts > stored_c_overlap_nts:
                    reg2maxol_dic[reg_id] = c_overlap_nts
                    reg2annot_dic[reg_id][0] = annot_id
                    reg2annot_dic[reg_id][1] = tr_id
                    reg2annot_dic[reg_id][2] = border_dist
                    reg2annot_dic[reg_id][3] = us_ds_label
                    reg2annot_dic[reg_id][4] = annot_reg_len
                    reg2annot_dic[reg_id][5] = exon_intron_nr_c
                elif c_overlap_nts == stored_c_overlap_nts:  # if region overlaps with > 1 feature by same amount.
                    # if tid2tio_dic is not None:
                    best_tid = reg2annot_dic[reg_id][1]
                    new_best_tid = select_more_prominent_tid(tr_id, best_tid, tid2tio_dic)
                    # If current tr_id has better annotation, update region annotation.
                    if new_best_tid == tr_id:
                        reg2maxol_dic[reg_id] = c_overlap_nts
                        reg2annot_dic[reg_id][0] = annot_id
                        reg2annot_dic[reg_id][1] = new_best_tid
                        reg2annot_dic[reg_id][2] = border_dist
                        reg2annot_dic[reg_id][3] = us_ds_label
                        reg2annot_dic[reg_id][4] = annot_reg_len
                        reg2annot_dic[reg_id][5] = exon_intron_nr_c

            test_reg_id = "chr19:3980035-3980070(-)"
            if reg_id == test_reg_id:
                print("test_reg_id:", reg_id, "annot_id:", annot_id, "tr_id:", tr_id, "border_dist:", border_dist, "us_ds_label:", us_ds_label, "exon_intron_nr_c:", exon_intron_nr_c, "c_overlap_nts:", c_overlap_nts)
                print("set:", reg2annot_dic[reg_id])
    f.closed

    if reg_ids_dic is not None:
        for reg_id in reg_ids_dic:
            if reg_id not in reg2annot_dic:
                reg2annot_dic[reg_id] = ["intergenic", False, -1, "-", -1, "0-0"]

    return reg2annot_dic


"""
$ intersectBed -a regions.bed -b introns.bed -s -wo -f 0.5
chr1	1980	2020	s1	0	+	chr1	1000	2000	exon	0	+	20
chr1	1980	2020	s1	0	+	chr1	2000	3000	intron	0	+	20

"""



################################################################################

def get_mrna_region_annotations(overlap_annotations_bed,
                                reg_ids_dic=None):
    """
    Get mRNA region annotations from overlapping motif hit regions with 
    mRNA region annotations (i.e. 5'UTR, CDS, 3'UTR in transcript context).

    Format of motif hit regions BED file:
    ENST00000434296.2	1368	1375	HNRNPL,HNRNPL_3;1;method_id,data_id	0	+	10.6566	7.24e-05	-1.0	-1.0
    ENST00000434296.2	832	840	HNRNPL,HNRNPL_4;1;method_id,data_id	0	+	10.47	8.96e-05	-1.0	-1.0
    ENST00000434296.2	1432	1440	HNRNPL,HNRNPL_4;1;method_id,data_id	0	+	10.66	7.4e-05	-1.0	-1.0

    reg_ids_dic:
        If set, compare genomic region IDs with IDs in dictionary. If region ID 
        from dictionary not in overlap_annotations_bed, set label "intergenic".

    """

    assert os.path.exists(overlap_annotations_bed), "file %s does not exist" %(overlap_annotations_bed)

    reg2maxol_dic = {}
    reg2annot_dic = {}
    annot_col = 9
    c_ol_nt_col = 12

    motif_id_pattern = re.compile(r"^.+?:(.+?);")

    with open(overlap_annotations_bed) as f:
        for line in f:
            cols = line.strip().split("\t")
            tr_id = cols[0]
            motif_hit_s = int(cols[1])
            motif_hit_e = int(cols[2])
            motif_hit_id = cols[3]
            mrna_reg_s = cols[11]
            mrna_reg_e = cols[12]
            mrna_reg_id = cols[13]
            c_overlap_nts = int(cols[16])

            # col[3] has format: "rbp_id:motif_id;1;method_id:data_id". Extract motif_id from this string.
            # m = re.search("^.+?:(.+?);", reg_id)
            # assert m is not None, "Motif ID extraction failed for region ID \"%s\"" %(motif_hit_id)
            # motif_id = m.group(1)
            m = motif_id_pattern.search(reg_id)
            assert m is not None, "Motif ID extraction failed for region ID \"%s\"" %(reg_id)
            motif_id = m.group(1)

            # motif_id = reg_id.split(":")[1].split(";")[0]
            # motif_id = motif_hit_id.split(",")[1].split(";")[0]
            # mrna_reg_id has format: ENST00000434296;3utr. Extract region type from this string.
            mrna_reg_type = mrna_reg_id.split(";")[1]  # can be: 5'UTR, CDS, 3'UTR

            # ID that identifies motif hit.
            reg_s = str(motif_hit_s + 1)
            reg_e = str(motif_hit_e)
            reg_id = tr_id + ":" + reg_s + "-" + reg_e + "(+)," + motif_id

            if reg_id not in reg2maxol_dic:
                reg2maxol_dic[reg_id] = c_overlap_nts
                reg2annot_dic[reg_id] = [mrna_reg_type, tr_id]
            else:
                if c_overlap_nts > reg2maxol_dic[reg_id]:
                    reg2maxol_dic[reg_id] = c_overlap_nts
                    reg2annot_dic[reg_id][0] = mrna_reg_type
                    reg2annot_dic[reg_id][1] = tr_id
    f.closed

    if reg_ids_dic is not None:
        for reg_id in reg_ids_dic:
            if reg_id not in reg2annot_dic:
                reg2annot_dic[reg_id] = ["intergenic", False]

    return reg2annot_dic


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

    if re.search(r".+\.gz$", in_gtf):
        f = gzip.open(in_gtf, 'rt')
    else:
        f = open(in_gtf, "r")
    for line in f:
        # Skip header.
        if line.startswith("#"):
        # if re.search("^#", line):
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
        m = re.search(r'exon_number "(\d+?)"', infos)
        # Try GENCODE encoding.
        if not m:
            m = re.search(r'exon_number (\d+?);', infos)
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

def gtf_output_gene_regions_to_bed(in_gtf, out_bed,
                                   bed_col6_infos=1,
                                   gids_dic=False,
                                   chr_id_style=0):
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

    if re.search(r".+\.gz$", in_gtf):
        f = gzip.open(in_gtf, 'rt')
    else: 
        f = open(in_gtf, "r")
    for line in f:

        # Skip header.
        if line.startswith("#"):
            continue

        cols = line.strip().split("\t")
        feature = cols[2]
        if feature != "gene":
            continue

        chr_id = cols[0]
        feat_s = int(cols[3])  # 1-based index (see e.g. start_codon feature for proof).
        feat_e = int(cols[4])
        feat_pol = cols[6]
        infos = cols[8]

        chr_id = check_convert_chr_id(chr_id, id_style=chr_id_style)
        # If not one of standard chromosomes, continue.
        if not chr_id:
            continue

        assert feat_e >= feat_s, "feature end < feature start in GTF file \"%s\", line \"%s\". Since both coordinates are expected to have 1-based index, this should not happen" %(in_gtf, line)

        m = re.search(r'gene_id "(.+?)"', infos)
        assert m, "gene_id entry missing in GTF file \"%s\", line \"%s\"" %(in_gtf, line)
        gene_id = m.group(1)

        if gids_dic:
            if gene_id not in gids_dic:
                continue

        m = re.search(r'gene_name "(.+?)"', infos)
        gene_name = "-"  # optional.
        if m:
            gene_name = m.group(1)
        gene_biotype = "-"  # # optional.
        m = re.search(r'gene_biotype "(.+?)"', infos)
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

    return c_gene_regions


################################################################################

def select_more_prominent_tid(tid1, tid2, tid2tio_dic):
    """
    Given two transcript IDs tid1 and tid2, select the more prominent 
    transcript (i.e., with better experimental evidence).

    Features to compare:
    tidtio_dic[tid].mane_select  # 0 or 1
    tidtio_dic[tid].basic_tag  # 0 or 1
    tidtio_dic[tid].ensembl_canonical  # 0 or 1
    tidtio_dic[tid].tsl_id  # 1-5 or "NA"

    """
    assert tid1 in tid2tio_dic, "tid1 \"%s\" not in tid2tio_dic" %(tid1)
    assert tid2 in tid2tio_dic, "tid2 \"%s\" not in tid2tio_dic" %(tid2)

    # Comparison dictionary.
    id2sc = {}
    for i in range(5):
        pos = i + 1
        pos_str = "%i" %(pos)
        id2sc[pos_str] = pos
    id2sc["NA"] = 6

    if tid1 == tid2:
        return tid1
    if tid2tio_dic[tid1].mane_select > tid2tio_dic[tid2].mane_select:
        return tid1
    if tid2tio_dic[tid1].mane_select < tid2tio_dic[tid2].mane_select:
        return tid2
    if tid2tio_dic[tid1].basic_tag > tid2tio_dic[tid2].basic_tag:
        return tid1
    if tid2tio_dic[tid1].basic_tag < tid2tio_dic[tid2].basic_tag:
        return tid2
    if tid2tio_dic[tid1].ensembl_canonical > tid2tio_dic[tid2].ensembl_canonical:
        return tid1
    if tid2tio_dic[tid1].ensembl_canonical < tid2tio_dic[tid2].ensembl_canonical:
        return tid2
    if id2sc[tid2tio_dic[tid1].tsl_id] < id2sc[tid2tio_dic[tid2].tsl_id]:
        return tid1
    if id2sc[tid2tio_dic[tid1].tsl_id] > id2sc[tid2tio_dic[tid2].tsl_id]:
        return tid2
    # Unspliced transcript lengths.
    tid1_len = tid2tio_dic[tid1].tr_e - tid2tio_dic[tid1].tr_s + 1
    tid2_len = tid2tio_dic[tid2].tr_e - tid2tio_dic[tid2].tr_s + 1
    if tid1_len > tid2_len:
        return tid1
    if tid1_len < tid2_len:
        return tid2
    return tid1


################################################################################

def select_mpts_from_gene_infos(gid2gio_dic,
                                basic_tag=True,
                                ensembl_canonical_tag=False,
                                mane_select_tag=False,
                                only_tsl=False,
                                prior_basic_tag=True,
                                prior_mane_select=False,
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
    mane_select_tag:
        If True only report transcripts with "MANE_Select" tag.
    only_tsl:
        If True only report transcripts with TSL 1-5 (excluding "NA").
    tr_min_len:
        If length set, only report transcripts with length >= tr_min_len.
    prior_mane_select:
        If True, MANE_Select tag trumps all other tags. According to manual,
        the MANE select is a default transcript per human gene, present in RefSeq
        and Ensembl databases.
        
    >>> test_gtf = "test_data/test_mpt_selection.gtf"
    >>> gid2gio_dic = gtf_read_in_gene_infos(test_gtf)
    >>> tr_ids_dic = select_mpts_from_gene_infos(gid2gio_dic, basic_tag=False, ensembl_canonical_tag=False, only_tsl=False)
    >>> tr_ids_dic
    {'ENST00000357266': 'ENSG00000096060'}
    >>> test_gtf = "test_data/test_mpt_selection2.gtf"
    >>> gid2gio_dic = gtf_read_in_gene_infos(test_gtf)
    >>> tr_ids_dic = select_mpts_from_gene_infos(gid2gio_dic, basic_tag=False, ensembl_canonical_tag=False, only_tsl=False, prior_basic_tag=True)
    >>> tr_ids_dic
    {'ENST00000530167': 'ENSG00000188486'}
    >>> test_gtf = "test_data/test_mpt_selection3.gtf"
    >>> gid2gio_dic = gtf_read_in_gene_infos(test_gtf)
    >>> tr_ids_dic = select_mpts_from_gene_infos(gid2gio_dic, basic_tag=False, ensembl_canonical_tag=False, only_tsl=False, prior_basic_tag=True)
    >>> tr_ids_dic
    {'g1_t2': 'g1', 'g2_t1': 'g2'}

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
        mpt_ms = 0

        for idx, tr_id in enumerate(gene_info.tr_ids):
            # print("mpt_id:", mpt_id, "tr_id:", tr_id)
            tr_tsl = gene_info.tr_tsls[idx]  # 1-5 or NA
            tr_bt = gene_info.tr_basic_tags[idx]  # 0 or 1
            tr_ec = gene_info.tr_ensembl_canonical_tags[idx]  # 0 or 1
            tr_ms = gene_info.tr_mane_select_tags[idx]  # 0 or 1
            tr_length = gene_info.tr_lengths[idx]

            # print(tr_id, "BT:", tr_bt, "TSL:", tr_tsl)
            # print(tr_id, tr_bt, tr_ec, tr_length, tr_tsl)
            if basic_tag:
                if not tr_bt:
                    continue
            if ensembl_canonical_tag:
                if not tr_ec:
                    continue
            if mane_select_tag:
                if not tr_ms:
                    continue
            if only_tsl:
                if tr_tsl == "NA":
                    continue
            if tr_min_len:
                if tr_length < tr_min_len:
                    continue
            if prior_mane_select:
                if tr_ms > mpt_ms:
                    mpt_id = tr_id
                    mpt_tsl = tr_tsl
                    mpt_len = tr_length
                    mpt_bt = tr_bt
                    mpt_ec = tr_ec
                    mpt_ms = tr_ms
                    continue

            if id2sc[tr_tsl] < id2sc[mpt_tsl]:
                if prior_basic_tag:
                    if tr_bt >= mpt_bt:
                        mpt_id = tr_id
                        mpt_tsl = tr_tsl
                        mpt_len = tr_length
                        mpt_bt = tr_bt
                        mpt_ec = tr_ec
                        mpt_ms = tr_ms
                else:
                    mpt_id = tr_id
                    mpt_tsl = tr_tsl
                    mpt_len = tr_length
                    mpt_bt = tr_bt
                    mpt_ec = tr_ec
                    mpt_ms = tr_ms

            elif id2sc[tr_tsl] == id2sc[mpt_tsl]:
                # print("Now equal, comparing tr_id %s with mpt_id %s" %(tr_id, mpt_id))
                # If transcript has basic tag, use this.
                if tr_bt > mpt_bt:
                    mpt_id = tr_id
                    mpt_tsl = tr_tsl
                    mpt_len = tr_length
                    mpt_bt = tr_bt
                    mpt_ec = tr_ec
                    mpt_ms = tr_ms
                    continue
                # If transcript has Ensembl canonical tag, use this.
                if tr_ec > mpt_ec:
                    mpt_id = tr_id
                    mpt_tsl = tr_tsl
                    mpt_len = tr_length
                    mpt_bt = tr_bt
                    mpt_ec = tr_ec
                    mpt_ms = tr_ms
                    continue
                # If same basic/Ensembl canonical tag combination.
                if tr_ec == mpt_ec and tr_bt == mpt_bt:
                    if tr_length > mpt_len:
                        mpt_id = tr_id
                        mpt_tsl = tr_tsl
                        mpt_len = tr_length
                        mpt_bt = tr_bt
                        mpt_ec = tr_ec
                        mpt_ms = tr_ms
            else:
                # If transcript has worse TSL but basic tag and current MPT has not.
                if prior_basic_tag and tr_bt > mpt_bt:
                    mpt_id = tr_id
                    mpt_tsl = tr_tsl
                    mpt_len = tr_length
                    mpt_bt = tr_bt
                    mpt_ec = tr_ec
                    mpt_ms = tr_ms

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

def get_fid2desc_mapping(fid2desc_file):
    """
    Get RBP molecular function to description mapping.
    
    """

    assert os.path.exists(fid2desc_file), "file %s does not exist" %(fid2desc_file)

    # Get function ID to function description mapping.
    fid2desc_dic = {}
    desc2fid_dic = {}

    with open(fid2desc_file) as f:
        for line in f:
            if line.startswith('Function ID'):
                continue
            cols = line.strip().split("\t")
            fid = cols[0]
            disc = cols[1]
            fid2desc_dic[fid] = disc
            desc2fid_dic[disc] = fid
    f.closed

    return fid2desc_dic, desc2fid_dic


################################################################################

def get_fid_db_counts(name2ids_dic, name2fids_dic,
                      motif_level=False):
    """
    Get function ID database counts (i.e. number of their appearances 
    over all RBPs, so count appearance for every motif of an RBP!).

    >>> name2ids_dic = {'A1CF': ['A1CF_1', 'A1CF_2'], 'ACIN1': ['ACIN1_1']}
    >>> name2fids_dic = {'A1CF': ['RM', 'RSD', 'RE'], 'ACIN1': ['RSD']}
    >>> get_fid_db_counts(name2ids_dic, name2fids_dic, motif_level=True)
    {'RM': 2, 'RSD': 3, 'RE': 2}
    >>> get_fid_db_counts(name2ids_dic, name2fids_dic, motif_level=False)
    {'RM': 1, 'RSD': 2, 'RE': 1}
    
    """
    fid2dbc_dic = {}

    for rbp_name in name2ids_dic:
        if motif_level:
            for motif_id in name2ids_dic[rbp_name]:
                fids_list = name2fids_dic[rbp_name]
                for fid in fids_list:
                    if fid in fid2dbc_dic:
                        fid2dbc_dic[fid] += 1
                    else:
                        fid2dbc_dic[fid] = 1
        else:
            fids_list = name2fids_dic[rbp_name]
            for fid in fids_list:
                if fid in fid2dbc_dic:
                    fid2dbc_dic[fid] += 1
                else:
                    fid2dbc_dic[fid] = 1

    return fid2dbc_dic


################################################################################

def get_rbp_id_mappings(rbp2ids_file,
                        only_meme_xml=False):
    """
    Read in file mapping RBP names to motif IDs and motif types.
    Return dictionaries with:
    rbp_name -> motif_ids_list
    motif_id -> motif_type

    FILE FORMAT:

    RBP_motif_ID	RBP_name	Motif_type	Organism
    AGGF1_1	AGGF1	meme_xml	human
    AGGF1_2	AGGF1	meme_xml	human
    AKAP1_1	AKAP1	meme_xml	human
    BCCIP_1	BCCIP	meme_xml	human
    BUD13_1	BUD13	meme_xml	human
    ...
    RF00032	SLBP	cm	human

    RBPBench v0.3:
    Currently ignore Organism column / do not use this information.

    RBPBench v1.0 updated:
    RBP_motif_ID	RBP_name	Motif_type	Organism	Gene_ID	Function_IDs	Reference	Experiment	Comments
    A1CF_1	A1CF	meme_xml	human	ENSG00000148584	RM;RSD;RE	34086933	-	-
    A1CF_2	A1CF	meme_xml	human	ENSG00000148584	RM;RSD;RE	34086933	-	-
    ACIN1_1	ACIN1	meme_xml	human	-	-	34086933	-	-

    RBPBench v1.01 updated (more pubmed IDs + experiment infos):
    RBP_motif_ID	RBP_name	Motif_type	Organism	Gene_ID	Function_IDs	Reference	Experiment	Comments
    A1CF_1	A1CF	meme_xml	human	ENSG00000148584	RM;RSD;RE	31724725;10669759	RBNS_ENCODE;RBPDB	-
    A1CF_2	A1CF	meme_xml	human	ENSG00000148584	RM;RSD;RE	31724725	RBNS_ENCODE	-
    ACIN1_1	ACIN1	meme_xml	human	ENSG00000100813	-	27365209	iCLIP	-
    ACIN1_2	ACIN1	meme_xml	human	ENSG00000100813	-	27365209	iCLIP	-
        
    """
    name2ids_dic = {}
    name2gid_dic = {}
    id2type_dic = {}
    name2fids_dic = {}
    id2pids_dic = {}
    id2exp_dic = {}
    # id2org_dic = {}

    with open(rbp2ids_file) as f:
        for line in f:
            if line.startswith("RBP_motif_ID") or line.startswith("#"):
                continue

            cols = line.strip().split("\t")
            motif_id = cols[0]
            rbp_name = cols[1]
            motif_type = cols[2]

            id2type_dic[motif_id] = motif_type
            if only_meme_xml:
                if motif_type != "meme_xml":
                    continue

            if rbp_name in name2ids_dic:
                name2ids_dic[rbp_name].append(motif_id)
            else:
                name2ids_dic[rbp_name] = [motif_id]

            name2gid_dic[rbp_name] = "-"
            name2fids_dic[rbp_name] = []

            if len(cols) > 3:

                # organism = cols[3]
                # if organism != "human":
                #     rbp_name = rbp_name + "_" + organism
                # id2org_dic[motif_id] = organism
                gene_id = cols[4]
                function_ids = cols[5]
                pids = cols[6]
                exp = cols[7]

                fids_list = []
                if function_ids != "-":
                    fids_list = function_ids.split(";")
        
                fids_list.sort()

                name2gid_dic[rbp_name] = gene_id
                name2fids_dic[rbp_name] = fids_list

                # Split pids.
                pids_list = pids.split(';')
                exp_list = exp.split(';')

                id2pids_dic[motif_id] = pids_list
                id2exp_dic[motif_id] = exp_list

    f.closed

    assert name2ids_dic, "no RBP IDs read in from %s" %(rbp2ids_file)

    return name2ids_dic, id2type_dic, name2gid_dic, name2fids_dic, id2pids_dic, id2exp_dic


################################################################################

def get_uniq_gen_size(gen_sites_bed):
    """
    Get unique genomic space size, which the genomic sites inside
    gen_sites_bed cover.

    sort -k1,1 -k2,2n in_sites.filtered.bed | mergeBed -i stdin -s -c 4 -o distinct -delim ";"

    >>> gen_sites_bed = "test_data/test_gen_size.bed"
    >>> get_uniq_gen_size(gen_sites_bed)
    2500
    >>> in_bed = "test_data/test_ol.sorted.bed"    
    >>> get_uniq_gen_size(in_bed)
    2100

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

def bed_check_format(bed_file, asserts=True,
                     param_str=False):
    """
    Check whether given BED file is not empty + has >= 6 columns.

    >>> bed_file = "test_data/test2.bed"
    >>> bed_check_format(bed_file)
    True
    
    """

    pols = ["+", "-"]
    okay = False
    with open(bed_file) as f:
        for line in f:
            cols = line.strip().split("\t")
            if asserts:
                assert len(cols) >= 6, "invalid --in BED format. Please provide valid BED file (i.e., >= 6 column format)"
                reg_pol = cols[5]
                assert reg_pol in pols, "invalid polarity in --in BED column 6 (found: \"%s\", expected \"+\" or \"-\"). Please provide valid BED file (i.e., >= 6 column format), or use --unstranded option" %(reg_pol)
                okay = True
                break
            else:
                if len(cols) >= 6:
                    reg_pol = cols[5]
                    if reg_pol in pols:
                        okay = True
                        break
    f.closed
    if asserts:
        false_msg = "provided file \"%s\" not in BED format (file empty or not in 6-column format?)" %(bed_file)
        if param_str:
            false_msg = "provided file \"%s\" (via %s) not in BED format (file empty or not in 6-column format?)" %(bed_file, param_str)
        assert okay, false_msg
    return okay


################################################################################

def fasta_check_format(fasta_file):
    """
    Quick check if file is a valid FASTA file.

    >>> test_fa = "test_data/test.fa"
    >>> fasta_check_format(test_fa)
    True
    >>> test_fa = "test_data/test.bed"
    >>> fasta_check_format(test_fa)
    False

    """
    with open(fasta_file, 'r') as f:
        lines = f.readlines()
    if len(lines) == 0:
        return False
    if not lines[0].startswith(">"):
        return False
    
    return True


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
            elif len(cols) == 21:
                if cols[0] == "data_id" and cols[20] == "internal_id":
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

    if re.search(r"\w+:\d+-\d+\([+|-]\)", reg_id):
        m = re.search(r"(\w+):(\d+)-(\d+)\(", reg_id)
        core_id = "%s:%s-%s" %(m.group(1), m.group(2), m.group(3))
        return core_id
    else:
        assert False, "region ID has invalid format (%s). Please contact developers" %(reg_id)


################################################################################

def get_hit_id_elements(hit_id):
    """
    From  hit ID to ID elements list.

    Hit ID format:
    chr1:100-110(+)motif_id

    >>> hit_id = "chr1:100-110(+)motif_id"
    >>> get_hit_id_elements(hit_id)
    ['chr1', '100', '110', '+', 'motif_id']
    
    """

    if re.search(r"^\w+?:\d+-\d+\([+|-]\)\w+", hit_id):
        m = re.search(r"^(\w+?):(\d+)-(\d+)\(([+|-])\)(.+)", hit_id)
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

def bed_read_chr_ids_dic(in_bed):
    """
    Read in chromosome IDs from BED file.
    Mapping is chromosome ID -> row count.

    """

    chr_ids_dic = {}

    with open(in_bed) as f:
        for line in f:
            cols = line.strip().split("\t")
            chr_id = cols[0]
            if chr_id in chr_ids_dic:
                chr_ids_dic[chr_id] += 1
            else:
                chr_ids_dic[chr_id] = 1

    f.closed

    return chr_ids_dic


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
        if re.search(r"\w+:\d+-\d+\([+|-]\)", seq_id):
            m = re.search(r"\w+:\d+-\d+\(([+|-])\)", seq_id)
            reg2pol_dic[seq_id] = m.group(1)
        else:
            assert False, "region ID has invalid format (%s). Please contact developers" %(seq_id)
    return reg2pol_dic


################################################################################

def bed_check_ids_output_bed(in_bed, out_bed,
                             id_check=True,
                             new_header_id="reg",
                             make_uniq_headers=False):
    """
    Given in_bed BED file, check column 4 IDs, 
    output to out_bed with unique IDs. 
    Return regions in dictionary.
    
    """

    OUTBED = open(out_bed, "w")

    c_reg = 0

    bed_reg_dic = {}

    with open(in_bed) as f:
        for line in f:
            cols = line.strip().split("\t")
            chr_id = cols[0]
            reg_s = cols[1]
            reg_e = cols[2]
            reg_id = cols[3]
            reg_sc = cols[4]
            reg_pol = cols[5]

            c_reg += 1
            if make_uniq_headers:
                reg_id = new_header_id + "_" + str(c_reg)

            if id_check:
                assert reg_id not in bed_reg_dic, "non-unique region ID \"%s\" found in --in BED file. Please provide unique column 4 IDs or set --make-uniq-headers" %(reg_id)
            
            bed_reg_dic[reg_id] = [chr_id, reg_s, reg_e, reg_pol]

            OUTBED.write("%s\t%s\t%s\t%s\t%s\t%s\n" % (chr_id, reg_s, reg_e, reg_id, reg_sc, reg_pol))

    f.closed
    OUTBED.close()

    return bed_reg_dic


################################################################################

def merge_files(files_list, out_file):
    """
    Merge list of files into one output file.

    """
    assert files_list, "given files_list is empty"
    # Delete out_file if exists.
    if os.path.exists(out_file):
        os.remove(out_file)
    for f in files_list:
        assert os.path.isfile(f), "list file \"%s\" not found" % (f)
        assert f != out_file, "cat does not like to cat file into same file (%s)" %(check_cmd)
        check_cmd = "cat " + f + " >> " + out_file
        output = subprocess.getoutput(check_cmd)
        error = False
        if output:
            error = True
        assert error == False, "cat did not like your input (in_file: %s, out_file: %s):\n%s" %(f, out_file, output)


################################################################################

def bed_generate_random_negatives(in_bed, chr_sizes_file, out_bed,
                                  incl_bed=False,
                                  excl_bed=False,
                                  allow_overlaps=False,
                                  seed=None):
    """
    Shuffle given in_bed, generating random negative regions. Optionally,
    the regions to extract negatives from can be controlled by incl_bed
    and excl_bed.

    in_bed:
        .bed file containing regions to shuffle, i.e., generate same number
        of random negatives (with same size distribution too)
    chr_sizes_file:
        File that stores chromosome IDs and their sizes
    out_bed:
        Output random negative regions in out_bed
    incl_bed:
        Regions from which to extract random negatives
    excl_bed:
        Regions from which no random negatives should be extracted
    allow_overlaps:
        Allow random negatives to overlap with each other

    Returns:
    Function returns True if no error occured.
    If loci error occured, function returns False.
    Any other error will throw an assertion error.
    If it is not possible to get the number of random negatives with the given
    restrictions, bedtools shuffle will throw the following error:
    Error, line 3: tried 1000 potential loci for entry, but could not avoid
    excluded regions.  Ignoring entry and moving on.
    This error will be thrown for every failed attempt to find a random
    negative for a certain positive instance.


    Tool:    bedtools shuffle (aka shuffleBed)
    Version: v2.31.1
    Summary: Randomly permute the locations of a feature file among a genome.

    Usage:   bedtools shuffle [OPTIONS] -i <bed/gff/vcf> -g <genome>

    Options: 
        -excl	A BED/GFF/VCF file of coordinates in which features in -i
            should not be placed (e.g. gaps.bed).

        -incl	Instead of randomly placing features in a genome, the -incl
            options defines a BED/GFF/VCF file of coordinates in which 
            features in -i should be randomly placed (e.g. genes.bed). 
            Larger -incl intervals will contain more shuffled regions. 
            This method DISABLES -chromFirst. 
        -chrom	Keep features in -i on the same chromosome.
            - By default, the chrom and position are randomly chosen.
            - NOTE: Forces use of -chromFirst (see below).

        -seed	Supply an integer seed for the shuffling.
            - By default, the seed is chosen automatically.
            - (INTEGER)

        -f	Maximum overlap (as a fraction of the -i feature) with an -excl
            feature that is tolerated before searching for a new, 
            randomized locus. For example, -f 0.10 allows up to 10%
            of a randomized feature to overlap with a given feature
            in the -excl file. **Cannot be used with -incl file.**
            - Default is 1E-9 (i.e., 1bp).
            - FLOAT (e.g. 0.50)

        -chromFirst	
            Instead of choosing a position randomly among the entire
            genome (the default), first choose a chrom randomly, and then
            choose a random start coordinate on that chrom.  This leads
            to features being ~uniformly distributed among the chroms,
            as opposed to features being distribute as a function of chrom size.

        -bedpe	Indicate that the A file is in BEDPE format.

        -maxTries	
            Max. number of attempts to find a home for a shuffled interval
            in the presence of -incl or -excl.
            Default = 1000.
        -noOverlapping	
            Don't allow shuffled intervals to overlap.
        -allowBeyondChromEnd	
            Allow shuffled intervals to be relocated to a position
            in which the entire original interval cannot fit w/o exceeding
            the end of the chromosome.  In this case, the end coordinate of the
            shuffled interval will be set to the chromosome's length.
            By default, an interval's original length must be fully-contained
            within the chromosome.

    """
    # Check for bedtools.
    assert is_tool("bedtools"), "bedtools not in PATH"
    # Construct call.
    check_cmd = "bedtools shuffle "
    if excl_bed:
        check_cmd = check_cmd + "-excl " + excl_bed + " "
    if incl_bed:
        check_cmd = check_cmd + "-incl " + incl_bed + " "
    if not allow_overlaps:
        check_cmd = check_cmd + "-noOverlapping "
    if seed is not None:
        check_cmd = check_cmd + "-seed " + str(seed) + " "
    check_cmd = check_cmd + "-i " + in_bed + " -g " + chr_sizes_file + " > " + out_bed
    output = subprocess.getoutput(check_cmd)
    error = False
    if output:
        error = True
    # Look for "tried 1000 potential loci" error.
    if error:
        if re.search("potential loci", output):
            print("WARNING: number of extracted random negatives < requested number")
            return False
        else:
            assert False, "bedtools shuffle is complaining:\n%s\n%s" %(check_cmd, output)
    else:
        return True


################################################################################

def output_chromosome_lengths_file(len_dic, out_file,
                                   ids2print_dic=None):
    """
    Output chromosome lengths file with format:
    sequence_ID<tab>sequence_length
    
    """
    LOUT = open(out_file, "w")
    c_pr = 0
    for seq_id in len_dic:
        if ids2print_dic is not None:
            if seq_id in ids2print_dic:
                c_pr += 1
                LOUT.write("%s\t%i\n" %(seq_id, len_dic[seq_id]))
        else:
            c_pr += 1
            LOUT.write("%s\t%i\n" %(seq_id, len_dic[seq_id]))
    LOUT.close()
    assert c_pr, "nothing was printed out"


################################################################################

def genome_fasta_get_chr_sizes_file(in_genome_fa, out_chr_sizes_file,
                                    check_ids=True,
                                    seq_len_dic=None):
    """
    Extract sequence names and lengths in format:
    seq_id<tab>seq_len

    >>> test_fa = "test_data/test.fa"
    >>> exp_chr_len_file = "test_data/test.fa.chr_len.exp.txt"
    >>> out_chr_len_file = "test_data/test.fa.chr_len.tmp.txt"
    >>> genome_fasta_get_chr_sizes_file(test_fa, out_chr_len_file)
    >>> diff_two_files_identical(out_chr_len_file, exp_chr_len_file)
    True

    """

    seq_id = "id"
    seq_len = 0
    seen_ids_dic = {}

    OUTCHRLEN = open(out_chr_sizes_file, "w")

    header_pattern = re.compile(r"^>(\S+)")

    with open(in_genome_fa) as f:
        for line in f:
            if line.startswith(">"):
                m = header_pattern.match(line)
                new_id = m.group(1)
                if check_ids:
                    assert new_id not in seen_ids_dic, f'non-unique sequence ID "{new_id}" found in --in genome FASTA file. Please provide unique sequence IDs'
                seen_ids_dic[new_id] = 1
                if seq_len:
                    OUTCHRLEN.write(f"{seq_id}\t{seq_len}\n")
                    if seq_len_dic is not None:
                        seq_len_dic[seq_id] = seq_len
                seq_len = 0
                seq_id = new_id
            else:
                seq_len += len(line.strip())

    # Print last sequence length.
    if seq_len:
        OUTCHRLEN.write(f"{seq_id}\t{seq_len}\n")
        if seq_len_dic is not None:
            seq_len_dic[seq_id] = seq_len

    OUTCHRLEN.close()


################################################################################

def genome_fasta_get_chr_sizes(in_genome_fa,
                               check_ids=True):
    """
    Extract sequence names and lengths and return dictionary with
    seq_id -> seq_len mapping.

    >>> test_fa = "test_data/test.fa"
    >>> genome_fasta_get_chr_sizes(test_fa)
    {'seq1': 12, 'seq2': 20}
    
    """

    seq_len_dic = {}
    seq_id = "id"
    seq_len = 0

    header_pattern = re.compile(r"^>(\S+)")

    with open(in_genome_fa) as f:
        for line in f:
            if line.startswith(">"):
                m = header_pattern.match(line)
                new_id = m.group(1)
                if check_ids:
                    assert new_id not in seq_len_dic, "non-unique sequence ID \"%s\" found in --in genome FASTA file. Please provide unique sequence IDs" %(new_id)
                if seq_len:
                    seq_len_dic[seq_id] = seq_len
                seq_len = 0
                seq_id = new_id
            else:
                seq_len += len(line.strip())
    f.closed

    # Get last sequence length.
    if seq_len:
        seq_len_dic[seq_id] = seq_len

    return seq_len_dic


################################################################################

def merge_pos_bed_files(shuffle_list, bg_shuffle_in_bed,
                        core_neg_id="neg"):
    """
    Merge BED files in shuffle list and output to bg_shuffle_in_bed.
    Use core_neg_id to add to column 4 of negative regions.
    So new col4 ID format: pos1;neg1 pos1;neg2 ...

    """
    assert shuffle_list, "no BED files in shuffle_list"
    
    OUTBED = open(bg_shuffle_in_bed, "w")

    idx = 0
    for out_bed in shuffle_list:
        idx += 1
        with open(out_bed) as f:
            for line in f:
                cols = line.strip().split("\t")
                chr_id = cols[0]
                reg_s = cols[1]
                reg_e = cols[2]
                reg_id = cols[3]
                reg_sc = cols[4]
                reg_pol = cols[5]

                new_reg_id = reg_id + ";" + core_neg_id + str(idx)

                OUTBED.write("%s\t%s\t%s\t%s\t%s\t%s\n" % (chr_id, reg_s, reg_e, new_reg_id, reg_sc, reg_pol))

        f.closed
    
    OUTBED.close()


################################################################################

def bed_filter_extend_bed(in_bed, out_bed,
                          ext_up=0,
                          ext_down=0,
                          remove_dupl=True,
                          reg2sc_dic=None,
                          reg2pol_dic=None,
                          score_col=5,
                          score_thr=None,
                          score_rev_filter=False,
                          chr_ids_dic=None,
                          bed_chr_ids_dic=None,
                          use_region_ids=False,
                          chr_len_dic=False,
                          new_reg_ids=False,
                          core_reg_id="pos",
                          core_reg_dic=None,
                          transcript_sites=False,
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
    reg2sc_dic:
        region ID -> polarity/strand of region
    transcript_sites:
        If transcript sites, just use first three columns.
    chr_len_dic:
        Provide transcript lengths for each chromosome, to extend not beyond ends.
    core_reg_dic:
        Store core region infos for each region ID (if core_reg_dic not None).
        
    """

    OUTBED = open(out_bed, "w")

    assert score_col > 0, "invalid score column (< 1)"

    c_in = 0
    c_out = 0
    c_chr_filter = 0
    c_dupl_filter = 0
    c_sc_thr = 0

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

            # If input regions are transcript sites, we just need first three columns.
            if transcript_sites:
                reg_id = "%s:%i-%i" %(chr_id, reg_s, reg_e)
                reg_sc = float(cols[score_col-1])
                reg_pol = "+"

            c_in += 1

            # Check chr_id.
            if chr_ids_dic is not None:
                if chr_id not in chr_ids_dic:
                    c_chr_filter += 1
                    continue

            # Filter by site score.
            if score_thr is not None:
                if score_rev_filter:
                    if reg_sc > score_thr:
                        c_sc_thr += 1
                        continue
                else:
                    if reg_sc < score_thr:
                        c_sc_thr += 1
                        continue

            if new_reg_ids:
                reg_id = core_reg_id + str(c_in)

            if core_reg_dic is not None:
                core_reg_dic[reg_id] = [chr_id, reg_s, reg_e, reg_pol]

            # Record present IDs.
            if bed_chr_ids_dic is not None:
                bed_chr_ids_dic[chr_id] = 1

            if reg2pol_dic is not None:
                reg2pol_dic[reg_id] = reg_pol

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
                new_e = reg_e + ext_down
                # Bound checks.
                if new_s < 0:
                    new_s = 0
                if chr_len_dic:
                    if chr_id in chr_len_dic:
                        chr_len = chr_len_dic[chr_id]
                        if new_e > chr_len:
                            new_e = chr_len

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
                else:
                    if reg_id not in reg2sc_dic:
                        reg2sc_dic[reg_id] = reg_sc

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
                # Bound checks.
                if new_s < 0:
                    new_s = 0
                if chr_len_dic:
                    if chr_id in chr_len_dic:
                        chr_len = chr_len_dic[chr_id]
                        if new_e > chr_len:
                            new_e = chr_len

                if remove_dupl:
                    reg_str = "%s:%i-%i(%s)" % (chr_id, new_s, new_e, reg_pol)
                    if reg_str not in reg2sc_dic:
                        reg2sc_dic[reg_str] = reg_sc
                    else:
                        c_dupl_filter += 1
                        continue
                else:
                    if reg_id not in reg2sc_dic:
                        reg2sc_dic[reg_id] = reg_sc

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
    stats_dic["c_sc_thr"] = c_sc_thr

    return stats_dic


################################################################################

def seq_id_get_parts(seq_id):
    """
    Given sequence ID with format chr20:62139082-62139128(-),
    return single parts.

    """

    if re.search(r"\w+:\d+-\d+\([+|-]\)", seq_id):
        m = re.search(r"(\w+):(\d+)-(\d+)\(([+|-])\)", seq_id)
        chr_id = m.group(1)
        reg_s = int(m.group(2))
        reg_e = int(m.group(3))
        strand = m.group(4)

        return [chr_id, reg_s, reg_e, strand]

    else:
        assert False, "invalid sequence ID format (%s). Expected format: chr6:666-6666(-)" %(seq_id)


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

    if re.search(r"\w+:\d+-\d+\([+|-]\)", seq_name):
        m = re.search(r"(\w+):(\d+)-(\d+)\(([+|-])\)", seq_name)
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
        assert False, "invalid seq_name format given (%s). This error could be due to --meme-no-check set and MEME >= v5.5.4 installed. In this case, please set --meme-no-pgc" %(seq_name)


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
    if re.search(r"\w+:\d+-\d+\(", seq_name):
        m = re.search(r"\w+:(\d+)-(\d+)\(", seq_name)
        reg_s = int(m.group(1))
        reg_e = int(m.group(2))
        return reg_e - reg_s
    else:
        assert False, "invalid seq_name format given (%s)" %(seq_name)


################################################################################

class EnmoStats:
    """
    Store motif enrichment stats for each motif.
    
    """

    def __init__(self,
                 motif_id: str,
                 rbp_id: str,
                 c_pos_hit_regions = 0,
                 c_neg_hit_regions = 0,
                 c_pos_regions = 0,
                 c_neg_regions = 0,
                 c_pos_hits = 0,
                 c_neg_hits = 0,
                 con_table = False,  # Continency table.
                 fisher_pval = 1.0,
                 fisher_pval_corr = 1.0,
                 fisher_corr_mode = 1,  # 1: BH, 2: Bonferroni, 3: no correction
                 fisher_alt_hyp_mode = 1,  # Alternative hypothesis mode, 1: greater, 2: two-sided, 3: less
                 motif_type="meme_xml",
                 consensus_seq="-",  # Consensus sequence of sequence motif (for structure motif "-", for regex "regex_string").
                 logo_png_file = False) -> None:
        self.motif_id = motif_id
        self.rbp_id = rbp_id
        self.c_pos_hit_regions = c_pos_hit_regions
        self.c_neg_hit_regions = c_neg_hit_regions
        self.c_pos_regions = c_pos_regions
        self.c_neg_regions = c_neg_regions
        self.c_pos_hits = c_pos_hits
        self.c_neg_hits = c_neg_hits
        self.con_table = con_table
        self.fisher_pval = fisher_pval
        self.fisher_pval_corr = fisher_pval_corr
        self.fisher_corr_mode = fisher_corr_mode
        self.fisher_alt_hyp_mode = fisher_alt_hyp_mode
        self.motif_type = motif_type
        self.consensus_seq = consensus_seq
        self.logo_png_file = logo_png_file


################################################################################

class NemoStats:
    """
    Store motif enrichment stats for each motif.
    
    """

    def __init__(self,
                 motif_id: str,
                 rbp_id: str,
                 c_pos_hit_regions = 0,
                 c_neg_hit_regions = 0,
                 c_pos_regions = 0,
                 c_neg_regions = 0,
                 c_pos_hits = 0,
                 c_neg_hits = 0,
                 con_table = False,  # Continency table.
                 fisher_pval = 1.0,
                 fisher_pval_corr = 1.0,
                 fisher_corr_mode = 1,  # 1: BH, 2: Bonferroni, 3: no correction
                 fisher_alt_hyp_mode = 1,  # Alternative hypothesis mode, 1: greater, 2: two-sided, 3: less
                 motif_type="meme_xml",
                 consensus_seq="-",  # Consensus sequence of sequence motif (for structure motif "-", for regex "regex_string").
                 pos_set_avg_center_dist="-",
                 neg_set_avg_center_dist="-",
                 pos_set_max_center_dist="-",
                 pos_set_max_center_dist_c="-",
                 neg_set_max_center_dist="-",
                 neg_set_max_center_dist_c="-",
                 logo_png_file = False,
                 wrs_pval_two_sided=False,
                 wrs_test_stat_two_sided=False,
                 wrs_pval_greater=False,
                 wrs_test_stat_greater=False,
                 wrs_pval_less=False,
                 wrs_test_stat_less=False,         
                 dist_plot_counts_dic={}) -> None:
        self.motif_id = motif_id
        self.rbp_id = rbp_id
        self.c_pos_hit_regions = c_pos_hit_regions
        self.c_neg_hit_regions = c_neg_hit_regions
        self.c_pos_regions = c_pos_regions
        self.c_neg_regions = c_neg_regions
        self.c_pos_hits = c_pos_hits
        self.c_neg_hits = c_neg_hits
        self.con_table = con_table
        self.fisher_pval = fisher_pval
        self.fisher_pval_corr = fisher_pval_corr
        self.fisher_corr_mode = fisher_corr_mode
        self.fisher_alt_hyp_mode = fisher_alt_hyp_mode
        self.motif_type = motif_type
        self.consensus_seq = consensus_seq
        self.pos_set_avg_center_dist = pos_set_avg_center_dist
        self.neg_set_avg_center_dist = neg_set_avg_center_dist
        self.pos_set_max_center_dist = pos_set_max_center_dist
        self.pos_set_max_center_dist_c = pos_set_max_center_dist_c
        self.neg_set_max_center_dist = neg_set_max_center_dist
        self.neg_set_max_center_dist_c = neg_set_max_center_dist_c
        self.logo_png_file = logo_png_file
        self.wrs_pval_two_sided = wrs_pval_two_sided
        self.wrs_test_stat_two_sided = wrs_test_stat_two_sided
        self.wrs_pval_greater = wrs_pval_greater
        self.wrs_test_stat_greater = wrs_test_stat_greater
        self.wrs_pval_less = wrs_pval_less
        self.wrs_test_stat_less = wrs_test_stat_less
        self.dist_plot_counts_dic = dist_plot_counts_dic


################################################################################

class MotifInfos:
    """
    Stores database motif inofs.

    Actual motifs (MEME motif format sequence motif, covariance  structure model) 
    are currently provided in separate file(s).
    
    Mandatory infos:
    RBP_motif_ID
    RBP_name
    Motif_type	
    Optional infos:
    Organism
    Gene_ID
    Function_IDs
    Reference
    Experiment
    Comments

    """
    def __init__(self,
                 rbp_id: str,
                 motif_id: str,
                 motif_type = str,
                 organism = "-",
                 gene_id = "-",
                 reference = "-",
                 experiment = "-",
                 comments = "-",
                 function_ids = None) -> None:
        self.rbp_id = rbp_id
        self.motif_id = motif_id
        self.motif_type = motif_type
        self.organism = organism
        self.gene_id = gene_id
        self.reference = reference
        self.experiment = experiment
        self.comments = comments
        if function_ids is None:
            self.function_ids = []
        else:
            self.function_ids = function_ids


################################################################################

class MotifStats:
    """
    Stores motif hit stats.

    data_id, method_id, run_id, motif_db stored in RBPStats object, linked 
    by dictionary (common internal_id).
    
    hit_id : chr:s-e(+)motif_id

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
                 matched_seq: Optional[str] = None,
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
        self.matched_seq = matched_seq
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
            internal_id = cols[20]
            if internal_id == "internal_id":
                continue
            hit_id = "%s:%s-%s(%s)%s" %(cols[7], cols[8], cols[9], cols[10], cols[6])
            
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
            motif_stats.matched_seq = cols[19]
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
                 median_reg_len = 0.0,
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
            rbp_stats.median_reg_len = float(cols[7])
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
                 center_dist: Optional[int] = None,  # Distance of motif center position to region center (rbpbench nemo).
                 genome: Optional[str] = None) -> None:
        super().__init__(chr_id, start, end, strand, score, genome)
        self.motif_id = motif_id
        self.seq_name = seq_name
        self.pval = pval
        self.qval = qval
        self.matched_seq = matched_seq
        self.center_dist = center_dist
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
        return f"{self.chr_id}:{self.start}-{self.end}({self.strand}){self.motif_id}"


################################################################################

def get_regex_hits(regex, regex_id, seqs_dic,
                   step_size_one=False,
                   seq_based=False,
                   reg_dic=None,
                   use_motif_regex_id=False):
    """
    Given a regular expression (regex), get all hits in sequence dictionary.
    Store hits as FimoHit objects in list.

    use_motif_regex_id:
        If True, store regex_id as motif ID. If False, store regex as motif ID. 

    seq_based:
        If true, input to regex search was sequences, so use coordinates as they are (not genomic + relative).
        Also seq_name is the sequence ID, and does not have to have format like "chr6:35575787-35575923(-)"
    reg_dic:
        If provided, use this dictionary (mapping reg_id -> "chr1:100-200(+)", take value string) 
        to get genomic coordinates from get_genomic_coords_from_seq_name, instead of seq_name.
        This can be applied in cases where seq_name is some ID that does not contain region info, but
        motif hit coordinates are still relative to the region.

    """

    assert is_valid_regex(regex), "invalid regex given"

    hits_dic = search_regex_in_seqs_dic(regex, seqs_dic,
                                        step_size_one=step_size_one,
                                        case_sensitive=True)

    regex_hits_list = []

    motif_id = regex
    if use_motif_regex_id:
        motif_id = regex_id

    for hit in hits_dic:
        seq_name = hit
        for hit_info in hits_dic[hit]:
            start = hit_info[0]  # 0-based.
            end = hit_info[1]  # 1-based.
            matched_seq = hit_info[2]  # matched sequence.

            if seq_based:

                regex_hit = FimoHit(chr_id=seq_name, 
                                start=start+1, 
                                end=end,
                                strand="+", 
                                score=-1.0, 
                                motif_id=motif_id, 
                                seq_name=seq_name, 
                                pval=-1.0, 
                                qval=-1.0,
                                seq_s=start+1,
                                seq_e=end,
                                matched_seq=matched_seq)

            else:

                region_info_seq_name = seq_name

                if reg_dic is not None:
                    region_info_seq_name = reg_dic[seq_name]

                gen_motif_coords = get_genomic_coords_from_seq_name(region_info_seq_name, start, end,
                                                                    one_based_start=True)

                regex_hit = FimoHit(chr_id=gen_motif_coords[0], 
                                start=gen_motif_coords[1], 
                                end=gen_motif_coords[2],
                                strand=gen_motif_coords[3], 
                                score=-1.0, 
                                motif_id=motif_id, 
                                seq_name=seq_name, 
                                pval=-1.0, 
                                qval=-1.0,
                                seq_s=start+1,
                                seq_e=end,
                                matched_seq=matched_seq)

            regex_hits_list.append(regex_hit)

    return regex_hits_list


################################################################################

def bed_get_core_rel_reg_dic(core_reg_dic, in_bed):
    """
    Based on core region dictionary, containing chromosome coordinates, 
    get relative region dictionary (i.e. relative to sequence not to chromosome start).

    test3.bed:
    chr1	100	200	r1	0	+
    chr1	100	200	r2	0	-

    >>> core_reg_dic = {"r1": ["chr1", 120, 130, "+"], "r2": ["chr1", 120, 130, "-"], "r3": ["chr1", 200, 210, "+"], "r4": ["chr1", 200, 210, "-"]}
    >>> in_bed = "test_data/test3.bed"
    >>> core_rel_reg_dic = bed_get_core_rel_reg_dic(core_reg_dic, in_bed)
    >>> print(core_rel_reg_dic["r1"])
    ['r1', 20, 30, '+']
    >>> print(core_rel_reg_dic["r2"])
    ['r2', 70, 80, '+']
    >>> print(core_rel_reg_dic["r3"])
    ['r3', 0, 10, '+']
    >>> print(core_rel_reg_dic["r4"])
    ['r4', 90, 100, '+']
    
    """
    core_rel_reg_dic = {}

    with open(in_bed) as f:
        for line in f:
            cols = line.strip().split("\t")
            reg_chr = cols[0]
            reg_s = int(cols[1])  # 0-based.
            reg_e = int(cols[2])  # 1-based.
            reg_id = cols[3]
            reg_pol = cols[5]
            core_chr = cols[0]
            core_s = core_reg_dic[reg_id][1] # 0-based.
            core_e = core_reg_dic[reg_id][2] # 1-based.
            core_pol = core_reg_dic[reg_id][3]
            assert reg_chr == core_chr, "chromosome mismatch between core and region for ID %s in BED file %s" %(reg_id, in_bed)
            assert reg_pol == core_pol, "polarity mismatch between core and region for ID %s in BED file %s" %(reg_id, in_bed)

            core_len = core_e - core_s

            # + strand.
            rel_s = core_s - reg_s
            rel_e = rel_s + core_len
            if core_pol == "-":
                rel_s = reg_e - core_e
                rel_e = rel_s + core_len

            core_rel_reg_dic[reg_id] = [reg_id, rel_s, rel_e, "+"]

    f.closed

    return core_rel_reg_dic


################################################################################

def bed_get_region_str_len_dic(in_bed):
    """
    Read in BED file and store region ID -> region string in dictionary.
    Also return region ID -> region length mappings.
    Format of region string:
    chr1:100-200(+)
    
    """
    reg_dic = {}
    len_dic = {}

    with open(in_bed) as f:
        for line in f:
            cols = line.strip().split("\t")
            reg_id = cols[3]
            reg_len = int(cols[2]) - int(cols[1])
            reg_str = "%s:%s-%s(%s)" %(cols[0], cols[1], cols[2], cols[5])
            reg_dic[reg_id] = reg_str
            len_dic[reg_id] = reg_len
    f.closed

    return reg_dic, len_dic


################################################################################

def intervals_overlap(start1, end1, start2, end2):
    """
    Check if two intervals [start1, end1] and [start2, end2] overlap.
    All indices are 1-based.

    >>> intervals_overlap(1, 10, 5, 15)
    True
    >>> intervals_overlap(1, 10, 11, 15)
    False
    >>> intervals_overlap(5, 20, 10, 15)
    True
    >>> intervals_overlap(25, 30, 20, 25)
    True
    >>> intervals_overlap(40, 50, 20, 25)
    False

    """
    return start1 <= end2 and end1 >= start2


################################################################################

def filter_out_center_motif_hits(hits_list, core_rel_reg_dic,
                                 allow_overlaps=False):
    """
    Filter positive regions hits list (FIMO, CMSEARCH), based on given core region
    (via core_rel_reg_dic) with which the hit should not overlap.
    core_rel_reg_dic: core region dictionary with:
    region_ID -> [chr, seq_s, seq_e, pol] (seq_s 0-based, seq_e 1-based).
    seq_s : start on sequence snippet (0-based).
    seq_e : end on sequence snippet (1-based).
    If hit overlaps with core region, filter it out, otherwise keep it and 
    store in flt_hits_list.

    >>> fh1 = FimoHit("chr1", 140, 150, "+", 0.0, "motif1", "pos1", 0.0, seq_s=40, seq_e=50)
    >>> fh2 = FimoHit("chr1", 120, 130, "+", 0.0, "motif1", "pos1", 0.0, seq_s=20, seq_e=30)
    >>> hits_list = [fh1, fh2]
    >>> core_rel_reg_dic = {"pos1": ["pos1", 45, 55, "+"]}
    >>> flt_hits_list = filter_out_center_motif_hits(hits_list, core_rel_reg_dic)
    >>> print(flt_hits_list[0])
    chr1:120-130(+)motif1
    >>> print(flt_hits_list[0].center_dist)
    -26

    """

    flt_hits_list = []

    for hit in hits_list:

        motif_id = hit.motif_id
        seq_name = hit.seq_name
        hit_seq_s = hit.seq_s  # hit start on sequence snippet, already 1-based.
        hit_seq_e = hit.seq_e
        core_seq_s = core_rel_reg_dic[seq_name][1] + 1  # make 1-based.
        core_seq_e = core_rel_reg_dic[seq_name][2]

        assert hit_seq_s <= hit_seq_e, "hit_seq_s > hit_seq_e"
        assert core_seq_s <= core_seq_e, "core_seq_s > core_seq_e"

        if not allow_overlaps:
            if intervals_overlap(hit_seq_s, hit_seq_e, core_seq_s, core_seq_e):
                continue

        hit_center_pos = get_center_position(hit_seq_s-1, hit_seq_e)
        core_center_pos = get_center_position(core_seq_s-1, core_seq_e)

        # Can be negative if hit upstream of center, or positive if hit downstream.
        hit.center_dist = hit_center_pos - core_center_pos

        flt_hits_list.append(hit)

    return flt_hits_list


################################################################################

def filter_out_neg_center_motif_hits(neg_hits_list, core_rel_reg_dic,
                                     allow_overlaps=False):
    """
    Filter negative regions hits list (FIMO, CMSEARCH), based on given positive 
    core region (via core_rel_reg_dic) with which the hit should not overlap.
    core_rel_reg_dic: core region dictionary with:
    region_ID -> [chr, seq_s, seq_e, pol] (seq_s 0-based, seq_e 1-based).
    seq_s : start on sequence snippet (0-based).
    seq_e : end on sequence snippet (1-based).
    If hit overlaps with core region, filter it out, otherwise keep it and 
    store in flt_hits_list.
    negative seq_name format: pos1;neg1 pos1;neg2 etc.
    So first part (pos1) identifies which core region to use (here pos1) for
    filtering.

    Currently just assume both positive and negative regions have same lengths, 
    thus mask same relative core region.

    >>> fh1 = FimoHit("chr1", 140, 150, "+", 0.0, "motif1", "pos1;neg1", 0.0, seq_s=40, seq_e=50)
    >>> fh2 = FimoHit("chr1", 120, 130, "+", 0.0, "motif1", "pos1;neg1", 0.0, seq_s=20, seq_e=30)
    >>> hits_list = [fh1, fh2]
    >>> core_rel_reg_dic = {"pos1": ["pos1", 45, 55, "+"]}
    >>> flt_hits_list = filter_out_neg_center_motif_hits(hits_list, core_rel_reg_dic)
    >>> print(flt_hits_list[0])
    chr1:120-130(+)motif1
    
    """

    flt_neg_hits_list = []

    for hit in neg_hits_list:

        # motif_id = hit.motif_id
        seq_name = hit.seq_name
        hit_seq_s = hit.seq_s  # hit start on sequence snippet, already 1-based.
        hit_seq_e = hit.seq_e
        pos_seq_name = seq_name.split(";")[0]
        core_seq_s = core_rel_reg_dic[pos_seq_name][1] + 1  # make 1-based.
        core_seq_e = core_rel_reg_dic[pos_seq_name][2]

        assert hit_seq_s <= hit_seq_e, "hit_seq_s > hit_seq_e"
        assert core_seq_s <= core_seq_e, "core_seq_s > core_seq_e"

        if not allow_overlaps:
            if intervals_overlap(hit_seq_s, hit_seq_e, core_seq_s, core_seq_e):
                continue

        hit_center_pos = get_center_position(hit_seq_s-1, hit_seq_e)
        core_center_pos = get_center_position(core_seq_s-1, core_seq_e)

        # Can be negative if hit upstream of center, or positive if hit downstream.
        hit.center_dist =  hit_center_pos - core_center_pos

        flt_neg_hits_list.append(hit)

    return flt_neg_hits_list


################################################################################

def read_in_fimo_results(fimo_tsv,
                         seq_based=False,
                         reg_dic=None,
                         only_best_hits=False,
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

    With transcripts/sequences (rbpbench searchseq):
    motif_id	motif_alt_id	sequence_name	start	stop	strand	score	p-value	q-value	matched_sequence
    PUM1_1		ENST00000561978.1	1242	1248	+	9.84848	0.000296		
    PUM1_1		ENST00000561978.1	1305	1311	+	6.9798	0.000768		
    PUM1_2		ENST00000561978.1	1127	1135	+	-2.35354	0.000906		
    PUM1_2		ENST00000561978.1	1949	1957	+	6.29293	0.00025		
    PUM1_3		ENST00000561978.1	1414	1421	+	3.61616	0.000958		

    UPDATES:
    Do not read in q-value or matched_sequence, as these are not produced when 
    using run_fast_fimo().

    seq_based:
        If true, input to FIMO search was sequences, so use coordinates as they are (not genomic + relative).
        Also seq_name is the sequence ID, and does not have to have format like "chr6:35575787-35575923(-)"
        Note that if prediction is on subsequences of transcripts, seq_based should be False too.
    reg_dic:
        If provided, use this dictionary (mapping reg_id -> "chr1:100-200(+)", take value string) 
        to get genomic coordinates from get_genomic_coords_from_seq_name, instead of seq_name.
        This can be applied in cases where seq_name is some ID that does not contain region info, but
        motif hit coordinates are still relative to the region.
    only_best_hits:
        Keep only best hits for each motif ID and sequence/site combination. I.e., hit with lowest p-value.
        
    """

    fimo_hits_list = []

    with open(fimo_tsv) as f:
        for line in f:
            # if re.search("^#", line):
            if line.startswith("#"):
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

            if seq_based:

                fimo_hit = FimoHit(chr_id=seq_name, 
                                start=motif_s+1, 
                                end=motif_e,
                                strand="+", 
                                score=score, 
                                motif_id=motif_id, 
                                seq_name=seq_name, 
                                pval=pval, 
                                qval=qval,
                                seq_s=motif_s+1,
                                seq_e=motif_e,
                                matched_seq=matched_seq)
                
            else:

                region_info_seq_name = seq_name

                if reg_dic is not None:
                    region_info_seq_name = reg_dic[seq_name]

                gen_motif_coords = get_genomic_coords_from_seq_name(region_info_seq_name, motif_s, motif_e,
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

    if only_best_hits and fimo_hits_list:
        print("--greatest-hits enabled. Keep only best FIMO motif hits ... ")
        print("# of FIMO motif hits before best hit filtering: %i" %(len(fimo_hits_list)))
        fimo_hits_list = get_best_fimo_hits(fimo_hits_list)
        print("# of FIMO motif hits after best hit filtering:  %i" %(len(fimo_hits_list)))

    return fimo_hits_list


################################################################################

def get_best_fimo_hits(fimo_hits_list):
    """
    Filter FIMO hits list, keep only best hit for each motif ID and sequence combination.
    Best hit is hit with lowest p-value.
    
    """
    assert fimo_hits_list, "fimo_hits_list is empty"

    filt_fimo_hits_list = []
    best_list_idx_dic = {}  # "motif_id,seq_name" -> best list index.

    for idx, fimo_hit in enumerate(fimo_hits_list):
        seq_name = fimo_hit.seq_name
        motif_id = fimo_hit.motif_id
        p_val = fimo_hit.pval
        comb_id = "%s,%s" %(motif_id, seq_name)
        if comb_id not in best_list_idx_dic:
            best_list_idx_dic[comb_id] = idx
        else:
            best_idx = best_list_idx_dic[comb_id]
            if fimo_hits_list[best_idx].pval > p_val:
                best_list_idx_dic[comb_id] = idx

    for comb_id in best_list_idx_dic:
        best_idx = best_list_idx_dic[comb_id]
        filt_fimo_hits_list.append(fimo_hits_list[best_idx])

    return filt_fimo_hits_list


################################################################################

def get_best_cmsearch_hits(cmsearch_hits_list):
    """
    Filter CMSEARCH hits list, keep only best hit for each motif ID and sequence combination.
    Best hit is hit with highest bit score.
    
    """
    assert cmsearch_hits_list, "cmsearch_hits_list is empty"

    filt_cmsearch_hits_list = []
    best_list_idx_dic = {}  # "motif_id,seq_name" -> best list index.

    for idx, cmsh_hit in enumerate(cmsearch_hits_list):
        seq_name = cmsh_hit.seq_name
        motif_id = cmsh_hit.motif_id
        score = cmsh_hit.score  # compare cmsearch bit score (the higher the better).
        comb_id = "%s,%s" %(motif_id, seq_name)
        if comb_id not in best_list_idx_dic:
            best_list_idx_dic[comb_id] = idx
        else:
            best_idx = best_list_idx_dic[comb_id]
            if cmsearch_hits_list[best_idx].score < score:
                best_list_idx_dic[comb_id] = idx

    for comb_id in best_list_idx_dic:
        best_idx = best_list_idx_dic[comb_id]
        filt_cmsearch_hits_list.append(cmsearch_hits_list[best_idx])

    return filt_cmsearch_hits_list


################################################################################

def get_target_genes_with_rbp_hits(reg2annot_dic, tr2gid_dic, region_rbp_binds_dic,
                                   gid2tid_dic=None,
                                   goa_cooc_mode=3):
    """
    Get target genes dictionary with RBP hits.

    reg2annot_dic format:
    'chr20:62139082-62139128(-)': ['CDS', 'ENST00000367770.3'] 
    'chr20:62139082-62139128(+)': ['intergenic', False] 

    region_rbp_binds_dic format:
    'chr20:62139082-62139128(-)': [False, False, False]
    
    >>> reg2annot_dic = {'chr1:1000-2000(+)': ['CDS', 'ENST6666'], 'chr1:1000-2000(-)': ['intergenic', False], 'chr1:5000-6000(-)': ['intron', 'ENST6667'], 'chr1:5000-6000(+)': ['CDS', 'ENST6668']}
    >>> tr2gid_dic = {'ENST6666': 'GID1', 'ENST6667': 'GID2', 'ENST6668': 'GID3'}
    >>> region_rbp_binds_dic = {'chr1:1000-2000(+)': [True, True, True], 'chr1:1000-2000(-)': [True, True, True], 'chr1:5000-6000(-)': [False, True, True], 'chr1:5000-6000(+)': [False, False, False]}
    >>> get_target_genes_with_rbp_hits(reg2annot_dic, tr2gid_dic, region_rbp_binds_dic, goa_cooc_mode=3)
    {'GID1': 1}
    >>> get_target_genes_with_rbp_hits(reg2annot_dic, tr2gid_dic, region_rbp_binds_dic, goa_cooc_mode=2)
    {'GID1': 1, 'GID2': 1}
    >>> get_target_genes_with_rbp_hits(reg2annot_dic, tr2gid_dic, region_rbp_binds_dic, goa_cooc_mode=1)
    {'GID1': 1, 'GID2': 1, 'GID3': 1}
    
    """

    target_genes_dic = {}
    for reg_id in region_rbp_binds_dic:
        tr_id = reg2annot_dic[reg_id][1]
        if tr_id:
            gid = tr2gid_dic[tr_id]
            if gid2tid_dic is not None:
                gid2tid_dic[gid] = tr_id
            add_gene = False
            if goa_cooc_mode == 1:
                add_gene = True
            elif goa_cooc_mode == 2:
                # Check if region_rbp_binds_dic[reg_id] list has at least one True value.
                if any(region_rbp_binds_dic[reg_id]):
                    add_gene = True
            elif goa_cooc_mode == 3:
                # Check if region_rbp_binds_dic[reg_id] list has only True values.
                if all(region_rbp_binds_dic[reg_id]):
                    add_gene = True

            if add_gene:
                if gid not in target_genes_dic:
                    target_genes_dic[gid] = 1
                else:
                    target_genes_dic[gid] += 1

    return target_genes_dic


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

    # pattern_check = re.compile(r"\w+:\d+-\d+\([+|-]\).+")
    pattern_extract = re.compile(r"^(.+):(\d+)-(\d+)\(([+|-])\)(.+)")

    for fh_str in unique_motifs_dic[rbp_id]:
        # if re.search("\w+:\d+-\d+\([+|-]\).+", fh_str):
        #   m = re.search("(\w+):(\d+)-(\d+)\(([+|-])\)(.+)", fh_str)

        m = pattern_extract.search(fh_str)
        assert m, "invalid fh_str format given (%s)" % (fh_str)

        chr_id = m.group(1)
        reg_s = int(m.group(2))
        reg_e = int(m.group(3))
        strand = m.group(4)
        motif_id = m.group(5)

        if one_based_start:
            reg_s -= 1

        OUTMRBED.write("%s\t%i\t%i\t%s\t0\t%s\n" % (chr_id, reg_s, reg_e, motif_id, strand))

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
    chr_id:start-end(strand)motif_id

    """
    
    OUTMRBED = open(out_bed, "w")

    pattern_extract = re.compile(r"^(\w+?):(\d+)-(\d+)\(([+|-])\)(.+)")

    for fh_str in unique_motifs_dic:
        # if re.search("^\w+?:\d+-\d+\([+|-]\).+", fh_str):
        #     m = re.search("^(\w+?):(\d+)-(\d+)\(([+|-])\)(.+)", fh_str)
        m = pattern_extract.search(fh_str)
        assert m, "invalid fh_str format given (%s)" % (fh_str)

        chr_id = m.group(1)
        reg_s = int(m.group(2))
        reg_e = int(m.group(3))
        strand = m.group(4)
        motif_id = m.group(5)

        if one_based_start:
            reg_s -= 1

        OUTMRBED.write("%s\t%i\t%i\t%s\t0\t%s\n" % (chr_id, reg_s, reg_e, motif_id, strand))

    OUTMRBED.close()


################################################################################

def read_in_cmsearch_results(in_tab,
                             check=True,
                             seq_based=False,
                             reg_dic=None,
                             only_best_hits=False,
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

    
    reg_dic:
        If provided, use this dictionary (mapping reg_id -> "chr1:100-200(+)", take value string) 
        to get genomic coordinates from get_genomic_coords_from_seq_name, instead of seq_name.
        This can be applied in cases where seq_name is some ID that does not contain region info, but
        motif hit coordinates are still relative to the region.

    """

    if check:
        assert os.path.exists(in_tab), "cmsearch output file %s does not exist! Check cmsearch output for errors (likely a failed cmsearch call due to invalid parameter settings)" % (in_tab)

    if hits_list is None:
        hits_list = []

    c_hits = 0

    with open(in_tab) as f:
        for line in f:
            if line.startswith("#"):
            # if re.search("^#", line):
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

            if seq_based:

                cmsearch_hit = CmsearchHit(chr_id=seq_name, 
                                    start=seq_s+1, 
                                    end=seq_e,
                                    strand="+",
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

            else:

                region_info_seq_name = seq_name

                if reg_dic is not None:
                    region_info_seq_name = reg_dic[seq_name]

                gen_motif_coords = get_genomic_coords_from_seq_name(region_info_seq_name, seq_s, seq_e,
                                                                    one_based_start=True)

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

    if only_best_hits and hits_list:
        print("--greatest-hits enabled. Keep only best CMSEARCH hits ... ")
        print("# of CMSEARCH hits before best hit filtering: %i" %(len(hits_list)))
        hits_list = get_best_cmsearch_hits(hits_list)
        print("# of CMSEARCH hits after best hit filtering:  %i" %(len(hits_list)))
        c_hits = len(hits_list)

    return hits_list, c_hits


################################################################################

def get_matched_seq(seq_id, seqs_dic, match_s, match_e):
    """
    Get motif hit sequence (matched sequence, subsequence) given:
    sequence dictionary seqs_dic, sequence name seq_id, 
    and start+end coordinates (match_s, match_e, both 1-based).

    >>> seqs_dic = {"seq1": "AAAAACGTAAAA"}
    >>> get_matched_seq("seq1", seqs_dic, 6, 8)
    'CGT'

    """
    assert seq_id in seqs_dic, "sequence ID \"%s\" not in seqs_dic" %(seq_id)
    # assert match_s <= match_e, "match_s > match_e for seq_id %s" %(seq_id)
    seq = seqs_dic[seq_id]
    # len_seq = len(seq)
    # assert match_s <= len_seq, "match_s > len_seq (%i > %i) for seq_id %s" %(seq_id)
    # assert match_e <= len_seq, "match_e > len_seq (%ifor seq_id %s" %(seq_id)
    matched_seq = seq[match_s-1:match_e]
    return matched_seq


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
                 center_dist: Optional[int] = None,  # Distance of motif center position to region center (rbpbench nemo).
                 matched_seq: Optional[str] = None,
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
        self.query_name = query_name
        self.center_dist = center_dist
        self.matched_seq = matched_seq

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
        return f"{self.chr_id}:{self.start}-{self.end}({self.strand}){self.motif_id}"


################################################################################

def insert_line_breaks(sequence,
                       line_len=60):
    """
    Insert line breaks for inserting sequence into hover box.
    
    >>> insert_line_breaks("AAAACCCCGG", line_len=4)
    'AAAA<br>CCCC<br>GG'

    """
    # return '<br>'.join(sequence[i:i+line_len] for i in range(0, len(sequence), 40))
    return '<br>'.join(sequence[i:i+line_len] for i in range(0, len(sequence), line_len))


################################################################################

def join_motif_hits(motif_hits_list,
                    motifs_per_line=3,
                    line_break_char="\n"):
    """
    Join motif hits list into string by groups of motifs_per_line.
    Define line break character via line_break_char.

    """

    return line_break_char.join(" ".join(motif_hits_list[i:i+motifs_per_line]) for i in range(0, len(motif_hits_list), motifs_per_line))


################################################################################

def remove_special_chars_from_str(check_str,
                                  reg_ex=r'[^A-Za-z0-9_-]+'):
    """
    Remove special characters from string.

    reg_ex:
        Regular expression defining what to keep.

    Special regex characters:
    . (dot)
    ^ (caret)
    $ (dollar sign)
    * (asterisk)
    + (plus sign)
    ? (question mark)
    {} (curly braces)
    [] (square brackets)
    () (parentheses)
    | (pipe)
    backslash

    >>> check_str = r"{_}[-](_)\\V/"
    >>> remove_special_chars_from_str(check_str)
    '_-_V'
    >>> check_str = ""
    >>> remove_special_chars_from_str(check_str)
    ''
    >>> check_str = "AC.+?GA[AC]A\\\\C(AAA)C;C.{2,8}AC"
    >>> remove_special_chars_from_str(check_str, reg_ex=r"[ ;\\()]")
    'AC.+?GA[AC]ACAAACC.{2,8}AC'

    """
    # To remove special regex chars: r"[.^$*+?{}[\]()|\]"

    check_str = check_str.replace(r"\t", "").replace(r"\n", "").replace("\\", "")
    clean_string = re.sub(reg_ex, '', check_str)
    return clean_string


################################################################################

def get_motif_id_from_str_repr(hit_str_repr):
    """
    From motif string representation:
    chr1:100-200(-)motif_id
    return motif_id

    ENST00000561978.1:1242-1248(+)PUM1_1

    >>> hit_str_repr = "chr6:66-666(-)satan6666"
    >>> get_motif_id_from_str_repr(hit_str_repr)
    'satan6666'
    >>> hit_str_repr = "chr6:66-666(-)CGGAC.{10,25}[CA]CA[CT]"
    >>> get_motif_id_from_str_repr(hit_str_repr)
    'CGGAC.{10,25}[CA]CA[CT]'

    """

    if re.search(r"^.+:\d+-\d+\([+|-]\).+", hit_str_repr):
        m = re.search(r"^.+?\)(.+)", hit_str_repr)
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

def read_file_content_into_str_var(file):
    """
    Read in file content and store in string variable.
    """
    str_var = ""
    with open(file, 'r') as f:
        str_var = f.read()
    f.closed
    return str_var


################################################################################

def create_color_scale(min_value, max_value, colors):
    """
    Create a color scale with equal intervals between colors.

    Parameters:
    min_value (float): The minimum value of the data.
    max_value (float): The maximum value of the data.
    colors (list of str): The colors to use in the color scale.

    Returns:
    list of [float, str]: The color scale.
    """

    # Calculate the range and interval
    range_ = max_value - min_value
    interval = range_ / (len(colors) - 1)

    # Create the color scale
    color_scale = [[min_value + i * interval, color] for i, color in enumerate(colors)]

    return color_scale


################################################################################

def rgb_to_hex(rgb):
    """
    Convert RGB to hex.

    """
    rgb = rgb.replace('rgb(', '').replace(')', '').split(',')
    return '#{:02x}{:02x}{:02x}'.format(int(rgb[0]), int(rgb[1]), int(rgb[2]))


################################################################################

def create_cooc_plot_plotly(df, pval_cont_lll, plot_out,
                            max_motif_dist=50,
                            min_motif_dist=0,
                            id1="RBP1",
                            id2="RBP2",
                            ids="RBPs",
                            include_plotlyjs="cdn",
                            full_html=False):
    """
    Plot co-occurrences as heat map with plotly.

    https://plotly.github.io/plotly.py-docs/generated/plotly.io.write_html.html

    """

    # Threshold value for grey coloring.
    # threshold = 0.2
    # # color_scale = ["grey", "darkblue", "blue", "purple", "red", "yellow"]
    # # color_scale = ["grey", "blue", "orange", "yellow"]
    # color_scale = ["grey", "darkblue", "purple", "yellow"]

    # plot = px.imshow(df, color_continuous_scale=color_scale, color_continuous_midpoint=threshold)
    # plot = px.imshow(df)

    # Define a custom color scale
    # color_scale = [[0, 'grey'], [0.01, 'darkblue'], [1, 'yellow']]

    # color_scale = [[0, "midnightblue"], [0, 'darkblue'], [0.25, 'blue'], [0.5, 'purple'], [0.75, 'red'], [1, 'yellow']]
    
    # Light blue to dark blue gradient.
    # colors = ['#E0FFFF', '#B0E0E6', '#87CEEB', '#4682B4', '#0000FF', '#00008B', '#000080']  
    # colors = ['darkblue', 'blue', 'purple', 'red', 'orange', 'yellow']
    # colors = ['green', 'red']

    min_val = 0.01
    max_val = 1  # df.max().max() # color scale values need to be 0 .. 1.

    # color_scale = create_color_scale(min_val, max_val, colors)
    # color_scale.insert(0, [0, "white"])  # lightgray ? dimgray (old color)


    # Define the Blues color scale.
    blues_colors = px.colors.sequential.Blues
    blues_colors_hex = [rgb_to_hex(color) for color in blues_colors]
    # blues_colors_hex: ['#f7fbff', '#deebf7', '#c6dbef', '#9ecae1', '#6baed6', '#4292c6', '#2171b5', '#08519c', '#08306b']

    color_scale = create_color_scale(min_val, max_val, blues_colors_hex)
    color_scale.insert(0, [0, "white"])  # lightgray ? dimgray was old color.

    # # Insert white at the beginning of the color scale
    # custom_color_scale = [[0, 'white']] + [[i / (len(blues_colors) - 1), color] for i, color in enumerate(blues_colors)]

    # Set the color scale range according to the data (minimum 0 .. 1).
    zmin = min(0, df.min().min())
    zmax = max(1, df.max().max())

    # color_scale = [[0, 'midnightblue'], [0.01, 'red'], [1, 'yellow']]
    # print("color_scale:")
    # print(color_scale)

    plot = px.imshow(df, color_continuous_scale=color_scale, zmin=zmin, zmax=zmax)
    # plot = px.imshow(df)

    hover_content = (
        '1) ' + id1 + ': %{x}<br>%{customdata[8]}'
        '2) ' + id2 + ': %{y}<br>%{customdata[9]}'
        '3) p-value: %{customdata[0]}<br>'
        '4) p-value after filtering: %{customdata[1]}<br>%{customdata[7]}'
        '5) ' + ids + ': %{customdata[2]}<br>'
        '6) Counts: %{customdata[3]}<br>'
        '7) Mean minimum motif distance (nt): %{customdata[4]}<br>'
        '8) Motif pairs within ' + str(max_motif_dist) + ' nt (%): %{customdata[5]}<br>'
        '9) Correlation: %{customdata[6]}<br>'
        '10) -log10(p-value after filtering): %{z}<extra></extra>'
    )

    plot.update(data=[{'customdata': pval_cont_lll,
                      'hovertemplate': hover_content}])
    plot.update_layout(plot_bgcolor='white')
    plot.write_html(plot_out,
                    full_html=full_html,
                    include_plotlyjs=include_plotlyjs)


################################################################################

def create_annot_comp_plot_plotly(dataset_ids_list, annots_ll, 
                                  annot_ids_list, annot2color_dic,
                                  plot_out,
                                  include_plotlyjs="cdn",
                                  full_html=False):

    """
    Create plotly 2d scatter plot of PCA reduced genomic annotations.

    annots_ll:
        List of list containing annotation frequencies (0 to 1). In order of dataset_ids_list.
    annot_ids_list:
        List of annotation IDs (e.g. "CDS", "UTR", "Intron", "Intergenic").
        In order of annots_ll, so that annots_ll[i][0] corresponds to annot_ids_list[0].

    """
    # Get highest percentage for every dataset.
    highest_perc_annot_list = []
    perc_annot_str_list = []
    for annots_l in annots_ll:
        perc_annot_string = "<br>"
        highest_perc_annot = "-"
        highest_freq = 0.0
        for idx, annot_freq in enumerate(annots_l):
            annot = annot_ids_list[idx]
            if annot_freq > highest_freq:
                highest_freq = annot_freq
                highest_perc_annot = annot
            # Make percentage out of frequency and round to two decimal places.
            annot_perc = round(annot_freq * 100, 2)
            perc_annot_string += annot + ": " + str(annot_perc) + '%<br>'
        highest_perc_annot_list.append(highest_perc_annot)
        perc_annot_str_list.append(perc_annot_string)
    

    pca = PCA(n_components=2)  # Reduce data to 2 dimensions.
    data_2d_pca = pca.fit_transform(annots_ll)
    df = pd.DataFrame(data_2d_pca, columns=['PC1', 'PC2'])
    df['Dataset ID'] = dataset_ids_list
    df['Highest % annotation'] = highest_perc_annot_list
    df['Annotation percentages'] = perc_annot_str_list

    explained_variance = pca.explained_variance_ratio_ * 100

    fig = px.scatter(
        df,  # Use the DataFrame directly
        x='PC1',
        y='PC2',
        color='Highest % annotation',
        color_discrete_map=annot2color_dic,
        # title='2D Visualization with Dataset IDs',
        labels={
            'PC1': f'PC1 ({explained_variance[0]:.2f}% variance)',
            'PC2': f'PC2 ({explained_variance[1]:.2f}% variance)'
        },
        hover_name='Dataset ID',
        hover_data=['Highest % annotation', 'Annotation percentages']
    )

    fig.update_traces(
        hovertemplate='<b>%{hovertext}</b><br>Annotation percentages: %{customdata[1]}<extra></extra>'
    )

    fig.update_traces(marker=dict(size=12, line=dict(width=0.75, color='white')))
    # fig.update_layout(margin=dict(l=0, r=0, b=0, t=0))
    fig.write_html(plot_out, full_html=full_html, include_plotlyjs=include_plotlyjs)


################################################################################

def create_pca_motif_sim_sig_plot_plotly(motif_ids_list, motif_sim_ll, 
                                         motif_sim_stats_dic,
                                         plot_out,
                                         include_plotlyjs="cdn",
                                         full_html=False):
    """

    motif_ids_list:
    List of motif IDs, in same order as motif_sim_ll.
    motif_sim_ll:
    List of lists containing motif similarities. In order of motif_ids_list, 
    so that motif_sim_ll[i][j] corresponds to motif_ids_list[i] and motif_ids_list[j].
    motif_sim_stats_dic:
    Dictionary with motif ID as key and list of [consensus_seq, con_table_str, pval, -log10pval] as value, 
    for hoverbox information.

    """

    assert motif_ids_list, "motif_ids_list empty"

    pca = PCA(n_components=2)
    data_2d_pca = pca.fit_transform(motif_sim_ll)

    df = pd.DataFrame(data_2d_pca, columns=['PC1', 'PC2'])
    
    df['Motif ID'] = motif_ids_list
    df['Consensus sequence'] = [motif_sim_stats_dic[motif_id][0] for motif_id in motif_ids_list]
    df['Counts'] = [motif_sim_stats_dic[motif_id][1] for motif_id in motif_ids_list]
    df['p-value'] = [motif_sim_stats_dic[motif_id][2] for motif_id in motif_ids_list]
    df['-log10(p-value)'] = [motif_sim_stats_dic[motif_id][3] for motif_id in motif_ids_list]

    # colors = ['#E0FFFF', '#B0E0E6', '#87CEEB', '#4682B4', '#0000FF', '#00008B', '#000080']
    # min_val = 0.01
    # max_val = 1
    # color_scale = create_color_scale(min_val, max_val, colors)
    # color_scale.insert(0, [0, "white"])
    # lp_min = df['-log10(p-value)'].min()
    # lp_max = df['-log10(p-value)'].max()

    # Define the custom blue color scale
    # colors = ['#E0FFFF', '#B0E0E6', '#87CEEB', '#4682B4', '#0000FF', '#00008B', '#000080']
    # color_scale = [[i / (len(colors) - 1), color] for i, color in enumerate(colors)]
    # Viridis, Magma
    # intron green: #2ca02c
    # 3utr blue: #1f77b4
    # cds orange: #ff7f0e

    # blues_colors_hex: ['#f7fbff', '#deebf7', '#c6dbef', '#9ecae1', '#6baed6', '#4292c6', '#2171b5', '#08519c', '#08306b']

    color_scale = ['#c6dbef', '#9ecae1', '#6baed6', '#4292c6', '#2171b5', '#08519c', '#08306b']

    explained_variance = pca.explained_variance_ratio_ * 100

    fig = px.scatter(
        df,
        x='PC1',
        y='PC2',
        color='-log10(p-value)',
        labels={
            'PC1': f'PC1 ({explained_variance[0]:.2f}% variance)',
            'PC2': f'PC2 ({explained_variance[1]:.2f}% variance)'
        },
        hover_name='Motif ID',
        hover_data=['Consensus sequence', 'Counts', 'p-value', '-log10(p-value)'],
        # color_continuous_scale='Blues'
        color_continuous_scale=color_scale # ylgnbu
        # range_color=[lp_min, lp_max]
    )

    fig.update_traces(
        hovertemplate='motif: <b>%{hovertext}</b><br>Consensus: %{customdata[0]}<br>Counts: %{customdata[1]}<br>p-value: %{customdata[2]}<br>-log10(p-value): %{customdata[3]}<extra></extra>'
    )

    fig.update_traces(marker=dict(size=12, line=dict(width=0.75, color='white')))
    fig.update_layout(coloraxis_colorbar_title=''
                      # plot_bgcolor='white',
                      # paper_bgcolor='white',
                      # xaxis=dict(showgrid=False, zeroline=False, gridcolor='lightgray', zerolinecolor='lightgray'),
                      # yaxis=dict(showgrid=False, zeroline=False, gridcolor='lightgray', zerolinecolor='lightgray')
                      )
    fig.write_html(plot_out, full_html=full_html, include_plotlyjs=include_plotlyjs)
    # fig.update_scenes(aspectmode='cube')
    # fig.update_traces(marker=dict(size=3))
    # fig.update_layout(margin=dict(l=0, r=0, b=0, t=0))
    # fig.write_html(plot_out, full_html=full_html, include_plotlyjs=include_plotlyjs)


################################################################################

def create_pca_motif_sim_dir_plot_plotly(motif_ids_list, motif_sim_ll, 
                                         motif_sim_stats_dic,
                                         plot_out,
                                         include_plotlyjs="cdn",
                                         full_html=False):
    """

    motif_ids_list:
    List of motif IDs, in same order as motif_sim_ll.
    motif_sim_ll:
    List of lists containing motif similarities. In order of motif_ids_list, 
    so that motif_sim_ll[i][j] corresponds to motif_ids_list[i] and motif_ids_list[j].
    motif_sim_stats_dic:
    Dictionary with motif ID as key and list of
    [conseq, con_table_str, pval, log_pval, wrs_pval_greater, wrs_pval_less, log_wrs_pval]
    as value, for hoverbox information.

    """

    assert motif_ids_list, "motif_ids_list empty"

    pca = PCA(n_components=2)
    data_2d_pca = pca.fit_transform(motif_sim_ll)

    df = pd.DataFrame(data_2d_pca, columns=['PC1', 'PC2'])
    
    df['Motif ID'] = motif_ids_list
    df['Consensus sequence'] = [motif_sim_stats_dic[motif_id][0] for motif_id in motif_ids_list]
    df['Counts'] = [motif_sim_stats_dic[motif_id][1] for motif_id in motif_ids_list]
    df['p-value'] = [motif_sim_stats_dic[motif_id][2] for motif_id in motif_ids_list]
    df['WRS p-value (upstream)'] = [motif_sim_stats_dic[motif_id][4] for motif_id in motif_ids_list]
    df['WRS p-value (downstream)'] = [motif_sim_stats_dic[motif_id][5] for motif_id in motif_ids_list]
    df['-log10(WRS p-value)'] = [motif_sim_stats_dic[motif_id][6] for motif_id in motif_ids_list]

    # Get maximum -log10(WRS p-value) for colorbar scaling.
    max_abs_log10_pval = df['-log10(WRS p-value)'].abs().max()

    # Define a custom diverging color scale
    custom_color_scale = [
        [0, '#ff7f0e'],  # counts higher in upstream context
        [0.5, 'white'],  # Color around zero
        [1, '#1f77b4']  # counts higher in downstream context
    ]

    # intron green: #2ca02c
    # 3utr blue: #1f77b4
    # cds orange: #ff7f0e

    # dark blue from Blues: #08306b
    # blues_colors_hex: ['#f7fbff', '#deebf7', '#c6dbef', '#9ecae1', '#6baed6', '#4292c6', '#2171b5', '#08519c', '#08306b']


    # colors = ['#E0FFFF', '#B0E0E6', '#87CEEB', '#4682B4', '#0000FF', '#00008B', '#000080']
    # min_val = 0.01
    # max_val = 1
    # color_scale = create_color_scale(min_val, max_val, colors)
    # color_scale.insert(0, [0, "white"])
    # lp_min = df['-log10(p-value)'].min()
    # lp_max = df['-log10(p-value)'].max()

    explained_variance = pca.explained_variance_ratio_ * 100

    fig = px.scatter(
        df,
        x='PC1',
        y='PC2',
        color='-log10(WRS p-value)',
        labels={
            'PC1': f'PC1 ({explained_variance[0]:.2f}% variance)',
            'PC2': f'PC2 ({explained_variance[1]:.2f}% variance)'
        },
        hover_name='Motif ID',
        hover_data=['Consensus sequence', 'Counts', 'p-value', 'WRS p-value (upstream)', 'WRS p-value (downstream)', '-log10(WRS p-value)'],
        color_continuous_scale=custom_color_scale,
        range_color=[-max_abs_log10_pval, max_abs_log10_pval]
    )

    fig.update_traces(
        hovertemplate='motif: <b>%{hovertext}</b><br>Consensus: %{customdata[0]}<br>Counts: %{customdata[1]}<br>p-value: %{customdata[2]}<br>WRS p-value (upstream): %{customdata[3]}<br>WRS p-value (downstream): %{customdata[4]}<br>-log10(WRS p-value): %{customdata[5]}<extra></extra>'
    )

    # fig.update_traces(marker=dict(size=11))
    # fig.update_layout(coloraxis_colorbar_title='')

    fig.update_traces(marker=dict(size=12, line=dict(width=0.75, color='white')))
    fig.update_layout(coloraxis_colorbar_title=''
                      # plot_bgcolor='white',
                      # paper_bgcolor='white',
                      # xaxis=dict(showgrid=False, zeroline=False, gridcolor='lightgray', zerolinecolor='lightgray'),
                      # yaxis=dict(showgrid=False, zeroline=False, gridcolor='lightgray', zerolinecolor='lightgray')
                      )

    fig.write_html(plot_out, full_html=full_html, include_plotlyjs=include_plotlyjs)


################################################################################

def create_seq_var_plot_plotly(seq_var_ll, seq_var_kmer_l, plot_out,
                               kmer_size=1,
                               top_bottom_n=5,
                               remove_zero_val=True,
                               include_plotlyjs="cdn",
                               full_html=False):
    """
    Create sequence variation plot, use PCA to reduce dimensions of single k-mer 
    coefficients of variations in the datasets.
    
    Formats:
    seq_var_ll = [dataset_id, average_cv, single_cv1, single_cv2 ... ]
    seq_var_kmer_l = [kmer1, kmer2, ...]

    remove_zero_val:
        In top CV rankings remove zero CV values.

    """

    assert seq_var_ll, "seq_var_ll empty"
    assert seq_var_kmer_l, "seq_var_kmer_l empty"
    assert top_bottom_n <= 10, "top_bottom_n > 10"
    assert kmer_size < 6, "kmer_size > 5"

    color_scale = ['#c6dbef', '#9ecae1', '#6baed6', '#4292c6', '#2171b5', '#08519c', '#08306b']

    single_var_ll = [seq_var_l[2:] for seq_var_l in seq_var_ll]

    kmer_cvs_str_list = []

    # Construct k-mer CV strings for hover box.
    for idx1, single_var_l in enumerate(single_var_ll):

        # Average CV.
        avg_cv = seq_var_ll[idx1][1]

        kmer2cv_dic = {}
        for idx2, kmer_cv in enumerate(single_var_l):
            kmer = seq_var_kmer_l[idx2]
            kmer2cv_dic[kmer] = kmer_cv

        kmer_cv_str = '<span style="font-family: \'Courier New\', monospace;">'

        if kmer_size < 3:
            kmer_cv_str += "Single k-mer CVs:<br>"
            for kmer, kmer_cv in sorted(kmer2cv_dic.items(), key=lambda x: x[1], reverse=False):
                kmer_cv_str += kmer + ": " + str(round(kmer_cv, 5)) + "<br>"
        else:
            kmer_cv_str += "Top %i k-mer CVs:<br>" %(top_bottom_n)
            top_idx = 0
            for kmer, kmer_cv in sorted(kmer2cv_dic.items(), key=lambda x: x[1], reverse=False):
                if remove_zero_val and kmer_cv == 0:
                    continue
                kmer_cv_str += kmer + ": " + str(round(kmer_cv, 5)) + "<br>"
                if top_idx == top_bottom_n-1:
                    break
                top_idx += 1
            kmer_cv_str += "Bottom %i k-mer CVs:<br>" %(top_bottom_n)
            for kmer, kmer_cv in sorted(kmer2cv_dic.items(), key=lambda x: x[1], reverse=True)[:top_bottom_n]:
                kmer_cv_str += kmer + ": " + str(round(kmer_cv, 5)) + "<br>"

        kmer_cv_str += "Average CV: " + str(round(avg_cv, 5)) + '</span><extra></extra>'

        kmer_cvs_str_list.append(kmer_cv_str)

    exp_n_cvs = 4 ** kmer_size
    assert len(single_var_ll[0]) == exp_n_cvs, "unexpected number of CVs in seq_var_ll"

    pca = PCA(n_components=2)
    data_2d_pca = pca.fit_transform(single_var_ll)
    df = pd.DataFrame(data_2d_pca, columns=['PC1', 'PC2'])

    df['Dataset ID'] = [seq_var_l[0] for seq_var_l in seq_var_ll]
    df['k-mer CVs'] = kmer_cvs_str_list
    df['Average CV'] = [seq_var_l[1] for seq_var_l in seq_var_ll]
    hover_data = ['Average CV', 'k-mer CVs']

    # # If mono-nucleotide k-mers, add k-mer CVs to hover box.
    # if kmer_size < 3:
    #     for idx, kmer in enumerate(seq_var_kmer_l):
    #         df[kmer] = [round(seq_var_l[2+idx], 5) for seq_var_l in seq_var_ll]
    #         hover_data.append(kmer)

    explained_variance = pca.explained_variance_ratio_ * 100
    color = 'Average CV'

    fig = px.scatter(
        df,
        x='PC1',
        y='PC2',
        color=color,
        # title='2D Visualization with Dataset IDs',
        labels={
            'PC1': f'PC1 ({explained_variance[0]:.2f}% variance)',
            'PC2': f'PC2 ({explained_variance[1]:.2f}% variance)'
        },
        hover_name='Dataset ID',
        color_continuous_scale=color_scale,
        hover_data=hover_data
    )

    # if kmer_size == 1:
    #     fig.update_traces(
    #         hovertemplate=(
    #             '<b>%{hovertext}</b><br>'
    #             'A CV: %{customdata[1]}<br>'
    #             'C CV: %{customdata[2]}<br>'
    #             'G CV: %{customdata[3]}<br>'
    #             'T CV: %{customdata[4]}<br>'
    #             'Average CV: %{customdata[0]}<extra></extra>'
    #         )
    #     )

    # if kmer_size < 3:
    #     # Construct k-mer CV string.
    #     kmer_cv_str = ""
    #     for index, kmer in enumerate(seq_var_kmer_l):
    #         kmer_cv_str += kmer + " CV: %{customdata[" + str(index+1) + "]}<br>"
    #     fig.update_traces(
    #         hovertemplate=(
    #             '<b>%{hovertext}</b><br><span style="font-family: \'Courier New\', monospace;">'
    #             + kmer_cv_str +
    #             'Average CV: %{customdata[0]}</span><extra></extra>'
    #         )
    #     )

    # else:

    fig.update_traces(
        hovertemplate=(
            '<b>%{hovertext}</b><br>%{customdata[1]}'
        )
    )

    fig.update_traces(marker=dict(size=12, line=dict(width=0.75, color='white')))

    fig.write_html(plot_out, full_html=full_html, include_plotlyjs=include_plotlyjs)


################################################################################

def create_seq_var_plot_plotly_v2(seq_var_ll, seq_var_kmer_l, plot_out,
                               kmer_size=1,
                               top_bottom_n=5,
                               remove_zero_val=True,
                               seq_len_stats_ll=False,
                               color_mode=1,
                               include_plotlyjs="cdn",
                               full_html=False):
    """
    Create sequence variation plot, use PCA to reduce dimensions of single k-mer 
    site reatios in the datasets.
    
    Formats:
    seq_var_ll = [dataset_id, avg_site_ratio. present_kmers_ratio, site_ratio_kmer1, site_ratio_kmer2 ... ]
    seq_var_kmer_l = [kmer1, kmer2, ...]

    remove_zero_val:
        If True, do not report zero site ratios/percentages in bottom rankings.

    """

    assert seq_var_ll, "seq_var_ll empty"
    assert seq_var_kmer_l, "seq_var_kmer_l empty"
    assert top_bottom_n <= 10, "top_bottom_n > 10"
    assert kmer_size < 6, "kmer_size > 5"

    color_scale = ['#c6dbef', '#9ecae1', '#6baed6', '#4292c6', '#2171b5', '#08519c', '#08306b']

    single_var_ll = [seq_var_l[3:] for seq_var_l in seq_var_ll]

    kmer_str_list = []
    
    # Construct strings for hover box.
    for idx1, single_var_l in enumerate(single_var_ll):

        # Average site ratio.
        dataset_id = seq_var_ll[idx1][0]
        avg_site_ratio = seq_var_ll[idx1][1]
        present_kmers_ratio = seq_var_ll[idx1][2]

        # Some sequence length stats.
        c_regions, median_len, min_len, max_len = False, False, False, False
        if seq_len_stats_ll:
            c_regions = seq_len_stats_ll[idx1][1]
            # mean_len = round(seq_len_stats_ll[idx1][2], 1)
            median_len = round(seq_len_stats_ll[idx1][3], 1)
            min_len = seq_len_stats_ll[idx1][6]
            max_len = seq_len_stats_ll[idx1][7]

        kmer2site_ratio_dic = {}
        for idx2, site_ratio in enumerate(single_var_l):
            kmer = seq_var_kmer_l[idx2]
            kmer2site_ratio_dic[kmer] = site_ratio

        kmer_str = '<span style="font-family: \'Courier New\', monospace;">'

        if kmer_size < 3:
            kmer_str += "%i-mer site %%:<br>" %(kmer_size)
            for kmer, site_ratio in sorted(kmer2site_ratio_dic.items(), key=lambda x: x[1], reverse=True):
                kmer_str += kmer + ": " + str(round(site_ratio*100, 3)) + "<br>"
        else:
            kmer_str += "Top %i %i-mer site %%:<br>" %(top_bottom_n, kmer_size)
            top_idx = 0
            for kmer, site_ratio in sorted(kmer2site_ratio_dic.items(), key=lambda x: x[1], reverse=True):
                # if remove_zero_val and site_ratio == 0:
                #     continue
                kmer_str += kmer + ": " + str(round(site_ratio*100, 3)) + "<br>"
                if top_idx == top_bottom_n-1:
                    break
                top_idx += 1
            kmer_str += "Bottom %i %i-mer site %%:<br>" %(top_bottom_n, kmer_size)
            bottom_idx = 0
            for kmer, site_ratio in sorted(kmer2site_ratio_dic.items(), key=lambda x: x[1], reverse=False):
                if remove_zero_val and site_ratio == 0:
                    continue
                kmer_str += kmer + ": " + str(round(site_ratio*100, 3)) + "<br>"
                if bottom_idx == top_bottom_n-1:
                    break
                bottom_idx += 1

        kmer_str += "Present %i-mer %%:      " %(kmer_size) + str(round(present_kmers_ratio*100, 3)) + "<br>"
        kmer_str += "Average %i-mer site %%: " %(kmer_size) + str(round(avg_site_ratio*100, 3))

        if c_regions:
            kmer_str += "<br># regions:     %i" %(c_regions)
        if median_len:
            kmer_str += "<br>Median length: %s" %(str(median_len))
        if min_len:
            kmer_str += "<br>Min length:    %i" %(min_len)
        if max_len:
            kmer_str += "<br>Max length:    %i" %(max_len)

        kmer_str += '</span><extra></extra>'

        kmer_str_list.append(kmer_str)

    exp_n_ratios = 4 ** kmer_size
    assert len(single_var_ll[0]) == exp_n_ratios, "unexpected number of site ratios in seq_var_ll"

    pca = PCA(n_components=2)
    data_2d_pca = pca.fit_transform(single_var_ll)
    df = pd.DataFrame(data_2d_pca, columns=['PC1', 'PC2'])

    df['Dataset ID'] = [seq_var_l[0] for seq_var_l in seq_var_ll]
    df['k-mer site %'] = kmer_str_list
    df['Mean site %'] = [seq_var_l[1]*100 for seq_var_l in seq_var_ll]
    df['Present k %'] = [seq_var_l[2]*100 for seq_var_l in seq_var_ll]
    hover_data = ['Mean site %', 'Present k %', 'k-mer site %']

    explained_variance = pca.explained_variance_ratio_ * 100

    if color_mode == 1:
        color = 'Mean site %'
    elif color_mode == 2:
        color = 'Present k %'
    else:
        assert False, "unexpected color_mode given"

    fig = px.scatter(
        df,
        x='PC1',
        y='PC2',
        color=color,
        labels={
            'PC1': f'PC1 ({explained_variance[0]:.2f}% variance)',
            'PC2': f'PC2 ({explained_variance[1]:.2f}% variance)'
        },
        hover_name='Dataset ID',
        color_continuous_scale=color_scale,
        hover_data=hover_data
    )

    fig.update_traces(
        hovertemplate=(
            '<b>%{hovertext}</b><br>%{customdata[2]}'
        )
    )

    fig.update_traces(marker=dict(size=12, line=dict(width=0.75, color='white')))

    fig.write_html(plot_out, full_html=full_html, include_plotlyjs=include_plotlyjs)


################################################################################

def create_eib_comp_plot_plotly(id2eib_stats_dic, id2eib_perc_dic, plot_out,
                                plot_3d=False,
                                id2hk_gene_stats_dic=False,
                                include_plotlyjs="cdn",
                                full_html=False):
    """
    Create exon-intron + border ratios comparison plot with plotly.

    Use PCA to reduce dimensions.
    
    Formats:
    id2eib_stats_dic[internal_id] = 
    [combined_id, ei_ol_stats.c_input_sites, exon_sites_ratio, intron_sites_ratio, us_ib_sites_ratio, 
    ds_ib_sites_ratio, us_ib_dist_sites_ratio, ds_ib_dist_sites_ratio, eib_sites_ratio, first_exon_sites_ratio, 
    last_exon_sites_ratio, single_exon_sites_ratio]
    id2eib_perc_dic[internal_id] = [exon_sites_perc, intron_sites_perc, us_ib_sites_perc, ds_ib_sites_perc, 
    us_ib_dist_sites_perc, ds_ib_dist_sites_perc, eib_sites_perc, first_exon_sites_perc, 
    last_exon_sites_perc, single_exon_sites_perc]

    """

    assert id2eib_stats_dic, "id2eib_stats_dic empty"
    assert id2eib_perc_dic, "id2eib_perc_dic empty"

    color_scale = ['#c6dbef', '#9ecae1', '#6baed6', '#4292c6', '#2171b5', '#08519c', '#08306b']

    eib_ratios_ll = []

    for internal_id in sorted(id2eib_stats_dic):
        combined_id = id2eib_stats_dic[internal_id][0]
        # c_input_sites = id2eib_stats_dic[internal_id][1]
        # exon_sites_ratio = id2eib_stats_dic[internal_id][2]
        # intron_sites_ratio = id2eib_stats_dic[internal_id][3]
        # us_ib_sites_ratio = id2eib_stats_dic[internal_id][4]
        # ds_ib_sites_ratio = id2eib_stats_dic[internal_id][5]
        # eib_sites_ratio = id2eib_stats_dic[internal_id][6]
        # Append elements 2-6.
        eib_ratios_ll.append(id2eib_stats_dic[internal_id][2:])

    if plot_3d:
        pca = PCA(n_components=3)
        data_3d_pca = pca.fit_transform(eib_ratios_ll)
        df = pd.DataFrame(data_3d_pca, columns=['PC1', 'PC2', 'PC3'])
    else:
        pca = PCA(n_components=2)
        data_2d_pca = pca.fit_transform(eib_ratios_ll)
        df = pd.DataFrame(data_2d_pca, columns=['PC1', 'PC2'])

    df['Dataset ID'] = [id2eib_stats_dic[internal_id][0] for internal_id in sorted(id2eib_stats_dic)]
    df['# input regions'] = [id2eib_stats_dic[internal_id][1] for internal_id in sorted(id2eib_stats_dic)]
    df['% exonic sites'] = [id2eib_perc_dic[internal_id][0] for internal_id in sorted(id2eib_perc_dic)]
    df['% intronic sites'] = [id2eib_perc_dic[internal_id][1] for internal_id in sorted(id2eib_perc_dic)]
    df['% us ib'] = [id2eib_perc_dic[internal_id][2] for internal_id in sorted(id2eib_perc_dic)]
    df['% ds ib'] = [id2eib_perc_dic[internal_id][3] for internal_id in sorted(id2eib_perc_dic)]
    df['% us dist ib'] = [id2eib_perc_dic[internal_id][4] for internal_id in sorted(id2eib_perc_dic)]
    df['% ds dist ib'] = [id2eib_perc_dic[internal_id][5] for internal_id in sorted(id2eib_perc_dic)]
    df['% eib'] = [id2eib_perc_dic[internal_id][6] for internal_id in sorted(id2eib_perc_dic)]
    df['% first exon'] = [id2eib_perc_dic[internal_id][7] for internal_id in sorted(id2eib_perc_dic)]
    df['% last exon'] = [id2eib_perc_dic[internal_id][8] for internal_id in sorted(id2eib_perc_dic)]
    df['% single exon'] = [id2eib_perc_dic[internal_id][9] for internal_id in sorted(id2eib_perc_dic)]
    # df['Exon sites perc'] = [id2eib_perc_dic[internal_id][0] for internal_id in sorted(id2eib_perc_dic)]
    # df['Intron sites perc'] = [id2eib_perc_dic[internal_id][1] for internal_id in sorted(id2eib_perc_dic)]
    # df['us ib sites perc'] = [id2eib_perc_dic[internal_id][2] for internal_id in sorted(id2eib_perc_dic)]
    # df['ds ib sites perc'] = [id2eib_perc_dic[internal_id][3] for internal_id in sorted(id2eib_perc_dic)]
    # df['eib sites perc'] = [id2eib_perc_dic[internal_id][4] for internal_id in sorted(id2eib_perc_dic)]

    # 11 features.
    hover_data= ['# input regions', '% exonic sites', '% intronic sites', '% us ib', '% ds ib', '% us dist ib', '% ds dist ib', '% eib', '% first exon', '% last exon', '% single exon']
    color = '% intronic sites'

    if id2hk_gene_stats_dic:
        # Genes (actually transcripts).
        df['# all genes'] = [id2hk_gene_stats_dic[internal_id][0] for internal_id in sorted(id2hk_gene_stats_dic)]
        df['# HK genes'] = [id2hk_gene_stats_dic[internal_id][1] for internal_id in sorted(id2hk_gene_stats_dic)]
        df['% HK genes'] = [id2hk_gene_stats_dic[internal_id][2] for internal_id in sorted(id2hk_gene_stats_dic)]
        hover_data += ['# all genes', '# HK genes', '% HK genes']
        color = '% HK genes'

    explained_variance = pca.explained_variance_ratio_ * 100

    if plot_3d:

        fig = px.scatter_3d(
            df,  # Use the DataFrame directly
            x='PC1',
            y='PC2',
            z='PC3',
            color=color,
            # title='2D Visualization with Dataset IDs',
            labels={
                'PC1': f'PC1 ({explained_variance[0]:.2f}% variance)',
                'PC2': f'PC2 ({explained_variance[1]:.2f}% variance)',
                'PC3': f'PC3 ({explained_variance[2]:.2f}% variance)'
            },
            hover_name='Dataset ID',
            color_continuous_scale=color_scale,
            hover_data=hover_data
        )

        if id2hk_gene_stats_dic:

            fig.update_traces(
                hovertemplate=(
                    '<b>%{hovertext}</b><br>'
                    'Input regions (#): %{customdata[0]}<br>'
                    'Exon: %{customdata[1]}%<br>'
                    'Intron: %{customdata[2]}%<br>'
                    'US intron border: %{customdata[3]}%<br>'
                    'DS intron border: %{customdata[4]}%<br>'
                    'US distant intron: %{customdata[5]}%<br>'
                    'DS distant intron: %{customdata[6]}%<br>'
                    'Exon-intron border: %{customdata[7]}%<extra></extra>'
                    'First exon: %{customdata[8]}%<br>'
                    'Last exon: %{customdata[9]}%<br>'
                    'Single exon: %{customdata[10]}%<br>'
                    'Occupied genes (#): %{customdata[11]}<br>'
                    'Occupied HK genes (#): %{customdata[12]}<br>'
                    'Occupied HK genes (%): %{customdata[13]}<extra></extra>'
                )
            )

        else:

            fig.update_traces(
                hovertemplate=(
                    '<b>%{hovertext}</b><br>'
                    'Input regions (#): %{customdata[0]}<br>'
                    'Exon: %{customdata[1]}%<br>'
                    'Intron: %{customdata[2]}%<br>'
                    'US intron border: %{customdata[3]}%<br>'
                    'DS intron border: %{customdata[4]}%<br>'
                    'US distant intron: %{customdata[5]}%<br>'
                    'DS distant intron: %{customdata[6]}%<br>'
                    'Exon-intron border: %{customdata[7]}%<extra></extra>'
                    'First exon: %{customdata[8]}%<br>'
                    'Last exon: %{customdata[9]}%<br>'
                    'Single exon: %{customdata[10]}%<extra></extra>'
                )
            )

        fig.update_traces(marker=dict(size=3, line=dict(width=0.5, color='white')))

    else:

        fig = px.scatter(
            df,
            x='PC1',
            y='PC2',
            color=color,
            # title='2D Visualization with Dataset IDs',
            labels={
                'PC1': f'PC1 ({explained_variance[0]:.2f}% variance)',
                'PC2': f'PC2 ({explained_variance[1]:.2f}% variance)'
            },
            hover_name='Dataset ID',
            color_continuous_scale=color_scale,
            hover_data=hover_data
        )

        if id2hk_gene_stats_dic:

            fig.update_traces(
                hovertemplate=(
                    '<b>%{hovertext}</b><br>'
                    'Input regions (#): %{customdata[0]}<br>'
                    'Exon: %{customdata[1]}%<br>'
                    'Intron: %{customdata[2]}%<br>'
                    'US intron border: %{customdata[3]}%<br>'
                    'DS intron border: %{customdata[4]}%<br>'
                    'US distant intron: %{customdata[5]}%<br>'
                    'DS distant intron: %{customdata[6]}%<br>'
                    'Exon-intron border: %{customdata[7]}%<br>'
                    'First exon: %{customdata[8]}%<br>'
                    'Last exon: %{customdata[9]}%<br>'
                    'Single exon: %{customdata[10]}%<br>'
                    'Occupied genes (#): %{customdata[11]}<br>'
                    'Occupied HK genes (#): %{customdata[12]}<br>'
                    'Occupied HK genes (%): %{customdata[13]}<extra></extra>'
                )
            )

        else:

            fig.update_traces(
                hovertemplate=(
                    '<b>%{hovertext}</b><br>'
                    'Input regions (#): %{customdata[0]}<br>'
                    'Exon: %{customdata[1]}%<br>'
                    'Intron: %{customdata[2]}%<br>'
                    'US intron border: %{customdata[3]}%<br>'
                    'DS intron border: %{customdata[4]}%<br>'
                    'US distant intron: %{customdata[5]}%<br>'
                    'DS distant intron: %{customdata[6]}%<br>'
                    'Exon-intron border: %{customdata[7]}%<br>'
                    'First exon: %{customdata[8]}%<br>'
                    'Last exon: %{customdata[9]}%<br>'
                    'Single exon: %{customdata[10]}%<extra></extra>'
                )
            )


        fig.update_traces(marker=dict(size=12, line=dict(width=0.75, color='white')))

    fig.write_html(plot_out, full_html=full_html, include_plotlyjs=include_plotlyjs)


################################################################################

def create_pca_reg_occ_plot_plotly(id2occ_list_dic, id2infos_dic, 
                                   plot_out,
                                   sparse_pca=False,
                                   id2hk_gene_stats_dic=False,
                                   add_motif_db_info=False,
                                   include_plotlyjs="cdn",
                                   full_html=False):
    """
    Create plotly 3d scatter plot of PCA reduced gene/ transcript occupancies.

    id2occ_list_dic:
        Format: internal_id -> [0, 1, 0, 0, 1] (transcript occupancy list).
    
    id2infos_dic:
        id2infos_dic[internal_id] = [rbp_id, data_id, method_id, motif_db_str, bed_file_path]
    
    """

    assert id2occ_list_dic, "id2occ_list_dic empty"

    occ_ll = []
    dataset_ids_list = []
    c_one_labels_list = []
    perc_one_labels_list = []
    c_total_number_genes = 0

    c_all_tr_list = []
    c_hk_tr_list = []
    perc_hk_tr_list = []

    for internal_id in sorted(id2occ_list_dic):

        if not c_total_number_genes:
            c_total_number_genes = len(id2occ_list_dic[internal_id])
        else:
            assert c_total_number_genes == len(id2occ_list_dic[internal_id]), "inconsistent number of genes/transcripts in occupancy lists"

        # Count number of 1's in occupancy list.
        c_one_labels = id2occ_list_dic[internal_id].count(1)
        c_one_labels_list.append(c_one_labels)
        # Percentage of 1's in occupancy list. Round to two decimal places.
        perc_one_labels = round((c_one_labels / c_total_number_genes) * 100, 2)
        perc_one_labels_list.append(perc_one_labels)

        rbp_id = id2infos_dic[internal_id][0]
        data_id = id2infos_dic[internal_id][1]
        method_id = id2infos_dic[internal_id][2]
        database_id = id2infos_dic[internal_id][3]

        combined_id = rbp_id + "," + method_id + "," + data_id
        if add_motif_db_info:
            combined_id = rbp_id + "," + database_id + "," + method_id + "," + data_id

        occ_ll.append(id2occ_list_dic[internal_id])
        dataset_ids_list.append(combined_id)

        if id2hk_gene_stats_dic:
            c_all_tr = id2hk_gene_stats_dic[internal_id][0]
            c_hk_tr = id2hk_gene_stats_dic[internal_id][1]
            perc_hk_tr = id2hk_gene_stats_dic[internal_id][2]
            c_all_tr_list.append(c_all_tr)
            c_hk_tr_list.append(c_hk_tr)
            perc_hk_tr_list.append(perc_hk_tr)

        # print("combined_id:", combined_id)
        # print(id2occ_list_dic[internal_id])

    hover_data = ['# occupied genes', '% occupied genes']
    color = '% occupied genes'

    if id2hk_gene_stats_dic:
        hover_data += ['# all transcripts', '# hk transcripts', '% HK genes']
        color = '% HK genes'

    color_scale = ['#c6dbef', '#9ecae1', '#6baed6', '#4292c6', '#2171b5', '#08519c', '#08306b']

    if sparse_pca:

        from sklearn.decomposition import SparsePCA

        pca = SparsePCA(n_components=3, random_state=0, ridge_alpha=0.01)
        data_3d_pca = pca.fit_transform(occ_ll)

        df = pd.DataFrame(data_3d_pca, columns=['PC1', 'PC2', 'PC3'])
        df['Dataset ID'] = dataset_ids_list
        df['# occupied genes'] = c_one_labels_list
        df['% occupied genes'] = perc_one_labels_list

        fig = px.scatter_3d(
            df,  # Use the DataFrame directly
            x='PC1',
            y='PC2',
            z='PC3',
            title='3D Visualization with Dataset IDs',
            hover_name='Dataset ID'
        )

    else:
        pca = PCA(n_components=3)  # Reduce data to 3 dimensions.
        data_3d_pca = pca.fit_transform(occ_ll)

        df = pd.DataFrame(data_3d_pca, columns=['PC1', 'PC2', 'PC3'])
        df['Dataset ID'] = dataset_ids_list
        df['# occupied genes'] = c_one_labels_list
        df['% occupied genes'] = perc_one_labels_list
        if id2hk_gene_stats_dic:
            df['# all transcripts'] = c_all_tr_list
            df['# hk transcripts'] = c_hk_tr_list
            df['% HK genes'] = perc_hk_tr_list

        explained_variance = pca.explained_variance_ratio_ * 100

        fig = px.scatter_3d(
            df,  # Use the DataFrame directly
            x='PC1',
            y='PC2',
            z='PC3',
            color=color,
            title='3D Visualization with Dataset IDs',
            labels={
                'PC1': f'PC1 ({explained_variance[0]:.2f}% variance)',
                'PC2': f'PC2 ({explained_variance[1]:.2f}% variance)',
                'PC3': f'PC3 ({explained_variance[2]:.2f}% variance)'
            },
            hover_name='Dataset ID',
            color_continuous_scale=color_scale,
            hover_data=hover_data
        )

    if id2hk_gene_stats_dic:

        fig.update_traces(
            hovertemplate='<b>%{hovertext}</b><br>Occupied genes (#): %{customdata[0]}<br>Occupied genes (%): %{customdata[1]}<br>Occupied HK genes (#): %{customdata[3]}<br>Occupied HK genes (%): %{customdata[4]}<extra></extra>'
        )

    else:

        fig.update_traces(
            hovertemplate='<b>%{hovertext}</b><br>Occupied genes (#): %{customdata[0]}<br>Occupied genes (%): %{customdata[1]}<extra></extra>'
        )


    fig.update_scenes(aspectmode='cube')

    fig.update_traces(marker=dict(size=3, line=dict(width=0.5, color='white')))
    fig.update_layout(margin=dict(l=0, r=0, b=0, t=0))
    fig.write_html(plot_out, full_html=full_html, include_plotlyjs=include_plotlyjs)


################################################################################

def calculate_k_site_ratios(reg2seq_dic, k=1,
                            nucleotides=['A', 'C', 'G', 'T'],
                            only_observed=True):
    """
    Calculate k-mer site ratios (i.e., ratios of sequences where k-mers 
    are present).
    Return single k-mer site ratios dictionary, the average site 
    ratio value over all k-mers, and the number of k-mers present ratio.

    only_observed:
        If True, calculate the average site ratio only for observed k-mers.

    >>> reg2seq_dic = {1: 'ACGT', 2: 'AC'}
    >>> kmer2sr_dic, avg_sr, pk_rat = calculate_k_site_ratios(reg2seq_dic, k=1)    
    >>> kmer2sr_dic, avg_sr, pk_rat
    ({'A': 1.0, 'C': 1.0, 'G': 0.5, 'T': 0.5}, 0.75, 1.0)
    >>> reg2seq_dic = {1: 'ACG'}
    >>> kmer2sr_dic, avg_sr, pk_rat = calculate_k_site_ratios(reg2seq_dic, k=1, only_observed=False)
    >>> kmer2sr_dic, avg_sr, pk_rat
    ({'A': 1.0, 'C': 1.0, 'G': 1.0, 'T': 0.0}, 0.75, 0.75)
    >>> kmer2sr_dic, avg_sr, pk_rat = calculate_k_site_ratios(reg2seq_dic, k=1, only_observed=True)
    >>> kmer2sr_dic, avg_sr, pk_rat
    ({'A': 1.0, 'C': 1.0, 'G': 1.0, 'T': 0.0}, 1.0, 0.75)
    >>> reg2seq_dic = {1: 'ACGT', 2: 'CC'}
    >>> kmer2sr_dic, avg_sr, pk_rat = calculate_k_site_ratios(reg2seq_dic, k=2, only_observed=True)
    >>> kmer2sr_dic, avg_sr, pk_rat
    ({'AA': 0.0, 'AC': 0.5, 'AG': 0.0, 'AT': 0.0, 'CA': 0.0, 'CC': 0.5, 'CG': 0.5, 'CT': 0.0, 'GA': 0.0, 'GC': 0.0, 'GG': 0.0, 'GT': 0.5, 'TA': 0.0, 'TC': 0.0, 'TG': 0.0, 'TT': 0.0}, 0.5, 0.25)
    >>> kmer2sr_dic, avg_sr, pk_rat = calculate_k_site_ratios(reg2seq_dic, k=2, only_observed=False)
    >>> kmer2sr_dic, avg_sr, pk_rat
    ({'AA': 0.0, 'AC': 0.5, 'AG': 0.0, 'AT': 0.0, 'CA': 0.0, 'CC': 0.5, 'CG': 0.5, 'CT': 0.0, 'GA': 0.0, 'GC': 0.0, 'GG': 0.0, 'GT': 0.5, 'TA': 0.0, 'TC': 0.0, 'TG': 0.0, 'TT': 0.0}, 0.125, 0.25)

    """
    assert reg2seq_dic, "reg2seq_dic empty"

    # Generate all possible k-mers for specified k.
    kmers = [''.join(p) for p in product(nucleotides, repeat=k)]
    # Calculate number of expected k-mers.
    exp_num_kmers = len(nucleotides) ** k

    # Number of sequences containing each k-mer.
    kmer2sitec_dic = {}
    for kmer in kmers:
        kmer2sitec_dic[kmer] = 0

    # Count k-mers for each sequence.
    seq_len_list = []
    for reg_id in sorted(reg2seq_dic):
        seq = reg2seq_dic[reg_id]
        # Number of k-mers in the sequence.
        length = len(seq) - k + 1
        kmer_c = 0
        if length > 0:
            # Count each k-mer in the sequence.
            counts = defaultdict(int)
            for i in range(length):
                kmer = seq[i:i + k]
                if kmer in kmers:  # Only valid k-mers.
                    counts[kmer] += 1

            # Count number of sequences containing each k-mer.
            for kmer in counts:
                if counts[kmer] > 0:
                    kmer2sitec_dic[kmer] += 1

    kmer2site_ratio_dic = {}
    site_ratios_list = []
    c_seqs = len(reg2seq_dic)
    c_zero_count_kmers = 0
    for kmer, sitec in kmer2sitec_dic.items():
        site_ratio = sitec / c_seqs
        if only_observed:
            if sitec > 0:  # Calculate average ratio only using k-mers present in sequences.
                site_ratios_list.append(site_ratio)
        else:
            site_ratios_list.append(site_ratio)
        kmer2site_ratio_dic[kmer] = site_ratio
        if sitec == 0:
            c_zero_count_kmers += 1

    assert len(kmer2site_ratio_dic) == exp_num_kmers, "unexpected number of k-mers (expected: %i, got: %i)" % (exp_num_kmers, len(kmer2site_ratio_dic))

    c_count_kmers = exp_num_kmers - c_zero_count_kmers
    # Ratio of k-mers present in sequences.
    present_kmers_ratio = c_count_kmers / exp_num_kmers

    # Get average site ratio over all k-mers.
    average_site_ratio = statistics.mean(site_ratios_list) if site_ratios_list else 0

    return kmer2site_ratio_dic, average_site_ratio, present_kmers_ratio


################################################################################

def calculate_k_nucleotide_cv(reg2seq_dic, k=1,
                              nucleotides=['A', 'C', 'G', 'T'],
                              kmer2stats_dic=None,
                              reg2sc_dic=False,
                              only_observed=True):
    """
    Calculate the coefficient of variation (CV) for each k-mer in the sequences.
    Return the single k-mer CVs and the average CV over all k-mers.

    only_observed:
        If True, calculate the average CV only for observed k-mers.
        If False, calculate the average CV for all possible k-mers.
    kmer2stats_dic:
        If supplied, store k-mer total count and ratio.
        kmer2stats_dic[kmer] = [count, ratio, site_ratio, corr]
        corr if reg2sc_dic is supplied.
    
    """
    # Generate all possible k-mers for specified k.
    kmers = [''.join(p) for p in product(nucleotides, repeat=k)]
    # Calculate number of expected k-mers.
    exp_num_kmers = len(nucleotides) ** k

    kmer_ratios = {kmer: [] for kmer in kmers}
    
    observed_kmer_dic = {}
    total_c = 0
    kmer2count_dic = {}
    for kmer in kmers:
        kmer2count_dic[kmer] = 0
    # Number of sequences containing each k-mer.
    kmer2sitec_dic = {}
    for kmer in kmers:
        kmer2sitec_dic[kmer] = 0

    # Calculate k-mer ratios for each sequence.
    reg_ids_list = []
    for reg_id in sorted(reg2seq_dic):
        seq = reg2seq_dic[reg_id]
        reg_ids_list.append(reg_id)
        length = len(seq) - k + 1  # Number of k-mers in the sequence.
        kmer_c = 0
        if length > 0:
            # Count each k-mer in the sequence.
            counts = defaultdict(int)
            for i in range(length):
                kmer = seq[i:i + k]
                if kmer in kmer_ratios:  # Only valid k-mers.
                    counts[kmer] += 1
                    kmer2count_dic[kmer] += 1
                    observed_kmer_dic[kmer] = 1
                    kmer_c += 1
                    total_c += 1
            
            # Store the ratio of each k-mer.
            for kmer in kmers:
                kmer_ratios[kmer].append(counts[kmer] / kmer_c)

            # Count number of sequences containing each k-mer.
            for kmer in counts:
                if counts[kmer] > 0:
                    kmer2sitec_dic[kmer] += 1
        else:
            # Store zeros for all k-mers.
            for kmer in kmers:
                kmer_ratios[kmer].append(0)

    for kmer in kmer_ratios:
        length = len(kmer_ratios[kmer])
        assert length == len(reg_ids_list), "unexpected number of ratios for k-mer %s (# ratios: %i, # region IDs: %i)" %(kmer, length, len(reg_ids_list))

    if kmer2stats_dic is not None:
        for kmer, count in kmer2count_dic.items():
            kmer_ratio = count / total_c
            site_ratio = (kmer2sitec_dic[kmer] / len(reg2seq_dic)) if len(reg2seq_dic) else 0
            kmer2stats_dic[kmer] = [count, kmer_ratio, site_ratio, 0]

    if reg2sc_dic is not None:
        scores_list = []
        for reg_id in reg_ids_list:
            assert reg_id in reg2sc_dic, "reg_id not in reg2sc_dic"
            scores_list.append(reg2sc_dic[reg_id])
        # Calculate the correlation between k-mer ratios and scores.
        scores = pd.Series(scores_list, index=reg_ids_list)
        df = pd.DataFrame(kmer_ratios, index=reg_ids_list)
        corrs = df.corrwith(scores, method='spearman')
        # Make dictionary with k-mer -> correlation.
        kmer2corr_dic = corrs.to_dict()
        for kmer in kmer2corr_dic:
            kmer2stats_dic[kmer][3] = round(kmer2corr_dic[kmer], 5)

    # Calculate the coefficient of variation (CV) for each k-mer.
    kmer_cv = {}
    for kmer, ratios in kmer_ratios.items():
        if ratios:  # Only for k-mers present in the sequences.
            mean = np.mean(ratios)
            std_dev = np.std(ratios)
            cv = (std_dev / mean) if mean != 0 else 0.0
            kmer_cv[kmer] = cv

    assert len(kmer_cv) == exp_num_kmers, "unexpected number of k-mers (expected: %i, got: %i)" % (exp_num_kmers, len(kmer_cv))

    if only_observed:
        # Calculate the average CV for observed k-mers only using observed_kmer_dic.
        observed_cvs = [cv for kmer, cv in kmer_cv.items() if observed_kmer_dic.get(kmer, 0)]
        average_cv = np.mean(observed_cvs) if observed_cvs else 0
    else:
        # Calculate the average CV for all k-mers.
        average_cv = np.mean(list(kmer_cv.values())) if kmer_cv else 0

    # Convert kmer_cv from numpy float to python float.
    kmer_cv = {k: float(v) for k, v in kmer_cv.items()}

    return kmer_cv, average_cv


################################################################################

def create_kmer_comp_plot_plotly(dataset_ids_list, kmer_list, kmer_freqs_ll,
                                 seq_len_stats_ll, plot_out, 
                                 seq_feat_ll=False,
                                 include_plotlyjs="cdn",
                                 full_html=False):
    
    """
    Create plotly 3d scatter plot of PCA reduced k-mer frequencies.

    seq_len_stats_ll + seq_feat_ll in same order as dataset_ids_list.

    kmer_list:
        1d list of k-mers, corresponding in order to k-mer frequency vectors in
        kmer_freqs_ll.

    """

    assert seq_len_stats_ll, "seq_len_stats_ll empty"

    kmer_len = len(kmer_list[0])
    n_top_kmers = 10
    top_kmer_str_list = []
    top_kmer_str_header = "Top " + str(n_top_kmers) + " k-mers"
    for kmer_freqs_l in kmer_freqs_ll:
        kmer2freq_dic = {}
        for idx, kmer_freq in enumerate(kmer_freqs_l):
            kmer = kmer_list[idx]
            kmer2freq_dic[kmer] = kmer_freq
        # Go over sorted kmer2freq_dic descending by value, outputting top n_top_kmers entries.
        top_kmer_str = "<br>"
        for kmer, freq in sorted(kmer2freq_dic.items(), key=lambda x: x[1], reverse=True)[:n_top_kmers]:
            # Convert frequency to percentage.
            perc = round(freq * 100, 2)
            top_kmer_str += kmer + ": " + str(perc) + "%<br>"
        top_kmer_str_list.append(top_kmer_str)

    # No .np conversion needed.
    # data = np.array(kmer_freqs_ll)
    pca = PCA(n_components=3)  # Reduce data to 3 dimensions.
    data_3d_pca = pca.fit_transform(kmer_freqs_ll)
    explained_variance = pca.explained_variance_ratio_ * 100

    df = pd.DataFrame(data_3d_pca, columns=['PC1', 'PC2', 'PC3'])
    df['Dataset ID'] = dataset_ids_list
    df[top_kmer_str_header] = top_kmer_str_list
    df['# input regions'] = [seq_len_stats[1] for seq_len_stats in seq_len_stats_ll]
    hover_data = [top_kmer_str_header, '# input regions']
    color = '# input regions'
    
    if seq_feat_ll:
        df['Mean complexity'] = [seq_feat_l[1] for seq_feat_l in seq_feat_ll]

        mono_nts_str_list = []
        for seq_feat_l in seq_feat_ll:
            # mono_nts_str = "A: %s%%<br>C: %s%%<br>G: %s%%<br>T: %s%%" %(seq_feat_l[1], seq_feat_l[2], seq_feat_l[3], seq_feat_l[4])
            mono_nts_str = "A: %s%%, C: %s%%, G: %s%%, T: %s%%" %(seq_feat_l[2], seq_feat_l[3], seq_feat_l[4], seq_feat_l[5])
            mono_nts_str_list.append(mono_nts_str)

        df['Mono-nucleotide percentages'] = mono_nts_str_list

        hover_data.append('Mean complexity')
        hover_data.append('Mono-nucleotide percentages')
        color = 'Mean complexity'

    color_scale = ['#c6dbef', '#9ecae1', '#6baed6', '#4292c6', '#2171b5', '#08519c', '#08306b']
    # color_scale = ['#9ecae1', '#6baed6', '#4292c6', '#2171b5', '#08519c', '#08306b']

    fig = px.scatter_3d(
        df,
        x='PC1',
        y='PC2',
        z='PC3',
        color=color,
        title='3D Visualization with Dataset IDs',
        labels={
            'PC1': f'PC1 ({explained_variance[0]:.2f}% variance)',
            'PC2': f'PC2 ({explained_variance[1]:.2f}% variance)',
            'PC3': f'PC3 ({explained_variance[2]:.2f}% variance)'
        },
        hover_name='Dataset ID',
        color_continuous_scale=color_scale,
        hover_data=hover_data
    )

    if seq_feat_ll:
        fig.update_traces(
            hovertemplate = (
                '<b>%{hovertext}</b><br>Top ' + str(n_top_kmers) + ' ' + str(kmer_len) + 
                '-mer percentages:<span style="font-family: \'Courier New\', monospace;">%{customdata[0]}</span>Input regions (#):<br>%{customdata[1]}<br>Mean sequence complexity:<br>%{customdata[2]}<br>Mono-nucleotide percentages:<br>%{customdata[3]}<extra></extra>'
            )
        )

    else:
        fig.update_traces(
            hovertemplate = (
                '<b>%{hovertext}</b><br>Top ' + str(n_top_kmers) + ' ' + str(kmer_len) + 
                '-mer percentages:<span style="font-family: \'Courier New\', monospace;">%{customdata[0]}</span>Input regions (#):<br>%{customdata[1]}<extra></extra>'
            )
        )


    # fig = px.scatter_3d(
    #     df,
    #     x='PC1',
    #     y='PC2',
    #     z='PC3',
    #     title='3D Visualization with Dataset IDs',
    #     labels={
    #         'PC1': f'PC1 ({explained_variance[0]:.2f}% variance)',
    #         'PC2': f'PC2 ({explained_variance[1]:.2f}% variance)',
    #         'PC3': f'PC3 ({explained_variance[2]:.2f}% variance)'
    #     },
    #     hover_name='Dataset ID',
    #     hover_data=[top_kmer_str_header]
    # )

    # fig.update_traces(
    #     hovertemplate='<b>%{hovertext}</b><br>Top ' + str(n_top_kmers) + ' ' + str(kmer_len) + '-mer percentages: %{customdata[0]}<extra></extra>'
    # )

    # fig.update_traces(hovertemplate='%{hovertext}')  # This sets the hover template to only show the hover text
    fig.update_scenes(aspectmode='cube')
    # fig.update_traces(marker=dict(size=3))
    fig.update_traces(marker=dict(size=3, line=dict(width=0.5, color='white')))
    fig.update_layout(margin=dict(l=0, r=0, b=0, t=0))
    fig.write_html(plot_out, full_html=full_html, include_plotlyjs=include_plotlyjs)


################################################################################

def get_gene_occ_cooc_tables(id2occ_list_dic, id2infos_dic,
                             optimal_ordering=True):
    """
    Calculate co-occurrence matrices for plotting cosine similarities between datasets
    (more precisely between their gene list binary vectors).
    So a high cosine similarity indicates that the gene lists of two datasets are
    similar / tend to have more co-occurrence or a lack of occurrence at same genes.

    id2occ_list_dic:
        Format: internal_id -> [0, 1, 0, 0, 1] (transcript/gene region occupancy list).
    id2infos_dic:
        id2infos_dic[internal_id] = [rbp_id, data_id, method_id, motif_db_str, bed_file_path]

    optimal_ordering:
        If True, the linkage matrix will be reordered so that the distance between successive 
        leaves is minimal. This results in a more intuitive tree structure when the data are 
        visualized. defaults to False, because this algorithm can be slow, particularly on 
        large datasets. See also the optimal_leaf_ordering function.
        https://docs.scipy.org/doc/scipy/reference/generated/scipy.cluster.hierarchy.linkage.html
        
    """

    assert id2occ_list_dic, "id2lst_dic empty"
    assert id2infos_dic, "id2infos_dic empty"

    from sklearn.metrics.pairwise import cosine_similarity
    from sklearn.metrics import pairwise_distances
    from scipy.cluster.hierarchy import linkage, leaves_list

    set_internal_ids_list = list(id2occ_list_dic.keys())
    set_plot_ids_list = []
    seen_ids_dic = {}
    for internal_id in set_internal_ids_list:
        set_plot_id = id2infos_dic[internal_id][0] + "," + id2infos_dic[internal_id][2] + "," + id2infos_dic[internal_id][1]
        if set_plot_id in seen_ids_dic:
            assert False, "duplicate dataset ID encountered (%s), which is incompatible with gene region occupancy heatmap plot. Please provide a unique combination of selected RBP ID, method ID, and data ID for input dataset or disable plot by setting --no-occ-heatmap" %(set_plot_id)
        seen_ids_dic[set_plot_id] = internal_id
        set_plot_ids_list.append(set_plot_id)

    len_set_list = len(set_plot_ids_list)

    # Make df out of gene region occupancy lists.
    df = pd.DataFrame.from_dict(id2occ_list_dic, orient='index')

    # Calculate the cosine similarity.
    similarity_matrix = cosine_similarity(df)

    # Convert the similarity matrix to a DataFrame
    df_sim = pd.DataFrame(similarity_matrix, index=df.index, columns=df.index)

    # Calculate a cosine distance matrix for clustering.
    dist_matrix = pairwise_distances(df, metric='cosine')

    # Perform hierarchical clustering.
    dist_matrix_condensed = dist_matrix[np.triu_indices(dist_matrix.shape[0], k=1)]
    linkage_matrix = linkage(dist_matrix_condensed, method='average', optimal_ordering=optimal_ordering)

    # Calculate leaves order.
    leaves_order = leaves_list(linkage_matrix)

    # print("df_sim:")
    # print(df_sim)

    # Assign new index IDs to df_sim, namely set_plot_ids_list.
    df_sim.index = set_plot_ids_list
    df_sim.columns = set_plot_ids_list

    # for idx, set_internal_id in enumerate(set_internal_ids_list):
    #     print(set_internal_id, set_plot_ids_list[idx])

    # print("df_sim new IDs:")
    # print(df_sim)

    # print("leaves_order:")
    # print(leaves_order)

    df_sim_ordered = df_sim.iloc[leaves_order, :].iloc[:, leaves_order]

    # Get newly ordered IDs list.
    set_plot_ids_list_ordered = [set_plot_ids_list[i] for i in leaves_order]

    # print("set_plot_ids_list_ordered:")
    # print(set_plot_ids_list_ordered)

    # print("df_sim_ordered:")
    # print(df_sim_ordered)

    # plot = px.imshow(df_sim_ordered)
    # plot_out = "test_clustered_heatmap.html"
    # plot.write_html(plot_out)

    # Create 3d list storing additional infos for co-occurrence heatmap plot.
    gene_cooc_lll = []

    for set_id in set_plot_ids_list:
        gene_cooc_lll.append([]*len_set_list)

    for i in range(len_set_list):
        for j in range(len_set_list):
            gene_cooc_lll[i].append(["-", "-", "-", "-"])

    seen_pairs_dic = {}

    for i,set_id_i in enumerate(set_plot_ids_list_ordered):
        for j,set_id_j in enumerate(set_plot_ids_list_ordered):

            pair = [i, j]
            pair.sort()
            pair_str = str(pair[0]) + "," + str(pair[1])
            # id_pair = [set_id_i, set_id_j]
            # id_pair.sort()
            # id_str = id_pair[0] + " " + id_pair[1]
            if pair_str in seen_pairs_dic:
                continue
            seen_pairs_dic[pair_str] = 1

            if i == j:  # Diagonal, i.e. set_id_i == set_id_j.
                gene_cooc_lll[i][j][0] = "1.0"
                gene_cooc_lll[i][j][1] = set_id_i
                gene_cooc_lll[i][j][2] = set_id_j
                set_12 = id2occ_list_dic[seen_ids_dic[set_id_i]]
                a = 0
                d = 0
                for occ in set_12:
                    if occ == 1:
                        a += 1
                    else:
                        d += 1
                table = [[a, 0], [0, d]]
                table_str = str(table)
                gene_cooc_lll[i][j][3] = table_str

            else:
                # Cosine similarity.
                cosine_sim = df_sim_ordered.loc[set_id_i, set_id_j]
                gene_cooc_lll[i][j][0] = str(cosine_sim)
                gene_cooc_lll[j][i][0] = str(cosine_sim)
                # Set 1+2 IDs.
                gene_cooc_lll[i][j][1] = set_id_i
                gene_cooc_lll[j][i][1] = set_id_i
                gene_cooc_lll[i][j][2] = set_id_j
                gene_cooc_lll[j][i][2] = set_id_j

                set1 = id2occ_list_dic[seen_ids_dic[set_id_i]]
                set2 = id2occ_list_dic[seen_ids_dic[set_id_j]]

                a = sum([1 for x, y in zip(set1, set2) if x == 1 and y == 1])
                b = sum([1 for x, y in zip(set1, set2) if x == 1 and y == 0])
                c = sum([1 for x, y in zip(set1, set2) if x == 0 and y == 1])
                d = sum([1 for x, y in zip(set1, set2) if x == 0 and y == 0])

                table = [[a, b], [c, d]]
                table_str = str(table)

                gene_cooc_lll[i][j][3] = table_str
                gene_cooc_lll[j][i][3] = table_str

    return df_sim_ordered, gene_cooc_lll


################################################################################

def create_gene_occ_cooc_plot_plotly(df_gene_occ, gene_cooc_lll, plot_out,
                                     include_plotlyjs="cdn",
                                     full_html=False):
    """
    Plot gene occupancy co-occurrences between batch input datasets 
    as heat map with plotly.

    """

    min_val = 0.01
    max_val = 1
    blues_colors = px.colors.sequential.Blues
    blues_colors_hex = [rgb_to_hex(color) for color in blues_colors]
    color_scale = create_color_scale(min_val, max_val, blues_colors_hex)
    color_scale.insert(0, [0, "white"])
    zmin = min(0, df_gene_occ.min().min())
    zmax = max(1, df_gene_occ.max().max())

    plot = px.imshow(df_gene_occ, color_continuous_scale=color_scale, zmin=zmin, zmax=zmax)

    plot.update(data=[{'customdata': gene_cooc_lll,
                       'hovertemplate': '1) Set1: %{customdata[1]}<br>2) Set2: %{customdata[2]}<br>3) cosine similarity: %{customdata[0]}<br>4) Counts: %{customdata[3]}<extra></extra>'}])
    plot.update_layout(plot_bgcolor='white')
    plot.update_layout(
        annotations=[
            dict(
                x=1.05,
                y=1.03,
                align="right",
                valign="top",
                text="Similarity",
                showarrow=False,
                xref="paper",
                yref="paper",
                xanchor="center",
                yanchor="top",
                font=dict(size=14)
            )
        ]
    )
    plot.write_html(plot_out,
                    full_html=full_html,
                    include_plotlyjs=include_plotlyjs)


################################################################################

def goa_generate_html_report(args, goa_results_df, 
                             goa_stats_dic, benchlib_path,
                             html_report_out="goa_results.rbpbench_goa.html"):
    """
    Generate GOA results HTML report (rbpbench goa mode).

    """
    out_folder = args.out_folder
    if args.plot_abs_paths:
        out_folder = os.path.abspath(out_folder)

    # Version string.
    version_str = "v" + args.version

    html_out = out_folder + "/" + "goa_results.rbpbench_goa.html"
    md_out = out_folder + "/" + "goa_results.rbpbench_goa.md"
    if html_report_out:
        html_out = html_report_out

    """
    Setup sorttable.js to make tables in HTML sortable.

    """
    sorttable_js_path = benchlib_path + "/content/sorttable.js"
    assert os.path.exists(sorttable_js_path), "sorttable.js not at %s" %(sorttable_js_path)
    sorttable_js_html = '<script src="' + sorttable_js_path + '" type="text/javascript"></script>'
    if args.sort_js_mode == 2:
        shutil.copy(sorttable_js_path, out_folder)
        sorttable_js_path = out_folder + "/sorttable.js"
        sorttable_js_html = '<script src="' + sorttable_js_path + '" type="text/javascript"></script>'
    elif args.sort_js_mode == 3:
        js_code = read_file_content_into_str_var(sorttable_js_path)
        sorttable_js_html = "<script>\n" + js_code + "\n</script>\n"

    # Logo path.
    logo_path_html = "logo.png"
    logo_path_out_folder = out_folder + "/" + logo_path_html
    if args.plot_abs_paths:
        logo_path_html = out_folder + "/" + logo_path_html

    logo_path = benchlib_path + "/content/logo.png"
    assert os.path.exists(logo_path), "logo.png not found in %s" %(logo_path)
    shutil.copy(logo_path, logo_path_out_folder)

    # HTML head section.
    html_head = """<!DOCTYPE html>
<html>
<head>
<title>RBPBench - GO Enrichment Analysis Report</title>

<style>
    th, td {
        border: 1px solid black;
        padding: 8px;
        text-align: center;
    }
    th {
        background-color: #f2f2f2;
    }
    .page {
        page-break-after: always;
    }
    .title-container {
        display: flex;
        align-items: center;
    }
    .title-container img {
        margin-right: 10px;
    }
</style>

</head>

<div class="title-container">
    <img src="%s" alt="Logo" width="175">
    <h1>GO Enrichment Analysis Report</h1>
</div>

<body>
""" %(logo_path_html)

    # HTML tail section.
    html_tail = """
%s
</body>
</html>
""" %(sorttable_js_html)

    # Markdown part.
    mdtext = """

List of available statistics and plots generated
by RBPBench (%s, rbpbench goa):

- [GO enrichment analysis results](#goa-results)""" %(version_str)
    mdtext += "\n"
    mdtext += "\n&nbsp;\n"

    mdtext += """
## GO enrichment analysis results ### {#goa-results}

"""
    c_goa_results = 0
    if isinstance(goa_results_df, pd.DataFrame) and not goa_results_df.empty:
        c_goa_results = len(goa_results_df)

    filter_purified_info = " GO terms with significantly higher and lower concentration ([e,p]) in study group are shown."
    filter_purified_info2 = "significant"
    if args.goa_filter_purified:
        filter_purified_info = " Only GO terms with significantly higher concentration in study group are shown."
        filter_purified_info2 = "significantly enriched"
        c_goa_results = len(goa_results_df[goa_results_df["enrichment"] == "e"])
    filter_further_info = ""
    if args.goa_max_child is not None: 
        filter_further_info += " Only GO terms with <= %i children are shown." %(args.goa_max_child)
    if args.goa_min_depth is not None:
        filter_further_info += " Only GO terms with >= %i depth are shown." %(args.goa_min_depth)
    if filter_further_info:
        filter_further_info += " Note that additional filters (children + depth) can result in an empty table. For all significant GO terms (i.e., unfiltered results) check *goa_results.tsv* output table."

    if c_goa_results > 0:

        mdtext += """
**Table:** GO enrichment analysis results. # of %s GO terms found: %i. Filter p-value threshold (on corrected p-value) = %s. # of target genes used for GOA: %i. # of background genes used for GOA: %i.
%s %s

""" %(filter_purified_info2, c_goa_results, str(goa_stats_dic["pval_thr"]), goa_stats_dic["c_target_genes_goa"], goa_stats_dic["c_background_genes_goa"], filter_purified_info, filter_further_info)

        mdtext += '<table style="max-width: 1200px; width: 100%; border-collapse: collapse; line-height: 0.9;">' + "\n"
        mdtext += "<thead>\n"
        mdtext += "<tr>\n"
        mdtext += "<th>GO</th>\n"
        mdtext += "<th>Term</th>\n"
        mdtext += "<th>Class</th>\n"
        mdtext += "<th>p-value</th>\n"
        mdtext += "<th>[e,p]</th>\n"
        mdtext += "<th>Depth</th>\n"
        mdtext += "<th># child</th>\n"
        mdtext += "<th># genes</th>\n"
        mdtext += "<th># study</th>\n"
        mdtext += "<th>% genes</th>\n"
        mdtext += "</tr>\n"
        mdtext += "</thead>\n"
        mdtext += "<tbody>\n"

        for index, row in goa_results_df.iterrows():

            go_id = row['GO']
            go_term = row['term']
            go_class = row['class']
            # go_p = row['p']
            go_p_corr = row['p_corr']
            go_enrichment = row['enrichment']
            go_depth = row['depth']
            go_n_genes = row['n_genes']
            go_n_study = row['n_study']
            go_perc_genes = row['perc_genes']
            go_n_children = row['n_children']

            if args.goa_filter_purified:
                if go_enrichment == "p":
                    continue

            if args.goa_max_child is not None:
                if go_n_children > args.goa_max_child:
                    continue
            if args.goa_min_depth is not None:
                if go_depth < args.goa_min_depth:
                    continue

            mdtext += '<tr>' + "\n"
            mdtext += "<td>" + go_id + "</td>\n"
            mdtext += "<td>" + go_term + "</td>\n"
            mdtext += "<td>" + go_class + "</td>\n"
            mdtext += "<td>" + str(go_p_corr) + "</td>\n"
            mdtext += "<td>" + go_enrichment + "</td>\n"
            mdtext += "<td>" + str(go_depth) + "</td>\n"
            mdtext += "<td>" + str(go_n_children) + "</td>\n"
            mdtext += "<td>" + str(go_n_genes) + "</td>\n"
            mdtext += "<td>" + str(go_n_study) + "</td>\n"
            mdtext += "<td>" + str(go_perc_genes) + "</td>\n"
            mdtext += '</tr>' + "\n"

        mdtext += '</tbody>' + "\n"
        mdtext += '</table>' + "\n"
        
        mdtext += "\n&nbsp;\n&nbsp;\n"
        mdtext += "\nColumn IDs have the following meanings: "
        mdtext += "**GO** -> gene ontology (GO) ID, "
        mdtext += "**Term** -> GO term / name, "
        mdtext += "**Class** -> GO term class (biological_process, molecular_function, or cellular_component), "
        mdtext += "**p-value** -> multiple testing corrected (BH) p-value, "
        mdtext += "**[e,p]** -> e: enriched, i.e., GO term with significantly higher concentration, p: purified, GO term with significantly lower concentration), "
        mdtext += "**Depth** -> depth / level of GO term in GO hierarchy (the higher number, the more specific), "
        mdtext += "**# child** -> number of GO term children, "
        mdtext += "**# genes** -> number of genes associated with GO term, "
        mdtext += "**# study** -> number of genes in study (i.e., target genes), "
        mdtext += "**% genes** -> percentage of study genes associated with GO term." + "\n"
        mdtext += "\n&nbsp;\n"

    else:

        if "c_target_genes_goa" in goa_stats_dic:

            mdtext += """

No %s GO terms found given p-value threshold of %s. # of target genes used for GOA: %i. # of background genes used for GOA: %i.

&nbsp;

""" %(filter_purified_info2, str(goa_stats_dic["pval_thr"]), goa_stats_dic["c_target_genes_goa"], goa_stats_dic["c_background_genes_goa"])

        else:

            mdtext += """

No significant GO terms found due to no GO IDs associated with target genes. # of initial target genes (i.e., genes overlapping with --in regions): %i.

&nbsp;

""" %(goa_stats_dic["c_target_genes_pre_filter"])

    # Convert mdtext to html.
    md2html = markdown(mdtext, extensions=['attr_list', 'tables'])

    # OUTMD = open(md_out,"w")
    # OUTMD.write("%s\n" %(mdtext))
    # OUTMD.close()

    html_content = html_head + md2html + html_tail

    OUTHTML = open(html_out,"w")
    OUTHTML.write("%s\n" %(html_content))
    OUTHTML.close()


################################################################################

def guess_chr_id_style(chr_ids_dic):
    """
    Guess chromosome style present in chr_ids_dic.

    chr_style:
        1: chr1, chr2, ..., chrX, chrM
        2: 1, 2, ... , X, MT

    >>> chr_ids_dic = {"chr1": 1, "chr2": 1, "chrX": 1, "chrM": 1, "chr11_KI270721v1_random": 1}
    >>> guess_chr_id_style(chr_ids_dic)
    1
    >>> chr_ids_dic = {"1": 1, "2": 1, "X": 1, "MT": 1}
    >>> guess_chr_id_style(chr_ids_dic)
    2
    
    """
    assert chr_ids_dic, "chr_ids_dic empty (no chromosome IDs found in --genome FASTA file?)"

    chr_style = 1
    
    for chr_id in chr_ids_dic:
        if chr_id.startswith("chr"):
            chr_style = 1
            break
        else:
            chr_style = 2
            break

    if chr_style == 1:
        for chr_id in chr_ids_dic:
            assert chr_id.startswith("chr"), "inconsistent chromosome IDs in --genome FASTA file (chr prefix expected but not present for ID %s)" %(chr_id)
    else:
        for chr_id in chr_ids_dic:
            assert not chr_id.startswith("chr"), "inconsistent chromosome IDs in --genome FASTA file (both chr prefix and no chr prefix present)"

    return chr_style


################################################################################

def batch_generate_html_report(args,
                               dataset_ids_list, 
                               kmer_list,
                               kmer_freqs_ll,
                               id2infos_dic, 
                               id2reg_annot_dic,
                               id2hit_reg_annot_dic,
                               benchlib_path,
                               seq_len_stats_ll,
                               seq_feat_ll=False,
                               seq_var_ll=False,
                               seq_var_kmer_l=False,
                               html_report_out="report.rbpbench_batch.html",
                               id2hk_gene_stats_dic=False,
                               id2motif_enrich_stats_dic=False,
                               id2regex_stats_dic=False,
                               regex_annot_dic=False,
                               id2occ_list_dic=False,
                               gene_occ_cooc_plot=False,
                               ei_ol_stats_dic=False,
                               plotly_full_html=False,
                               plotly_embed_style=1,
                               add_motif_db_info=False,
                               heatmap_cluster_olo=False,
                               goa_results_df=False,
                               goa_stats_dic=False,
                               dsid2add_annot_stats_dic=False,
                               plots_subfolder="html_report_plots"):
    """
    Create plots for RBPBench batch run results.

    """

    # Minimum dataset_ids_list and kmer_freqs_ll needed for plotting anything.
    assert dataset_ids_list, "no dataset IDs found for report creation"
    assert kmer_freqs_ll, "no k-mer frequencies found for report creation"

    # Use absolute paths?
    out_folder = args.out_folder
    if args.plot_abs_paths:
        out_folder = os.path.abspath(out_folder)

    # Version string.
    version_str = "v" + args.version

    plots_folder = plots_subfolder
    plots_out_folder = out_folder + "/" + plots_folder
    if args.plot_abs_paths:
        plots_folder = plots_out_folder

    # Delete folder if already present.
    if os.path.exists(plots_out_folder):
        shutil.rmtree(plots_out_folder)
    os.makedirs(plots_out_folder)

    html_out = out_folder + "/" + "report.rbpbench_batch.html"
    md_out = out_folder + "/" + "report.rbpbench_batch.md"
    if html_report_out:
        html_out = html_report_out

    """
    Setup sorttable.js to make tables in HTML sortable.

    """
    sorttable_js_path = benchlib_path + "/content/sorttable.js"
    assert os.path.exists(sorttable_js_path), "sorttable.js not at %s" %(sorttable_js_path)
    sorttable_js_html = '<script src="' + sorttable_js_path + '" type="text/javascript"></script>'
    if args.sort_js_mode == 2:
        shutil.copy(sorttable_js_path, plots_out_folder)
        sorttable_js_path = plots_folder + "/sorttable.js"
        sorttable_js_html = '<script src="' + sorttable_js_path + '" type="text/javascript"></script>'
    elif args.sort_js_mode == 3:
        js_code = read_file_content_into_str_var(sorttable_js_path)
        sorttable_js_html = "<script>\n" + js_code + "\n</script>\n"

    """
    Setup plotly .js to support plotly plots.

    https://plotly.com/javascript/getting-started/#download
    plotly-latest.min.js
    Packaged version: plotly-2.20.0.min.js
    
    """

    include_plotlyjs = "cdn"
    # plotly_full_html = False
    plotly_js_html = ""
    plotly_js_path = benchlib_path + "/content/plotly-2.20.0.min.js"
    assert os.path.exists(plotly_js_path), "plotly .js %s not found" %(plotly_js_path)
    if args.plotly_js_mode == 2:
        include_plotlyjs = plotly_js_path
    elif args.plotly_js_mode == 3:
        shutil.copy(plotly_js_path, plots_out_folder)
        include_plotlyjs = "plotly-2.20.0.min.js" # Or plots_folder + "/plotly-2.20.0.min.js" ?
    elif args.plotly_js_mode == 4:
        include_plotlyjs = True
        # plotly_full_html = False # Don't really need full html (head body ..) in plotly html.
    elif args.plotly_js_mode == 5:
        plotly_js_web = "https://cdn.plot.ly/plotly-2.25.2.min.js"
        plotly_js_html = '<script src="' + plotly_js_web + '"></script>' + "\n"
        include_plotlyjs = False
        # plotly_full_html = True
    elif args.plotly_js_mode == 6:
        shutil.copy(plotly_js_path, plots_out_folder)
        plotly_js = plots_folder + "/plotly-2.20.0.min.js"
        plotly_js_html = '<script src="' + plotly_js + '"></script>' + "\n"
        include_plotlyjs = False
    elif args.plotly_js_mode == 7:
        js_code = read_file_content_into_str_var(plotly_js_path)
        plotly_js_html = "<script>\n" + js_code + "\n</script>\n"
        include_plotlyjs = False
        # plotly_full_html = True

    annot_dic = {}
    annot2color_dic = {}
    if id2reg_annot_dic:  # if --gtf provided.

        # All annotations.
        annot_dic = {}
        c_annot = 0
        for internal_id in id2reg_annot_dic:
            for annot in id2reg_annot_dic[internal_id]:
                c_annot += 1
                if annot not in annot_dic:
                    annot_dic[annot] = 1
                else:
                    annot_dic[annot] += 1

        # Annotation color assignments.
        hex_colors = get_hex_colors_list(min_len=len(annot_dic))

        idx = 0
        for annot in sorted(annot_dic, reverse=False):
            hc = hex_colors[idx]
            # print("Assigning hex color %s to annotation %s ... " %(hc, annot))
            annot2color_dic[annot] = hex_colors[idx]
            idx += 1

    # Logo path.
    logo_path_html = plots_folder + "/logo.png"
    if not os.path.exists(logo_path_html):
        logo_path = benchlib_path + "/content/logo.png"
        assert os.path.exists(logo_path), "logo.png not found in %s" %(logo_path)
        shutil.copy(logo_path, plots_out_folder)

    # HTML head section.
    html_head = """<!DOCTYPE html>
<html>
<head>
<title>RBPBench - Batch Report</title>
%s

<style>
    th, td {
        border: 1px solid black;
        padding: 8px;
        text-align: center;
    }
    th {
        background-color: #f2f2f2;
    }
    .page {
        page-break-after: always;
    }
    .title-container {
        display: flex;
        align-items: center;
    }
    .title-container img {
        margin-right: 10px; /* Adjust the spacing as needed */
    }
</style>

</head>

<div class="title-container">
    <img src="%s" alt="Logo" width="175">
    <h1>Batch Report</h1>
</div>

<body>
""" %(plotly_js_html, logo_path_html)

    # HTML tail section.
    html_tail = """
%s
</body>
</html>
""" %(sorttable_js_html)

    # Markdown part.
    mdtext = """

List of available statistics and plots generated
by RBPBench (%s, rbpbench batch):

- [Input datasets sequence length statistics](#seq-len-stats)""" %(version_str)
    
    mdtext += "\n"

    if ei_ol_stats_dic:
        mdtext += "- [Input datasets exon-intron overlap statistics](#ei-ol-stats)\n"
        mdtext += "- [Input datasets exon-intron overlap comparative plot](#ei-ol-plot)\n"

    if id2motif_enrich_stats_dic: # if not empty.
        mdtext += "- [Input datasets RBP region score motif enrichment statistics](#motif-enrich-stats)\n"

    if id2regex_stats_dic:  # if not empty.
        mdtext += "- [Regular expression motif enrichment statistics](#regex-enrich-stats)\n"
        mdtext += "- [Regular expression RBP motif co-occurrence statistics](#regex-rbp-cooc-stats)\n"

    if seq_feat_ll:
        mdtext += "- [Input datasets nucleotide percentages statistics](#nt-perc-stats)\n"

    mdtext += "- [Input datasets k-mer frequencies comparative plot](#kmer-comp-plot)\n"

    if seq_var_ll:
        mdtext += "- [Input datasets k-mer variation comparative plot](#seq-var-plot)\n"

    if id2reg_annot_dic:  # if --gtf provided.
        mdtext += "- [Input datasets occupied gene regions comparative plot](#occ-comp-plot)\n"
        if gene_occ_cooc_plot:
            mdtext += "- [Input datasets occupied gene regions similarity heat map](#cooc-heat-map)\n"
        mdtext += "- [Input datasets genomic region annotations comparative plot](#annot-comp-plot)\n"

    if id2reg_annot_dic:  # if --gtf provided.
        data_idx = 0
        internal_ids_list = []
        for internal_id in id2infos_dic:
            rbp_id = id2infos_dic[internal_id][0]
            data_id = id2infos_dic[internal_id][1]
            method_id = id2infos_dic[internal_id][2]
            database_id = id2infos_dic[internal_id][3]
            if add_motif_db_info:
                combined_id = rbp_id + "," + database_id + "," + method_id + "," + data_id
            else:
                combined_id = rbp_id + "," + method_id + "," + data_id
            internal_ids_list.append(internal_id)
            mdtext += "- [%s region annotations](#annot-plot-%i)\n" %(combined_id, data_idx)
            data_idx += 1
        # mdtext += "\n&nbsp;\n"

    if dsid2add_annot_stats_dic:
        mdtext += "- [Additional region annotation statistics](#add-annot-stats)\n"

    if args.run_goa:
        mdtext += "- [GO enrichment analysis results](#goa-results)\n"

    add_head_info = ""
    if args.bed_sc_thr is not None:
        add_head_info = " BED score threshold (--bed-sc-thr) = %s" %(str(args.bed_sc_thr))
        if args.bed_sc_thr_rev_filter:
            add_head_info += " (reverse filtering applied, i.e., the lower the better)."
        else:
            add_head_info += "."

    mdtext += "\nUsed motif database = %s. FIMO p-value threshold (--fimo-pval) = %s.%s Region extension (upstream, downstream) = (%i, %i).\n" %(args.motif_db_str, str(args.fimo_pval), add_head_info, args.ext_up, args.ext_down)
    mdtext += "\n&nbsp;\n"


    """
    Input sequence stats table.

    """

    dataset_id_format = "rbp_id,method_id,data_id"
    if add_motif_db_info:
        dataset_id_format = "rbp_id,motif_database_id,method_id,data_id"

    seq_stats_info = ""
    if args.unstranded and not args.unstranded_ct:
        seq_stats_info = "--unstranded option selected, i.e., both strands of a region are included in the length statistics, but the two strands are counted as one region."
    elif args.unstranded and args.unstranded_ct:
        seq_stats_info = "--unstranded and --unstranded-ct options selected, i.e., both strands of a region are included in the length statistics, and each strand counts as separate region."

    mdtext += """
## Input datasets sequence length statistics ### {#seq-len-stats}

**Table:** Sequence length statistics of input datasets.
Sequence length in nt.
Input dataset ID format: %s. %s

""" %(dataset_id_format, seq_stats_info)

    # return [nr_seqs, seq_len_mean, seq_len_median, seq_len_q1, seq_len_q3, seq_len_min, seq_len_max]

    # mdtext += '| Dataset ID | # input regions | mean length | median length | Q1 percentile | Q3 percentile | min length | max length |' + " \n"
    # mdtext += '| :-: | :-: | :-: | :-: | :-: | :-: | :-: | :-: |' + " \n"

    mdtext += '<table style="max-width: 1200px; width: 100%; border-collapse: collapse; line-height: 0.8;">' + "\n"
    mdtext += "<thead>\n"
    mdtext += "<tr>\n"
    mdtext += "<th>Dataset ID</th>\n"
    mdtext += "<th># input regions</th>\n"
    mdtext += "<th>mean length</th>\n"
    mdtext += "<th>median length</th>\n"
    mdtext += "<th>Q1 percentile</th>\n"
    mdtext += "<th>Q3 percentile</th>\n"
    mdtext += "<th>min length</th>\n"
    mdtext += "<th>max length</th>\n"
    mdtext += "</tr>\n"
    mdtext += "</thead>\n"
    mdtext += "<tbody>\n"

    for seq_len_stats in seq_len_stats_ll:
        # mdtext += "| %s | %i | %.1f | %.1f | %.1f | %.1f | %i | %i |\n" %(seq_len_stats[0], seq_len_stats[1], seq_len_stats[2], seq_len_stats[3], seq_len_stats[4], seq_len_stats[5], seq_len_stats[6], seq_len_stats[7])
        mdtext += "<tr>\n"
        mdtext += "<td>%s</td>\n" %(seq_len_stats[0])
        mdtext += "<td>%i</td>\n" %(seq_len_stats[1])
        mdtext += "<td>%.1f</td>\n" %(seq_len_stats[2])
        mdtext += "<td>%.1f</td>\n" %(seq_len_stats[3])
        mdtext += "<td>%.1f</td>\n" %(seq_len_stats[4])
        mdtext += "<td>%.1f</td>\n" %(seq_len_stats[5])
        mdtext += "<td>%i</td>\n" %(seq_len_stats[6])
        mdtext += "<td>%i</td>\n" %(seq_len_stats[7])
        mdtext += "</tr>\n"
    
    mdtext += "</tbody>\n"
    mdtext += "</table>\n"

    # mdtext += "\n&nbsp;\n&nbsp;\n"
    mdtext += "\n&nbsp;\n"

    #     mdtext += "| :-: | :-: | :-: | :-: | :-: |\n"
    # for rbp_id, wc_pval in sorted(pval_dic.items(), key=lambda item: item[1], reverse=False):
    #     wc_pval_str = convert_sci_not_to_decimal(wc_pval)  # Convert scientific notation to decimal string for sorting to work.
    #     c_hit_reg = search_rbps_dic[rbp_id].c_hit_reg
    #     perc_hit_reg = search_rbps_dic[rbp_id].perc_hit_reg
    #     c_uniq_motif_hits = search_rbps_dic[rbp_id].c_uniq_motif_hits
    #     mdtext += "| %s | %i | %.2f | %i | %s |\n" %(rbp_id, c_hit_reg, perc_hit_reg, c_uniq_motif_hits, wc_pval)
    # mdtext += "\n&nbsp;\n&nbsp;\n"
    # mdtext += "\nColumn IDs have the following meanings: "
    # mdtext += "**RBP ID** -> RBP ID from database or user-defined (typically RBP name), "
    # mdtext += '**# hit regions** -> number of input genomic regions with motif hits (after filtering and optional extension), '
    # mdtext += '**% hit regions** -> percentage of hit regions over all regions (i.e. how many input regions contain >= 1 RBP binding motif), '
    # mdtext += '**# motif hits** -> number of unique motif hits in input regions (removed double counts), '
    # mdtext += '**p-value** -> Wilcoxon rank-sum test p-value.' + "\n"
    # mdtext += "\n&nbsp;\n"

    id2eib_stats_dic = {}  # For plotly plot.
    id2eib_perc_dic = {}  # For plotly plot.

    if ei_ol_stats_dic:  # True if --gtf provided.

        ib_len =  args.gtf_intron_border_len
        eib_len = args.ei_border_len

        mdtext += """
## Input datasets exon-intron overlap statistics ### {#ei-ol-stats}

**Table:** Exon, intron + border region overlap statistics for each input dataset.
Considered intron border region length = %i nt (change via --gtf-intron-border-len). 
Considered exon-intron border region = +/- %i nt relative to border.
%s

""" %(ib_len, eib_len, seq_stats_info)


        mdtext += '<table style="max-width: 1200px; width: 100%; border-collapse: collapse; line-height: 0.8;">' + "\n"
        mdtext += "<thead>\n"
        mdtext += "<tr>\n"
        mdtext += "<th>Dataset ID</th>\n"
        mdtext += "<th># input regions</th>\n"
        mdtext += "<th>% exon regions</th>\n"
        mdtext += "<th>% intron regions</th>\n"
        mdtext += "<th>% us intron border</th>\n"
        mdtext += "<th>% ds intron border</th>\n"
        mdtext += "<th>% distant us intron</th>\n"
        mdtext += "<th>% distant ds intron</th>\n"
        mdtext += "<th>% exon-intron border</th>\n"
        mdtext += "<th>% first exon</th>\n"
        mdtext += "<th>% last exon</th>\n"
        mdtext += "<th>% single exon</th>\n"
        mdtext += "</tr>\n"
        mdtext += "</thead>\n"
        mdtext += "<tbody>\n"

        # OUTEIB = open(eib_ol_table_out,"w")
        # OUTEIB.write("dataset_id\tperc_exons\tperc_introns\tperc_us_ib\tperc_ds_ib\tperc_eib\tc_regions\tc_exon_sites\tc_intron_sites\tc_us_ib_sites\tc_ds_ib_sites\tc_eib_sites\tc_tr_ids\tc_tr_ids_with_sites\n")

        for internal_id in ei_ol_stats_dic:
            combined_id = internal_id
            rbp_id = id2infos_dic[internal_id][0]
            data_id = id2infos_dic[internal_id][1]
            method_id = id2infos_dic[internal_id][2]
            database_id = id2infos_dic[internal_id][3]
            if add_motif_db_info:
                combined_id = rbp_id + "," + database_id + "," + method_id + "," + data_id
            else:
                combined_id = rbp_id + "," + method_id + "," + data_id

            # ei_ol_stats = ei_ol_stats_dic[internal_id]
            # exon_sites_perc = 0.0
            # exon_sites_ratio = 0.0
            # if ei_ol_stats.c_exon_sites and ei_ol_stats.c_input_sites:
            #     exon_sites_perc = round(ei_ol_stats.c_exon_sites / ei_ol_stats.c_input_sites * 100, 1)
            #     exon_sites_ratio = round(ei_ol_stats.c_exon_sites / ei_ol_stats.c_input_sites, 3)
            # intron_sites_perc = 0.0
            # intron_sites_ratio = 0.0
            # if ei_ol_stats.c_intron_sites and ei_ol_stats.c_input_sites:
            #     intron_sites_perc = round(ei_ol_stats.c_intron_sites / ei_ol_stats.c_input_sites * 100, 1)
            #     intron_sites_ratio = round(ei_ol_stats.c_intron_sites / ei_ol_stats.c_input_sites, 3)
            # us_ib_sites_perc = 0.0
            # us_ib_sites_ratio = 0.0
            # if ei_ol_stats.c_us_ib_sites and ei_ol_stats.c_input_sites:
            #     us_ib_sites_perc = round(ei_ol_stats.c_us_ib_sites / ei_ol_stats.c_input_sites * 100, 1)
            #     us_ib_sites_ratio = round(ei_ol_stats.c_us_ib_sites / ei_ol_stats.c_input_sites, 3)
            # ds_ib_sites_perc = 0.0
            # ds_ib_sites_ratio = 0.0
            # if ei_ol_stats.c_ds_ib_sites and ei_ol_stats.c_input_sites:
            #     ds_ib_sites_perc = round(ei_ol_stats.c_ds_ib_sites / ei_ol_stats.c_input_sites * 100, 1)
            #     ds_ib_sites_ratio = round(ei_ol_stats.c_ds_ib_sites / ei_ol_stats.c_input_sites, 3)
            # eib_sites_perc = 0.0
            # eib_sites_ratio = 0.0
            # if ei_ol_stats.c_eib_sites and ei_ol_stats.c_input_sites:
            #     eib_sites_perc = round(ei_ol_stats.c_eib_sites / ei_ol_stats.c_input_sites * 100, 1)
            #     eib_sites_ratio = round(ei_ol_stats.c_eib_sites / ei_ol_stats.c_input_sites, 3)

            ei_ol_stats = ei_ol_stats_dic[internal_id]
            exon_sites_perc = 0.0
            exon_sites_ratio = 0.0
            if ei_ol_stats.c_exon_sites and ei_ol_stats.c_input_sites:
                exon_sites_perc = round(ei_ol_stats.c_exon_sites / ei_ol_stats.c_input_sites * 100, 1)
                exon_sites_ratio = round(ei_ol_stats.c_exon_sites / ei_ol_stats.c_input_sites, 3)
            intron_sites_perc = 0.0
            intron_sites_ratio = 0.0
            if ei_ol_stats.c_intron_sites and ei_ol_stats.c_input_sites:
                intron_sites_perc = round(ei_ol_stats.c_intron_sites / ei_ol_stats.c_input_sites * 100, 1)
                intron_sites_ratio = round(ei_ol_stats.c_intron_sites / ei_ol_stats.c_input_sites, 3)
            us_ib_sites_perc = 0.0
            us_ib_sites_ratio = 0.0
            if ei_ol_stats.c_us_ib_sites and ei_ol_stats.c_input_sites:
                us_ib_sites_perc = round(ei_ol_stats.c_us_ib_sites / ei_ol_stats.c_input_sites * 100, 1)
                us_ib_sites_ratio = round(ei_ol_stats.c_us_ib_sites / ei_ol_stats.c_input_sites, 3)
            ds_ib_sites_perc = 0.0
            ds_ib_sites_ratio = 0.0
            if ei_ol_stats.c_ds_ib_sites and ei_ol_stats.c_input_sites:
                ds_ib_sites_perc = round(ei_ol_stats.c_ds_ib_sites / ei_ol_stats.c_input_sites * 100, 1)
                ds_ib_sites_ratio = round(ei_ol_stats.c_ds_ib_sites / ei_ol_stats.c_input_sites, 3)
            us_ib_dist_sites_perc = 0.0
            us_ib_dist_sites_ratio = 0.0
            if ei_ol_stats.c_us_ib_dist_sites and ei_ol_stats.c_input_sites:
                us_ib_dist_sites_perc = round(ei_ol_stats.c_us_ib_dist_sites / ei_ol_stats.c_input_sites * 100, 1)
                us_ib_dist_sites_ratio = round(ei_ol_stats.c_us_ib_dist_sites / ei_ol_stats.c_input_sites, 3)
            ds_ib_dist_sites_perc = 0.0
            ds_ib_dist_sites_ratio = 0.0
            if ei_ol_stats.c_ds_ib_dist_sites and ei_ol_stats.c_input_sites:
                ds_ib_dist_sites_perc = round(ei_ol_stats.c_ds_ib_dist_sites / ei_ol_stats.c_input_sites * 100, 1)
                ds_ib_dist_sites_ratio = round(ei_ol_stats.c_ds_ib_dist_sites / ei_ol_stats.c_input_sites, 3)
            first_exon_sites_perc = 0.0
            first_exon_sites_ratio = 0.0
            if ei_ol_stats.c_first_exon_sites and ei_ol_stats.c_input_sites:
                first_exon_sites_perc = round(ei_ol_stats.c_first_exon_sites / ei_ol_stats.c_input_sites * 100, 1)
                first_exon_sites_ratio = round(ei_ol_stats.c_first_exon_sites / ei_ol_stats.c_input_sites, 3)
            last_exon_sites_perc = 0.0
            last_exon_sites_ratio = 0.0
            if ei_ol_stats.c_last_exon_sites and ei_ol_stats.c_input_sites:
                last_exon_sites_perc = round(ei_ol_stats.c_last_exon_sites / ei_ol_stats.c_input_sites * 100, 1)
                last_exon_sites_ratio = round(ei_ol_stats.c_last_exon_sites / ei_ol_stats.c_input_sites, 3)
            single_exon_sites_perc = 0.0
            single_exon_sites_ratio = 0.0
            if ei_ol_stats.c_single_exon_sites and ei_ol_stats.c_input_sites:
                single_exon_sites_perc = round(ei_ol_stats.c_single_exon_sites / ei_ol_stats.c_input_sites * 100, 1)
                single_exon_sites_ratio = round(ei_ol_stats.c_single_exon_sites / ei_ol_stats.c_input_sites, 3)
            eib_sites_perc = 0.0
            eib_sites_ratio = 0.0
            if ei_ol_stats.c_eib_sites and ei_ol_stats.c_input_sites:
                eib_sites_perc = round(ei_ol_stats.c_eib_sites / ei_ol_stats.c_input_sites * 100, 1)
                eib_sites_ratio = round(ei_ol_stats.c_eib_sites / ei_ol_stats.c_input_sites, 3)

            id2eib_stats_dic[internal_id] = [combined_id, ei_ol_stats.c_input_sites, exon_sites_ratio, intron_sites_ratio, us_ib_sites_ratio, ds_ib_sites_ratio, us_ib_dist_sites_ratio, ds_ib_dist_sites_ratio, eib_sites_ratio, first_exon_sites_ratio, last_exon_sites_ratio, single_exon_sites_ratio]
            id2eib_perc_dic[internal_id] = [exon_sites_perc, intron_sites_perc, us_ib_sites_perc, ds_ib_sites_perc, us_ib_dist_sites_perc, ds_ib_dist_sites_perc, eib_sites_perc, first_exon_sites_perc, last_exon_sites_perc, single_exon_sites_perc]

            mdtext += "<tr>\n"
            mdtext += "<td>%s</td>\n" %(combined_id)
            mdtext += "<td>%i</td>\n" %(ei_ol_stats.c_input_sites)
            mdtext += "<td>%.1f</td>\n" %(exon_sites_perc)
            mdtext += "<td>%.1f</td>\n" %(intron_sites_perc)
            mdtext += "<td>%.1f</td>\n" %(us_ib_sites_perc)
            mdtext += "<td>%.1f</td>\n" %(ds_ib_sites_perc)
            mdtext += "<td>%.1f</td>\n" %(us_ib_dist_sites_perc)
            mdtext += "<td>%.1f</td>\n" %(ds_ib_dist_sites_perc)
            mdtext += "<td>%.1f</td>\n" %(eib_sites_perc)
            mdtext += "<td>%.1f</td>\n" %(first_exon_sites_perc)
            mdtext += "<td>%.1f</td>\n" %(last_exon_sites_perc)
            mdtext += "<td>%.1f</td>\n" %(single_exon_sites_perc)
            mdtext += "</tr>\n"

        mdtext += "</tbody>\n"
        mdtext += "</table>\n"

        mdtext += "\n&nbsp;\n&nbsp;\n"
        mdtext += "\nColumn IDs have the following meanings: "
        mdtext += "**Dataset ID** -> Dataset ID for input dataset with following format: %s, " %(dataset_id_format)
        mdtext += '**# input regions** -> number of considered input regions from input dataset, '
        mdtext += '**% exon regions** -> % of input regions overlapping with exon regions (== exonic regions), '
        mdtext += '**% intron regions** -> % of input regions overlapping with intron regions (== intronic regions), '
        mdtext += '**%% us intron border** -> %% of intronic regions closer to intron upstream borders + within first %i nt of intron, ' %(ib_len)
        mdtext += '**%% ds intron border** -> %% of intronic regions closer to intron downstream borders + within first %i nt of intron, ' %(ib_len)
        mdtext += '**%% distant us intron** -> %% of intronic regions closer to intron upstream borders and > %i nt away from intron borders, ' %(ib_len)
        mdtext += '**%% distant ds intron** -> %% of intronic regions closer to intron downstream borders and > %i nt away from intron borders, ' %(ib_len)
        mdtext += '**%% exon-intron border** -> %% of input regions overlapping with exon-intron borders (+/- %i nt of exon-intron borders). ' %(eib_len)
        mdtext += "**% first exon** -> % of input regions overlapping with transcript first exons, "
        mdtext += '**% last exon** -> % of input regions overlapping with transcript last exons, '
        mdtext += '**% single exon** -> % of input regions overlapping with single exon transcripts. '
        mdtext += "**NOTE** that for upstream/downstream intron end overlaps, the input region is always assigned "
        mdtext += "to the closest intron border (based on distance between input region center position and intron ends), "
        mdtext += "independent of intron length (i.e., intron length can also be < %i nt)." %(ib_len)
        mdtext += "\n&nbsp;\n"

        # id2eib_stats_dic[internal_id] = [combined_id, ei_ol_stats.c_input_sites, exon_sites_ratio, intron_sites_ratio, us_ib_sites_ratio, ds_ib_sites_ratio, eib_sites_ratio]

        mdtext += """
## Input datasets exon-intron overlap comparative plot ### {#ei-ol-plot}

"""
        if len(dataset_ids_list) > 3:

            eib_comp_plot_plotly =  "eib_comparative_plot.plotly.html"
            eib_comp_plot_plotly_out = plots_out_folder + "/" + eib_comp_plot_plotly

            create_eib_comp_plot_plotly(id2eib_stats_dic, id2eib_perc_dic, 
                                        eib_comp_plot_plotly_out,
                                        plot_3d=False,
                                        id2hk_gene_stats_dic=False,
                                        include_plotlyjs=include_plotlyjs,
                                        full_html=plotly_full_html)

            plot_path = plots_folder + "/" + eib_comp_plot_plotly

            if args.plotly_js_mode in [5, 6, 7]:
                # Read in plotly code.
                # mdtext += '<div style="width: 1200px; height: 1200px; align-items: center;">' + "\n"
                js_code = read_file_content_into_str_var(eib_comp_plot_plotly_out)
                js_code = js_code.replace("height:100%; width:100%;", "height:1000px; width:1200px;")
                mdtext += js_code + "\n"
                # mdtext += '</div>'
            else:
                if plotly_embed_style == 1:
                    # mdtext += '<div class="container-fluid" style="margin-top:40px">' + "\n"
                    mdtext += "<div>\n"
                    mdtext += '<iframe src="' + plot_path + '" width="1200" height="1000"></iframe>' + "\n"
                    mdtext += '</div>'
                elif plotly_embed_style == 2:
                    mdtext += '<object data="' + plot_path + '" width="1200" height="1000"> </object>' + "\n"

            mdtext += r"""

**Figure:** Exon, intron + border region overlap statistics visualized as 2D PCA plot.
The closer two datasets (i.e., the dots representing them), the more similar the datasets are w.r.t. their exon-intron overlap statistics.
Hover box shows: 
dataset ID (bold-faced, format: %s), 
**Input regions** -> \# of dataset regions,
**Exon** -> %% of input regions overlapping with exon regions (== exonic regions),
**Intron** -> %% of input regions overlapping with intron regions (== intronic regions),
**US intron border** -> %% of intronic regions closer to intron upstream borders + within first %i nt of intron,
**DS intron border** -> %% of intronic regions closer to intron downstream borders + within first %i nt of intron,
**US distant intron** -> %% of intronic regions closer to intron upstream borders and > %i nt away from intron borders,
**DS distant intron** -> %% of intronic regions closer to intron downstream borders and > %i nt away from intron borders,
**Exon-intron border** -> %% of input regions overlapping with exon-intron borders (+/- %i nt of exon-intron borders),
**First exon** -> %% of input regions overlapping with transcript first exons,
**Last exon** -> %% of input regions overlapping with transcript last exons,
**Single exon** -> %% of input regions overlapping with single exon transcripts.
**NOTE** that for upstream/downstream intron end overlaps, the input region is always assigned
to the closest intron border (based on distance between input region center position and intron ends),
independent of intron length (i.e., intron length can also be < %i nt).
&nbsp;

""" %(dataset_id_format, ib_len, ib_len, ib_len, ib_len, eib_len, ib_len)

        else:
            mdtext += """

No plot generated since < 4 datasets were provided.

&nbsp;

"""


    """

    Input datasets RBP region score motif enrichment statistics.

    Format:
    id2motif_enrich_stats_dic[internal_id] = [c_reg_with_hits, perc_reg_with_hits, c_uniq_motif_hits, wc_pval]


    """

    # Inform about set alterntive hypothesis for Wilcoxon rank sum test.
    wrs_mode_info1 = "Wilcoxon rank-sum test alternative hypothesis is set to 'greater', i.e., low p-values mean hit-containing regions have significantly higher scores."
    wrs_mode_info2 = "higher"
    if args.wrs_mode == 2:
        wrs_mode_info1 = "Wilcoxon rank-sum test alternative hypothesis is set to 'less', i.e., low p-values mean hit-containing regions have significantly lower scores."
        wrs_mode_info2 = "lower"

    if id2motif_enrich_stats_dic:

        mdtext += """
## Input datasets region score motif enrichment statistics ### {#motif-enrich-stats}

**Table:** Input datasets region score motif enrichment statistics for all input datasets.
For each input dataset, consisting of a set of genomic regions with associated scores (set BED score column via --bed-score-col),
RBPbench checks whether regions with RBP motif hits have significantly different scores compared to regions without hits.
%s
In other words, a low test p-value for a given dataset indicates 
that %s-scoring regions are more likely to contain RBP motif hits.
**NOTE** that if scores associated with input genomic regions are all the same, p-values become meaningless 
(i.e., they result in p-values of 1.0).
Likewise, the p-value becomes non-informative if most or all input regions have RBP motif hits (i.e., very high hit region percentages).
By default, BED genomic regions input file column 5 is used as the score column (change with --bed-score-col).

""" %(wrs_mode_info1, wrs_mode_info2)

        mdtext += '<table style="max-width: 1000px; width: 100%; border-collapse: collapse; line-height: 0.8;">' + "\n"
        mdtext += "<thead>\n"
        mdtext += "<tr>\n"
        mdtext += "<th>Dataset ID</th>\n"
        mdtext += "<th># hit regions</th>\n"
        mdtext += "<th>% hit regions</th>\n"
        mdtext += "<th># motif hits</th>\n"
        mdtext += "<th>p-value</th>\n"
        mdtext += "</tr>\n"
        mdtext += "</thead>\n"
        mdtext += "<tbody>\n"

        pval_dic = {}
        for internal_id in id2motif_enrich_stats_dic:
            wc_pval = id2motif_enrich_stats_dic[internal_id][3]
            pval_dic[internal_id] = wc_pval

        for internal_id, wc_pval in sorted(pval_dic.items(), key=lambda item: item[1], reverse=False):

            rbp_id = id2infos_dic[internal_id][0]
            data_id = id2infos_dic[internal_id][1]
            method_id = id2infos_dic[internal_id][2]
            database_id = id2infos_dic[internal_id][3]
            combined_id = rbp_id + "," + method_id + "," + data_id

            c_hit_reg = id2motif_enrich_stats_dic[internal_id][0]
            perc_hit_reg = id2motif_enrich_stats_dic[internal_id][1]
            c_uniq_motif_hits = id2motif_enrich_stats_dic[internal_id][2]

            mdtext += '<tr>' + "\n"
            mdtext += "<td>" + combined_id + "</td>\n"
            mdtext += "<td>" + str(c_hit_reg) + "</td>\n"
            mdtext += "<td>%.2f" %(perc_hit_reg) + "</td>\n"
            mdtext += "<td>" + str(c_uniq_motif_hits) + "</td>\n"
            mdtext += "<td>" + str(wc_pval) + "</td>\n"
            mdtext += '</tr>' + "\n"

        mdtext += "</tbody>\n"
        mdtext += "</table>\n"
        
        mdtext += "\n&nbsp;\n&nbsp;\n"
        mdtext += "\nColumn IDs have the following meanings: "
        mdtext += "**Dataset ID** -> Dataset ID with following format: %s, " %(dataset_id_format)
        mdtext += '**# regions** -> number of input genomic regions in dataset (after filtering and optional extension), '
        mdtext += '**# hit regions** -> number of input genomic regions with motif hits (after filtering and optional extension), '
        mdtext += '**% hit regions** -> percentage of motif hit regions over all regions (i.e., how many input regions contain >= 1 RBP motif hit), '
        mdtext += '**# motif hits** -> number of unique motif hits in input regions (removed double counts), '
        mdtext += '**p-value** -> Wilcoxon rank-sum test p-value.' + "\n"
        mdtext += "\n&nbsp;\n"


    """
    regex motif enrichment statistics.
    
    """

    if id2regex_stats_dic:


        mdtext += """
## Regular expression region score motif enrichment statistics ### {#regex-enrich-stats}

**Table:** Regular expression (regex) '%s' region score motif enrichment statistics for all input datasets.
For each input dataset, consisting of a set of genomic regions with associated scores (set BED score column via --bed-score-col),
RBPbench checks whether regions containing regex hits have significantly different scores compared to regions without regex hits.
%s
In other words, a low test p-value for a given dataset indicates 
that %s-scoring regions are more likely to contain regex hits.
**NOTE** that if scores associated to input genomic regions are all the same, p-values become meaningless 
(i.e., they result in p-values of 1.0).
By default, BED genomic regions input file column 5 is used as the score column (change with --bed-score-col).

""" %(args.regex, wrs_mode_info1, wrs_mode_info2)

        # mdtext += '| Dataset ID  | # regions | # hit regions | % hit regions | # regex hits | p-value |' + " \n"
        # mdtext += '| :-: | :-: | :-: | :-: | :-: | :-: |' + " \n"

        mdtext += '<table style="max-width: 1000px; width: 100%; border-collapse: collapse; line-height: 0.8;">' + "\n"
        mdtext += "<thead>\n"
        mdtext += "<tr>\n"
        mdtext += "<th>Dataset ID</th>\n"
        mdtext += "<th># regions</th>\n"
        mdtext += "<th># hit regions</th>\n"
        mdtext += "<th>% hit regions</th>\n"
        mdtext += "<th># regex hits</th>\n"
        mdtext += "<th>p-value</th>\n"
        mdtext += "</tr>\n"
        mdtext += "</thead>\n"
        mdtext += "<tbody>\n"

        for internal_id in id2regex_stats_dic:
            rbp_id = id2infos_dic[internal_id][0]
            data_id = id2infos_dic[internal_id][1]
            method_id = id2infos_dic[internal_id][2]
            database_id = id2infos_dic[internal_id][3]

            combined_id = rbp_id + "," + method_id + "," + data_id
            if add_motif_db_info:
                combined_id = rbp_id + "," + database_id + "," + method_id + "," + data_id

            c_regex_hit_reg = id2regex_stats_dic[internal_id][0]
            c_regex_no_hit_reg = id2regex_stats_dic[internal_id][1]
            c_uniq_regex_hits = id2regex_stats_dic[internal_id][2]
            wc_pval = id2regex_stats_dic[internal_id][3]
            c_all_set_reg = c_regex_hit_reg + c_regex_no_hit_reg
            perc_hit_reg = (c_regex_hit_reg / c_all_set_reg) * 100
            perc_hit_reg = str(round(perc_hit_reg, 2))

            # mdtext += "| %s | %i | %i | %s | %i | %s |\n" %(combined_id, c_all_set_reg, c_regex_hit_reg, perc_hit_reg, c_uniq_regex_hits, str(wc_pval))

            mdtext += "<tr>\n"
            mdtext += "<td>%s</td>\n" %(combined_id)
            mdtext += "<td>%i</td>\n" %(c_all_set_reg)
            mdtext += "<td>%i</td>\n" %(c_regex_hit_reg)
            mdtext += "<td>%s</td>\n" %(perc_hit_reg)
            mdtext += "<td>%i</td>\n" %(c_uniq_regex_hits)
            mdtext += "<td>%s</td>\n" %(str(wc_pval))
            mdtext += "</tr>\n"

        mdtext += "</tbody>\n"
        mdtext += "</table>\n"

        mdtext += "\n&nbsp;\n&nbsp;\n"
        mdtext += "\nColumn IDs have the following meanings: "
        mdtext += "**Dataset ID** -> Dataset ID with following format: %s, " %(dataset_id_format)
        mdtext += '**# regions** -> number of input genomic regions in dataset (after filtering and optional extension), '
        mdtext += '**# hit regions** -> number of input genomic regions with regex hits (after filtering and optional extension), '
        mdtext += '**% hit regions** -> percentage of regex hit regions over all regions (i.e., how many input regions contain >= 1 regex motif hit), '
        mdtext += '**# regex hits** -> number of unique regex motif hits in input regions (removed double counts), '
        mdtext += '**p-value** -> Wilcoxon rank-sum test p-value.' + "\n"
        mdtext += "\n&nbsp;\n"


    """
    regex RBP motif co-occurrence statistics.
    
    """

    if id2regex_stats_dic:

        # Inform about set alterntive hypothesis for Fisher exact test on regex RBP motif co-occurrences.
        fisher_mode_info = "Fisher exact test alternative hypothesis is set to 'greater', i.e., low p-values mean that regex and RBP motifs have significantly high co-occurrence."
        if args.fisher_mode == 2:
            fisher_mode_info = "Fisher exact test alternative hypothesis is set to 'two-sided', i.e., low p-values mean that regex and RBP motifs have significantly high or low co-occurrence."
        elif args.fisher_mode == 3:
            fisher_mode_info = "Fisher exact test alternative hypothesis is set to 'less', i.e., low p-values mean that regex and RBP motifs have significantly low co-occurrence."

        mdtext += """
## Regular expression RBP motif co-occurrence statistics ### {#regex-rbp-cooc-stats}

**Table:** Motif hit co-occurrences (Fisher's exact test p-values) between 
regular expression (regex) '%s' and RBP motif hits for each input dataset.
Fisher's exact test p-value is calculated based on contingency table of co-occurrence 
counts (i.e., number of genomic regions with/without shared motif hits) 
between regex and RBP motif(s) for each dataset.
%s

""" %(args.regex, fisher_mode_info)

        # mdtext += '| Dataset ID  | contingency table | avg min distance | perc close hits |  p-value |' + " \n"
        # mdtext += '| :-: | :-: | :-: | :-: | :-: |' + " \n"

        mdtext += '<table style="max-width: 1200px; width: 100%; border-collapse: collapse; line-height: 0.8;">' + "\n"
        mdtext += "<thead>\n"
        mdtext += "<tr>\n"
        mdtext += "<th>Dataset ID</th>\n"
        mdtext += "<th>contingency table</th>\n"
        mdtext += "<th>avg min distance</th>\n"
        mdtext += "<th>% close hits</th>\n"
        mdtext += "<th>p-value</th>\n"
        mdtext += "</tr>\n"
        mdtext += "</thead>\n"
        mdtext += "<tbody>\n"

        for internal_id in id2regex_stats_dic:
            rbp_id = id2infos_dic[internal_id][0]
            data_id = id2infos_dic[internal_id][1]
            method_id = id2infos_dic[internal_id][2]
            database_id = id2infos_dic[internal_id][3]

            avg_min_dist = id2regex_stats_dic[internal_id][4]
            perc_close_hits = id2regex_stats_dic[internal_id][5]
            cont_table = id2regex_stats_dic[internal_id][6]
            fisher_pval = id2regex_stats_dic[internal_id][7]

            combined_id = rbp_id + "," + method_id + "," + data_id
            if add_motif_db_info:
                combined_id = rbp_id + "," + database_id + "," + method_id + "," + data_id

            # mdtext += "| %s | %s | %s | %s | %s |\n" %(combined_id, cont_table, avg_min_dist, perc_close_hits, fisher_pval)

            mdtext += "<tr>\n"
            mdtext += "<td>%s</td>\n" %(combined_id)
            mdtext += "<td>%s</td>\n" %(cont_table)
            mdtext += "<td>%s</td>\n" %(avg_min_dist)
            mdtext += "<td>%s</td>\n" %(perc_close_hits)
            mdtext += "<td>%s</td>\n" %(fisher_pval)
            mdtext += "</tr>\n"

        mdtext += "</tbody>\n"
        mdtext += "</table>\n"
        
        mdtext += "\n&nbsp;\n&nbsp;\n"

        perc_sign = "%"

        mdtext += """

Column IDs have the following meanings: 
**Dataset ID** -> Dataset ID with following format: %s,
**contingency table** -> contingency table of co-occurrence counts (i.e., number of genomic regions with/without shared regex + RBP motif hits), 
with format [[A, B], [C, D]], where 
A: regex AND RBP, 
B: NOT regex AND RBP
C: regex AND NOT RBP
D: NOT regex AND NOT RBP. 
**avg min distance** -> Mean minimum distance of regex and RBP motif hits (mean over all regions containing regex + RBP motif hits).
**%s close hits** -> Over all regions containing regex and RBP motif hit pairs, percentage of regions where regex + RBP motif hits are within %i nt distance (set via --max-motif-dist).
**p-value** -> Fisher's exact test p-value (calculated based on contingency table).

&nbsp;

""" %(dataset_id_format, perc_sign, args.max_motif_dist)


    """
    Input datasets nucleotide percentages table.

    """

    if seq_feat_ll:

        comp_info = "Mono-nucleotide contents (A, C, G, T) are used to calculate sequence complexity (change via --seq-comp-k). Equal A,C,G,T contents result in a complexity value of 1.0, while, e.g., an AA content of 100.0% results in a complexity of 0.0."
        if args.seq_comp_k == 2:
            comp_info = "Di-nucleotide contents (AA, AC, ..., TG, TT) are used to calculate sequence complexity (change via --seq-comp-k). Equal AA,AC,... contents result in a complexity value of 1.0, while, e.g., an AA content of 100.0% results in a complexity of 0.0."

        mdtext += """
## Input datasets nucleotide percentages statistics ### {#nt-perc-stats}

**Table:** Input datasets nucleotide percentages statistics for all input datasets.
Percentages of mono- and di-nucleotide contents (AC, AG, AT, CG, CT, GT) are given, as well as mean sequence complexity (Shannon entropy).

"""

        # mdtext += '| Dataset ID  | contingency table | avg min distance | perc close hits |  p-value |' + " \n"
        # mdtext += '| :-: | :-: | :-: | :-: | :-: |' + " \n"

        mdtext += '<table style="max-width: 1200px; width: 100%; border-collapse: collapse; line-height: 0.8;">' + "\n"
        mdtext += "<thead>\n"
        mdtext += "<tr>\n"
        mdtext += "<th>Dataset ID</th>\n"
        mdtext += "<th>% A</th>\n"
        mdtext += "<th>% C</th>\n"
        mdtext += "<th>% G</th>\n"
        mdtext += "<th>% T</th>\n"
        mdtext += "<th>% AC</th>\n"
        mdtext += "<th>% AG</th>\n"
        mdtext += "<th>% AT</th>\n"
        mdtext += "<th>% CG</th>\n"
        mdtext += "<th>% CT</th>\n"
        mdtext += "<th>% GT</th>\n"
        mdtext += "<th>Mean complexity</th>\n"
        mdtext += "</tr>\n"
        mdtext += "</thead>\n"
        mdtext += "<tbody>\n"

        for seq_feat_l in seq_feat_ll:
            dataset_id = seq_feat_l[0]
            seq_comp = str(seq_feat_l[1])
            perc_a = str(seq_feat_l[2])
            perc_c = str(seq_feat_l[3])
            perc_g = str(seq_feat_l[4])
            perc_t = str(seq_feat_l[5])
            perc_ac = str(seq_feat_l[6])
            perc_ag = str(seq_feat_l[7])
            perc_at = str(seq_feat_l[8])
            perc_cg = str(seq_feat_l[9])
            perc_ct = str(seq_feat_l[10])
            perc_gt = str(seq_feat_l[11])

            mdtext += "<tr>\n"
            mdtext += "<td>%s</td>\n" %(dataset_id)
            mdtext += "<td>%s</td>\n" %(perc_a)
            mdtext += "<td>%s</td>\n" %(perc_c)
            mdtext += "<td>%s</td>\n" %(perc_g)
            mdtext += "<td>%s</td>\n" %(perc_t)
            mdtext += "<td>%s</td>\n" %(perc_ac)
            mdtext += "<td>%s</td>\n" %(perc_ag)
            mdtext += "<td>%s</td>\n" %(perc_at)
            mdtext += "<td>%s</td>\n" %(perc_cg)
            mdtext += "<td>%s</td>\n" %(perc_ct)
            mdtext += "<td>%s</td>\n" %(perc_gt)
            mdtext += "<td>%s</td>\n" %(seq_comp)
            mdtext += "</tr>\n"

        mdtext += "</tbody>\n"
        mdtext += "</table>\n"
        
        mdtext += "\n&nbsp;\n&nbsp;\n"

        perc_sign = "%"

        mdtext += """

Column IDs have the following meanings: 
**Dataset ID** -> Dataset ID with following format: %s,
**%% A** -> percentage of A nucleotides in input regions,
**%% C** -> percentage of C nucleotides in input regions,
**%% G** -> percentage of G nucleotides in input regions,
**%% T** -> percentage of T nucleotides in input regions,
**%% AC** -> AC content percentage in input regions,
**%% AG** -> AG content percentage in input regions,
**%% AT** -> AT content percentage in input regions,
**%% CG** -> CG content percentage in input regions,
**%% CT** -> CT content percentage in input regions,
**%% GT** -> GT content percentage in input regions,
**Mean complexity** -> Mean sequence complexity (Shannon entropy) of input regions. 
%s

&nbsp;

""" %(dataset_id_format, comp_info)


    """
    Input datasets k-mer frequencies comparative plot.

    """

    mdtext += """
## Input datasets k-mer frequencies comparative plot ### {#kmer-comp-plot}

"""

    if len(dataset_ids_list) > 3:

        kmer_comp_plot_plotly =  "kmer_comparative_plot.plotly.html"
        kmer_comp_plot_plotly_out = plots_out_folder + "/" + kmer_comp_plot_plotly

        create_kmer_comp_plot_plotly(dataset_ids_list, kmer_list, kmer_freqs_ll, 
                                     seq_len_stats_ll,
                                     kmer_comp_plot_plotly_out,
                                     seq_feat_ll=seq_feat_ll,
                                     include_plotlyjs=include_plotlyjs,
                                     full_html=plotly_full_html)

        plot_path = plots_folder + "/" + kmer_comp_plot_plotly

        if args.plotly_js_mode in [5, 6, 7]:
            # Read in plotly code.
            # mdtext += '<div style="width: 1200px; height: 1200px; align-items: center;">' + "\n"
            js_code = read_file_content_into_str_var(kmer_comp_plot_plotly_out)
            js_code = js_code.replace("height:100%; width:100%;", "height:1000px; width:1200px;")
            mdtext += js_code + "\n"
            # mdtext += '</div>'
        else:
            if plotly_embed_style == 1:
                # mdtext += '<div class="container-fluid" style="margin-top:40px">' + "\n"
                mdtext += "<div>\n"
                mdtext += '<iframe src="' + plot_path + '" width="1200" height="1000"></iframe>' + "\n"
                mdtext += '</div>'
            elif plotly_embed_style == 2:
                mdtext += '<object data="' + plot_path + '" width="1200" height="1000"> </object>' + "\n"

        mdtext += """

**Figure:** Comparison of input datasets, using k-mer frequencies (k = %i, change via --kmer-size) of input region sequences 
(3-dimensional PCA) as features, to show similarities of input datasets based on their sequence k-mer frequencies 
(points close to each other have similar k-mer frequencies).
Input dataset IDs (show via hovering over data points) have following format: %s.

&nbsp;

""" %(args.kmer_size, dataset_id_format)

    else:
        mdtext += """

No plot generated since < 4 datasets were provided.

&nbsp;

"""

    """
    Input datasets sequence / k-mer variation comparative plot.

    """

    if seq_var_ll:

        assert seq_var_kmer_l, "no k-mer list provided for sequence variation plot"

        top_bottom_n = 10
        remove_zero_val = True
        color_mode = 1

        mdtext += """
## Input datasets k-mer variation comparative plot ### {#seq-var-plot}

"""
        if len(dataset_ids_list) > 3:

            seq_var_plot_plotly =  "seq_variation_plot.plotly.html"
            seq_var_plot_plotly_out = plots_out_folder + "/" + seq_var_plot_plotly

            if args.seq_var_feat_mode == 1:

                create_seq_var_plot_plotly_v2(seq_var_ll, seq_var_kmer_l,
                                        seq_var_plot_plotly_out,
                                        kmer_size=args.seq_var_kmer_size,
                                        top_bottom_n=top_bottom_n,
                                        color_mode=args.seq_var_color_mode,
                                        remove_zero_val=remove_zero_val,
                                        seq_len_stats_ll=seq_len_stats_ll,
                                        include_plotlyjs=include_plotlyjs,
                                        full_html=plotly_full_html)

            elif args.seq_var_feat_mode == 2:

                create_seq_var_plot_plotly(seq_var_ll, seq_var_kmer_l,
                                        seq_var_plot_plotly_out,
                                        kmer_size=args.seq_var_kmer_size,
                                        top_bottom_n=top_bottom_n,
                                        remove_zero_val=remove_zero_val,
                                        include_plotlyjs=include_plotlyjs,
                                        full_html=plotly_full_html)

            else:
                assert False, "invalid seq_var_feat_mode given"

            plot_path = plots_folder + "/" + seq_var_plot_plotly

            if args.plotly_js_mode in [5, 6, 7]:
                js_code = read_file_content_into_str_var(seq_var_plot_plotly_out)
                js_code = js_code.replace("height:100%; width:100%;", "height:1000px; width:1200px;")
                mdtext += js_code + "\n"
            else:
                if plotly_embed_style == 1:
                    mdtext += "<div>\n"
                    mdtext += '<iframe src="' + plot_path + '" width="1200" height="1000"></iframe>' + "\n"
                    mdtext += '</div>'
                elif plotly_embed_style == 2:
                    mdtext += '<object data="' + plot_path + '" width="1200" height="1000"> </object>' + "\n"

            if args.seq_var_feat_mode == 1:

                color_info = (
                    "Points are colored by the mean site percentage (change via --seq-var-color-mode), "
                    "which is the average of the site percentages of all present k-mers in the dataset. "
                    "A higher mean percentage for a dataset means that the present k-mers are more evenly distributed over the dataset sequences. "
                    "This can stem from a larger sequence lengths, or in general a more diverse set of sequences, "
                    "possibly reflecting RBP binding preferences."
                )                

                if args.seq_var_color_mode == 2:
                    color_info = (
                        "Points are colored by the percentage of present k-mers in the dataset (change via --seq-var-color-mode). "
                        "A percentage of 100%% means that all k-mers are present, while lower percentages can occur e.g. if k-mer size is high, "
                        "the number of dataset sequences is low, or if the sequences are less diverse (e.g. containing many repeat regions)."
                    )

                mdtext += r"""

**Figure:** Sequence k-mer variations for each input dataset visualized as 2D PCA plot.
For each k-mer, the percentage of sites containing the k-mer is used as feature (k = %i, change via --seq-var-kmer-size).
The closer two datasets in the PCA plot (i.e., the dots representing them), the more similar their k-mer site percentage profiles.
The hover box shows the k-mers ranked by site percentages (k-mers not present in any sites are not listed). It can be assumed 
that k-mers with high site percentages contribute more to the RBP binding (possibly also as direct binding motifs), 
while low percentage k-mers likely should also have low affinity to the RBP.
%s
Hover box content:
Dataset ID (bold-faced, format: %s),
sorted k-mer site percentages, 
**Present k-mer %%** -> percentage of all possible k-mers present in the dataset.
**Average k-mer %%** -> mean site percentage of all k-mers present in the dataset.
**# regions** -> number of input regions (i.e., sequences) in dataset.
**Median length** -> median length of input regions.
**Mean length** -> mean length of input regions.
**Min length** -> minimum input region length.
**Max length** -> maximum input region length.

&nbsp;

""" %(args.seq_var_kmer_size, color_info, dataset_id_format)

            elif args.seq_var_feat_mode == 2:

                mdtext += r"""

**Figure:** Sequence k-mer variations for each input dataset visualized as 2D PCA plot.
For each dataset, the variation of each k-mer is calculated over all sequences (k = %i, change via --seq-var-kmer-size).
The k-mer ratio is calculated for each sequence (number of times k-mer occurs / number of k-mers in sequence), 
resulting in a list of ratios for each k-mer.
The ratios are then used to calculate the standard deviation and mean for each k-mer.
Based on these, the coefficient of variation (CV) is calculated for each k-mer, 
which is defined as the standard deviation of ratios divided by the mean ratio.
For every dataset, the list of single k-mer CVs is used as PCA input for dimension reduction.
Thus, the closer two datasets in the PCA plot (i.e., the dots representing them), 
the more similar the k-mer variations in their sequences.
A CV of 0.0 for a given k-mer would mean that all dataset sequences have exactly the same ratio of the k-mer.
A CV of > 1.0 for given k-mer means that the standard deviation is larger than the mean, 
indicating a higher variation of the k-mer ratio in the sequences (although larger CVs are expected for larger k-mer sizes).
The lower a k-mer CV in a dataset, the more even its k-mer ratio over the dataset sequences.
In general, we can assume that k-mers with high CVs are less important for the RBP binding, since
they do not occur evenly over the dataset. In contrast, k-mers with low CVs are more 
evenly distributed over the dataset sequences, and thus should be more important for RBP binding 
(assuming a sufficient dataset quality + an affinity of the RBP towards specific k-mer sequences). 
The average CV (mean over all single k-mer CVs) is used for coloring.
**Hover box content:**
dataset ID (bold-faced, format: %s),
single k-mer CVs -> single k-mer CVs sorted by ascending CV (for k-mer sizes > 2 only top and bottom %i CVs are shown),
Average CV -> average CV ratio of the dataset.

&nbsp;

""" %(args.seq_var_kmer_size, dataset_id_format, top_bottom_n)


        else:
            mdtext += """

No plot generated since < 4 datasets were provided.

&nbsp;

"""








    """
    Input datasets gene / transcript region occupancies comparative plot.

    """
    c_total_number_genes = 0  # total number of gene regions considered for PCA plots.

    if id2reg_annot_dic:  # if --gtf provided.

        mdtext += """
## Input datasets occupied gene regions comparative plot ### {#occ-comp-plot}

"""

        if len(dataset_ids_list) > 3:

            occ_comp_plot_plotly =  "gene_reg_occ_comparative_plot.plotly.html"
            occ_comp_plot_plotly_out = plots_out_folder + "/" + occ_comp_plot_plotly

            # Get number of gene regions considered.
            for internal_id in id2occ_list_dic:
                c_total_number_genes = len(id2occ_list_dic[internal_id])
                break

            create_pca_reg_occ_plot_plotly(id2occ_list_dic, id2infos_dic,
                                           occ_comp_plot_plotly_out,
                                           sparse_pca=False,
                                           id2hk_gene_stats_dic=id2hk_gene_stats_dic,
                                           add_motif_db_info=add_motif_db_info,
                                           include_plotlyjs=include_plotlyjs,
                                           full_html=plotly_full_html)

            plot_path = plots_folder + "/" + occ_comp_plot_plotly

            if args.plotly_js_mode in [5, 6, 7]:
                # Read in plotly code.
                # mdtext += '<div style="width: 1200px; height: 1200px; align-items: center;">' + "\n"
                js_code = read_file_content_into_str_var(occ_comp_plot_plotly_out)
                js_code = js_code.replace("height:100%; width:100%;", "height:1000px; width:1200px;")
                mdtext += js_code + "\n"
                # mdtext += '</div>'
            else:
                if plotly_embed_style == 1:
                    # mdtext += '<div class="container-fluid" style="margin-top:40px">' + "\n"
                    mdtext += "<div>\n"
                    mdtext += '<iframe src="' + plot_path + '" width="1200" height="1000"></iframe>' + "\n"
                    mdtext += '</div>'
                elif plotly_embed_style == 2:
                    mdtext += '<object data="' + plot_path + '" width="1200" height="1000"> </object>' + "\n"

            mdtext += """

**Figure:** Comparison of input datasets, using gene region occupancy of input regions (3-dimensional PCA) as features, 
to show similarities of input datasets based on the gene regions they occupy (points close to each other have similar occupancy profiles).
Note that only gene regions covered by regions in any of the input datasets are considered (# of gene regions: %i).
Datasets are colored by the percentage of the total %i gene regions they occupy.
Input dataset IDs (show via hovering over data points) have following format: %s.

&nbsp;

""" %(c_total_number_genes, c_total_number_genes, dataset_id_format)

        else:
            mdtext += """

No plot generated since < 4 datasets were provided.

&nbsp;

"""



    """
    Input datasets gene / transcript region occupancies co-occurrence heat map plot.

    """

    if id2reg_annot_dic and gene_occ_cooc_plot:

        mdtext += """
## Input datasets occupied gene regions similarity heat map ### {#cooc-heat-map}

"""

        # Get tables for plotting heatmap.
        df_gene_cooc, gene_cooc_lll = get_gene_occ_cooc_tables(id2occ_list_dic, id2infos_dic,
                                                               optimal_ordering=heatmap_cluster_olo)

        cooc_plot_plotly =  "gene_occ_cooc_heatmap.plotly.html"
        cooc_plot_plotly_out = plots_out_folder + "/" + cooc_plot_plotly

        create_gene_occ_cooc_plot_plotly(df_gene_cooc, gene_cooc_lll, cooc_plot_plotly_out,
                                        include_plotlyjs=include_plotlyjs,
                                        full_html=plotly_full_html)


        plot_path = plots_folder + "/" + cooc_plot_plotly

        if args.plotly_js_mode in [5, 6, 7]:
            # Read in plotly code.
            # mdtext += '<div style="width: 1200px; height: 1200px; align-items: center;">' + "\n"
            js_code = read_file_content_into_str_var(cooc_plot_plotly_out)
            js_code = js_code.replace("height:100%; width:100%;", "height:1200px; width:1200px;")
            mdtext += js_code + "\n"
            # mdtext += '</div>'
        else:
            if plotly_embed_style == 1:
                # mdtext += '<div class="container-fluid" style="margin-top:40px">' + "\n"
                mdtext += "<div>\n"
                mdtext += '<iframe src="' + plot_path + '" width="1200" height="1200"></iframe>' + "\n"
                mdtext += '</div>'
            elif plotly_embed_style == 2:
                mdtext += '<object data="' + plot_path + '" width="1200" height="1200"> </object>' + "\n"


        mdtext += """
**Figure:** Heat map comparing gene region occupancy profiles between input datasets. Similarity between the occupied 
gene regions of two datasets is calculated using the cosine similarity of the gene region binary vectors. A cosine similarity 
of 1 indicates that two datasets have identical occupancy profiles (like datasets compared with themselves on diagonal), 
while a cosine similarity of 0 means that the two vectors are completely different (i.e., if dataset 1 covers a gene region, 
dataset 2 does not and vice versa). 
Datasets are ordered based on hierarchical clustering of the cosine similarity values, 
so that similar datasets are placed close to each other.
Hover box: 
**1)** dataset 1 ID.
**2)** dataset 2 ID.
**3)** Cosine similarity value between the two datasets.
**4)** Counts[]: contingency table of co-occurrence counts (i.e., number of gene regions covered/not covered by the two datasets), 
with format [[A, B], [C, D]], where 
A: dataset 1 AND dataset 2, 
B: NOT dataset1 AND dataset2,
C: dataset1 AND NOT dataset2,
D: NOT dataset1 AND NOT dataset2.
Gene regions are labelled 1 or 0 (1 if a genomic region from the dataset overlaps with the gene region, otherwise 0), 
resulting in a vector of 1s and 0s for each RBP, which is then used to construct the contigency table.
Note that only gene regions covered by regions in any of the input datasets are considered (# of gene regions: %i).

&nbsp;

""" %(c_total_number_genes)




    """
    Input datasets region annotations comparative plot.

    """

    # Format: id2reg_annot_dic[internal_id][annot] = count
    if id2reg_annot_dic:  # if --gtf provided.


        mdtext += """
## Input datasets genomic region annotations comparative plot ### {#annot-comp-plot}

"""
        # If --regex + --gtf given, regex_annot_dic contains regex hit region annotations.
        regex_info = ""
        if regex_annot_dic:
            for annot in regex_annot_dic:
                if annot not in annot_dic:
                    annot_dic[annot] = 1

        annot_ids_list = []
        for annot in sorted(annot_dic, reverse=True):
            annot_ids_list.append(annot)

        annot_dataset_ids_list = []
        annot_freqs_ll = []
        for internal_id in id2reg_annot_dic:
            annot_freqs_list = []
            sum_annot = 0
            rbp_id, data_id, method_id, motif_db_str, bed_file_path = id2infos_dic[internal_id]

            if add_motif_db_info:
                dataset_id = rbp_id + "," + motif_db_str + "," + method_id + "," + data_id
            else:
                dataset_id = rbp_id + "," + method_id + "," + data_id

            for annot in sorted(annot_dic, reverse=True):
                c_annot = 0
                if annot in id2reg_annot_dic[internal_id]:
                    c_annot = id2reg_annot_dic[internal_id][annot]
                annot_freqs_list.append(c_annot)
                sum_annot += c_annot
            
            # Normalize counts in list by sum.
            annot_freqs_list = [x/sum_annot for x in annot_freqs_list]
            # Append to list of lists.
            annot_freqs_ll.append(annot_freqs_list)
            annot_dataset_ids_list.append(dataset_id)

        if regex_annot_dic:
            dataset_id = "regex:" + args.regex
            regex_info = " Genomic region annotations for regex hit regions (regex: %s) in all input datasets are also included." %(args.regex)

            annot_freqs_list = []
            sum_annot = 0

            for annot in sorted(annot_dic, reverse=True):
                c_annot = 0
                if annot in regex_annot_dic:
                    c_annot = regex_annot_dic[annot]
                annot_freqs_list.append(c_annot)
                sum_annot += c_annot
            
            annot_freqs_list = [x/sum_annot for x in annot_freqs_list]
            annot_freqs_ll.append(annot_freqs_list)
            annot_dataset_ids_list.append(dataset_id)

        if len(annot_dataset_ids_list) > 2:

            annot_comp_plot_plotly =  "gen_annot_comparative_plot.plotly.html"
            annot_comp_plot_plotly_out = plots_out_folder + "/" + annot_comp_plot_plotly

            create_annot_comp_plot_plotly(annot_dataset_ids_list, annot_freqs_ll,
                                          annot_ids_list, annot2color_dic,
                                          annot_comp_plot_plotly_out,
                                          include_plotlyjs=include_plotlyjs,
                                          full_html=plotly_full_html)

            plot_path = plots_folder + "/" + annot_comp_plot_plotly

            if args.plotly_js_mode in [5, 6, 7]:
                # Read in plotly code.
                # mdtext += '<div style="width: 1200px; height: 1200px; align-items: center;">' + "\n"
                js_code = read_file_content_into_str_var(annot_comp_plot_plotly_out)
                js_code = js_code.replace("height:100%; width:100%;", "height:1000px; width:1200px;")
                mdtext += js_code + "\n"
                # mdtext += '</div>'
            else:
                if plotly_embed_style == 1:
                    # mdtext += '<div class="container-fluid" style="margin-top:40px">' + "\n"
                    mdtext += "<div>\n"
                    mdtext += '<iframe src="' + plot_path + '" width="1200" height="1000"></iframe>' + "\n"
                    mdtext += '</div>'
                elif plotly_embed_style == 2:
                    mdtext += '<object data="' + plot_path + '" width="1200" height="1000"> </object>' + "\n"

            mdtext += """

**Figure:** Comparison of input datasets, using genomic region annotations from GTF file of input regions as features, 
to show similarities between input datasets based on similar genomic region occupancy (see detailed annotations for each input dataset below).
Input dataset IDs (show via hovering over data points) have following format: %s.%s

&nbsp;

""" %(dataset_id_format, regex_info)

        else:
            mdtext += """

No plot generated since < 4 datasets were provided.

&nbsp;

"""


    """
    Region annotation plots for each dataset.

    """

    if id2reg_annot_dic:  # if --gtf provided.

        # Generate plotting sections.
        for idx, internal_id in enumerate(internal_ids_list):
            rbp_id = id2infos_dic[internal_id][0]
            data_id = id2infos_dic[internal_id][1]
            method_id = id2infos_dic[internal_id][2]
            database_id = id2infos_dic[internal_id][3]

            if add_motif_db_info:
                combined_id = rbp_id + "," + database_id + "," + method_id + "," + data_id
            else:
                combined_id = rbp_id + "," + method_id + "," + data_id

            annot_stacked_bars_plot =  "annotation_stacked_bars_plot.%i.png" %(idx)
            annot_stacked_bars_plot_out = plots_out_folder + "/" + annot_stacked_bars_plot

            mdtext += """
## %s region annotations ### {#annot-plot-%i}

""" %(combined_id, idx)

            # Check if no regions in dataset (i.e. no annotations).
            if not c_annot:
                mdtext += """

No plot generated since no regions for plotting.

&nbsp;

"""

            else:

                # print("rbp2regidx_dic:")
                # print(rbp2regidx_dic)
                # print("reg_ids_list:")
                # print(reg_ids_list)
                # print("reg2annot_dic:")
                # print(reg2annot_dic)

                create_batch_annotation_stacked_bars_plot(internal_id, id2infos_dic, id2reg_annot_dic, id2hit_reg_annot_dic,
                                                          annot2color_dic, annot_stacked_bars_plot_out,
                                                          plot_pdf=args.plot_pdf)

                plot_path = plots_folder + "/" + annot_stacked_bars_plot

                mdtext += '<img src="' + plot_path + '" alt="Annotation stacked bars plot"' + "\n"
                # mdtext += 'title="Annotation stacked bars plot" width="800" />' + "\n"
                mdtext += 'title="Annotation stacked bars plot" />' + "\n"
                mdtext += """
**Figure:** Genomic region annotations for input dataset **%s** (dataset ID format: %s). 
Input regions are overlapped with genomic regions from GTF file and genomic region feature with highest overlap 
is assigned to each input region. "intergenic" feature means none of the used GTF region features overlap with the input region 
(minimum overlap amount controlled by --gtf-feat-min-overlap). 
By default, RBPBench for each gene in the GTF file selects the region features of the most prominent transcript (see manual for details).
**%s**: annotations for input regions containing %s motif hits. 
**All**: annotations for all input regions (with and without motif hits).

&nbsp;

""" %(combined_id, dataset_id_format, rbp_id, rbp_id)


    """
    Additional region annotations statistics.
    

    dsid2add_annot_stats_dic["general"] = {"c_genes": add_annot_stats_dic["c_genes"], 
                                           "c_promoters": add_annot_stats_dic["c_promoters"], 
                                           "c_filt_min_tr_len": add_annot_stats_dic["c_filt_min_tr_len"], 
                                           "c_filt_mrna_only": add_annot_stats_dic["c_filt_mrna_only"]}

    """

    if dsid2add_annot_stats_dic:
        
        add_annot_bed_info = ""
        add_annot_not_ol = " "
        if args.add_annot_bed:
            if args.add_annot_comp:
                add_annot_not_ol = " NOT "
            add_annot_bed_info = 'Also included are the percentages of input regions' + add_annot_not_ol + 'overlapping with additionally provided regions (via --add-annot-bed, ID: ' + '"' + args.add_annot_id + '").'

        only_mrna_info = ""
        if args.prom_mrna_only:
            only_mrna_info = "Only mRNA transcripts were used for promoter region extraction (# transcripts removed = %i)." %(dsid2add_annot_stats_dic["general"]["c_filt_mrna_only"])
        min_tr_len_info = ""
        if args.prom_min_tr_len:
            min_tr_len_info = "Minimum transcript length for promoter region extraction = %i nt (# transcripts removed = %i)." %(args.prom_min_tr_len, dsid2add_annot_stats_dic["general"]["c_filt_min_tr_len"])

        mdtext += """
## Additional region annotation statistics ### {#add-annot-stats}

**Table:** Percentages of input regions that overlap with additional region annotations are shown for each input dataset.
This includes percentages of input regions located outside of gene regions (annotated in provided GTF, # of considered
gene regions = %i) and input regions overlapping with putative promoter regions (taking the regions %i nt upstream 
to %i nt downstream of the transcript start sites). %s %s # of considered promoter regions = %i. %s High percentages 
of input regions located outside gene regions or 
inside promoter regions can point at dataset issues (assuming RBPs bind primarily to gene/transcript regions) 
or distinct protein functions (e.g., RBPs moonlighting as transcription factors). Note that depending 
on the methods used for dataset generation, input regions outside of gene regions might also have
been removed already.

""" %(dsid2add_annot_stats_dic["general"]["c_genes"], dsid2add_annot_stats_dic["general"]["prom_ext_up"], dsid2add_annot_stats_dic["general"]["prom_ext_down"], only_mrna_info, min_tr_len_info, dsid2add_annot_stats_dic["general"]["c_promoters"], add_annot_bed_info)

        mdtext += '<table style="max-width: 1000px; width: 100%; border-collapse: collapse; line-height: 0.8;">' + "\n"
        mdtext += "<thead>\n"
        mdtext += "<tr>\n"
        mdtext += "<th>Dataset ID</th>\n"
        mdtext += "<th># input regions</th>\n"
        mdtext += "<th>% outside genes</th>\n"
        mdtext += "<th>% inside promoters</th>\n"
        if args.add_annot_bed:
            mdtext += "<th>% " + args.add_annot_id + "</th>\n"
        mdtext += "</tr>\n"
        mdtext += "</thead>\n"
        mdtext += "<tbody>\n"

        for internal_id in sorted(dsid2add_annot_stats_dic):
            if internal_id == "general":
                continue
            rbp_id = id2infos_dic[internal_id][0]
            data_id = id2infos_dic[internal_id][1]
            method_id = id2infos_dic[internal_id][2]
            database_id = id2infos_dic[internal_id][3]
            combined_id = rbp_id + "," + method_id + "," + data_id  # Dataset ID

            c_outside_genes = dsid2add_annot_stats_dic[internal_id]["c_outside_genes"]
            c_inside_prom = dsid2add_annot_stats_dic[internal_id]["c_inside_prom"]
            c_add_annot = dsid2add_annot_stats_dic[internal_id]["c_add_annot"]
            c_in_regions = dsid2add_annot_stats_dic[internal_id]["c_regions"]

            perc_outside_genes = 0.0
            if c_outside_genes > 0:
                perc_outside_genes = round( (c_outside_genes / c_in_regions) * 100, 1)
            perc_inside_prom = 0.0
            if c_inside_prom > 0:
                perc_inside_prom = round( (c_inside_prom / c_in_regions) * 100, 1)

            mdtext += '<tr>' + "\n"
            mdtext += "<td>" + str(combined_id) + "</td>\n"
            mdtext += "<td>" + str(c_in_regions) + "</td>\n"
            mdtext += "<td>" + str(perc_outside_genes) + "</td>\n"
            mdtext += "<td>" + str(perc_inside_prom) + "</td>\n"
            if args.add_annot_bed:
                perc_add_annot = 0.0
                if c_add_annot > 0:
                    perc_add_annot = round( (c_add_annot / c_in_regions) * 100, 1)
                mdtext += "<td>" + str(perc_add_annot) + "</td>\n"
            mdtext += '</tr>' + "\n"

        mdtext += '</tbody>' + "\n"
        mdtext += '</table>' + "\n"

        mdtext += "\n&nbsp;\n&nbsp;\n"
        mdtext += "\nColumn IDs have the following meanings: "
        mdtext += "**Dataset ID** -> Dataset ID for input dataset with following format: %s, " %(dataset_id_format)
        mdtext += '**# input regions** -> number of considered input regions from input dataset, '
        mdtext += "**% outside genes** -> percentage of input regions located outside of gene regions, "
        mdtext += '**% inside promoters** -> percentage of input regions overlapping with promoter regions, '
        if args.add_annot_bed:
            mdtext += "**% " + args.add_annot_id + "** -> percentage of input regions" + add_annot_not_ol + 'overlapping with additionally provided regions (via --add-annot-bed, ID: "' + args.add_annot_id + '"), '
        mdtext += '**% input regions** -> percentage of input regions overlapping or not overlapping with region annotations.'  + "\n"
        mdtext += "\n&nbsp;\n"


    """
    GOA results.

    """

    if args.run_goa:

        mdtext += """
## GO enrichment analysis results ### {#goa-results}

"""
        c_goa_results = 0
        if isinstance(goa_results_df, pd.DataFrame) and not goa_results_df.empty:
            c_goa_results = len(goa_results_df)

        filter_purified_info = " GO terms with significantly higher and lower concentration ([e,p]) in study group are shown."
        filter_purified_info2 = "significant"
        if args.goa_filter_purified:
            filter_purified_info = " Only GO terms with significantly higher concentration in study group are shown."
            filter_purified_info2 = "significantly enriched"
            c_goa_results = len(goa_results_df[goa_results_df["enrichment"] == "e"])
        filter_further_info = ""
        if args.goa_max_child is not None: 
            filter_further_info += " Only GO terms with <= %i children are shown." %(args.goa_max_child)
        if args.goa_min_depth is not None:
            filter_further_info += " Only GO terms with >= %i depth are shown." %(args.goa_min_depth)
        if filter_further_info:
            filter_further_info += " Note that additional filters (children + depth) can result in an empty table. For all significant GO terms (i.e., unfiltered results) check *goa_results.tsv* output table."
        filter_only_cooc_info = "Only target genes are considered which are covered by regions from all input datasets."
        if args.goa_only_cooc:
            filter_only_cooc_info = " Only target genes are considered which are covered by regions with motif hits from all input datasets (--goa-only-cooc enabled)."

        if c_goa_results > 0:

            mdtext += """
**Table:** GO enrichment analysis results. # of %s GO terms found: %i. Filter p-value threshold (on corrected p-value) = %s. # of target genes used for GOA: %i. # of background genes used for GOA: %i.
%s %s %s

""" %(filter_purified_info2, c_goa_results, str(goa_stats_dic["pval_thr"]), goa_stats_dic["c_target_genes_goa"], goa_stats_dic["c_background_genes_goa"], filter_only_cooc_info, filter_purified_info, filter_further_info)

            mdtext += '<table style="max-width: 1200px; width: 100%; border-collapse: collapse; line-height: 0.9;">' + "\n"
            mdtext += "<thead>\n"
            mdtext += "<tr>\n"
            mdtext += "<th>GO</th>\n"
            mdtext += "<th>Term</th>\n"
            mdtext += "<th>Class</th>\n"
            mdtext += "<th>p-value</th>\n"
            mdtext += "<th>[e,p]</th>\n"
            mdtext += "<th>Depth</th>\n"
            mdtext += "<th># child</th>\n"
            mdtext += "<th># genes</th>\n"
            mdtext += "<th># study</th>\n"
            mdtext += "<th>% genes</th>\n"
            mdtext += "</tr>\n"
            mdtext += "</thead>\n"
            mdtext += "<tbody>\n"

            for index, row in goa_results_df.iterrows():

                go_id = row['GO']
                go_term = row['term']
                go_class = row['class']
                # go_p = row['p']
                go_p_corr = row['p_corr']
                go_enrichment = row['enrichment']
                go_depth = row['depth']
                go_n_genes = row['n_genes']
                go_n_study = row['n_study']
                go_perc_genes = row['perc_genes']
                go_n_children = row['n_children']

                if args.goa_filter_purified:
                    if go_enrichment == "p":
                        continue

                if args.goa_max_child is not None:
                    if go_n_children > args.goa_max_child:
                        continue
                if args.goa_min_depth is not None:
                    if go_depth < args.goa_min_depth:
                        continue

                mdtext += '<tr>' + "\n"
                mdtext += "<td>" + go_id + "</td>\n"
                mdtext += "<td>" + go_term + "</td>\n"
                mdtext += "<td>" + go_class + "</td>\n"
                mdtext += "<td>" + str(go_p_corr) + "</td>\n"
                mdtext += "<td>" + go_enrichment + "</td>\n"
                mdtext += "<td>" + str(go_depth) + "</td>\n"
                mdtext += "<td>" + str(go_n_children) + "</td>\n"
                mdtext += "<td>" + str(go_n_genes) + "</td>\n"
                mdtext += "<td>" + str(go_n_study) + "</td>\n"
                mdtext += "<td>" + str(go_perc_genes) + "</td>\n"
                mdtext += '</tr>' + "\n"

            mdtext += '</tbody>' + "\n"
            mdtext += '</table>' + "\n"
            
            mdtext += "\n&nbsp;\n&nbsp;\n"
            mdtext += "\nColumn IDs have the following meanings: "
            mdtext += "**GO** -> gene ontology (GO) ID, "
            mdtext += "**Term** -> GO term / name, "
            mdtext += "**Class** -> GO term class (biological_process, molecular_function, or cellular_component), "
            mdtext += "**p-value** -> multiple testing corrected (BH) p-value, "
            mdtext += "**[e,p]** -> e: enriched, i.e., GO term with significantly higher concentration, p: purified, GO term with significantly lower concentration), "
            mdtext += "**Depth** -> depth / level of GO term in GO hierarchy (the higher number, the more specific), "
            mdtext += "**# child** -> number of GO term children, "
            mdtext += "**# genes** -> number of genes associated with GO term, "
            mdtext += "**# study** -> number of genes in study (i.e., target genes), "
            mdtext += "**% genes** -> percentage of study genes associated with GO term." + "\n"
            mdtext += "\n&nbsp;\n"

        else:

            if "c_target_genes_goa" in goa_stats_dic:

                mdtext += """

No %s GO terms found given p-value threshold of %s. # of target genes used for GOA: %i. # of background genes used for GOA: %i.

&nbsp;

""" %(filter_purified_info2, str(goa_stats_dic["pval_thr"]), goa_stats_dic["c_target_genes_goa"], goa_stats_dic["c_background_genes_goa"])

            else:

                mdtext += """

No significant GO terms found due to no GO IDs associated with target genes. # of initial target genes (i.e., genes overlapping with --in regions): %i.

&nbsp;

""" %(goa_stats_dic["c_target_genes_pre_filter"])




    # Convert mdtext to html.
    md2html = markdown(mdtext, extensions=['attr_list', 'tables'])

    # OUTMD = open(md_out,"w")
    # OUTMD.write("%s\n" %(mdtext))
    # OUTMD.close()

    html_content = html_head + md2html + html_tail

    OUTHTML = open(html_out,"w")
    OUTHTML.write("%s\n" %(html_content))
    OUTHTML.close()


################################################################################

def create_mrna_region_occ_plot(motif_ids_list, mrna_reg_occ_dic, 
                                annot2color_dic, plot_out,
                                same_y_scale=True,
                                plot_pdf=False,
                                rbp_id=False):
    """
    Create mRNA region occupancy stacked line plot for rbp_id and associated 
    motif IDs.
    
    mrna_reg_occ_dic:
        mRNA region occupancy dictionary for rbp_id and motif IDs.
        rbp_id/motif_id -> mrna_region -> positional counts list

    motif_ids_list:
        List of motif IDs to plot mRNA occupancy profiles for.

    """

    datasets = {}
    max_pos_count = 0  # Get maximum positional count.
    if rbp_id and len(motif_ids_list) > 1:
        datasets[rbp_id] = mrna_reg_occ_dic[rbp_id]
        for mrna_reg in mrna_reg_occ_dic[rbp_id]:
            max_pos_count = max(max_pos_count, max(mrna_reg_occ_dic[rbp_id][mrna_reg]))
    for motif_id in motif_ids_list:
        datasets[motif_id] = mrna_reg_occ_dic[motif_id]
        for mrna_reg in mrna_reg_occ_dic[motif_id]:
            max_pos_count = max(max_pos_count, max(mrna_reg_occ_dic[motif_id][mrna_reg]))

    # Number of datasets
    num_datasets = len(datasets)

    # Create a figure and subplots
    fig, axs = plt.subplots(nrows=num_datasets, sharex=True, figsize=(12, 2 * num_datasets))

    # Check if axs is an array or not (not if there's only one subplot).
    if num_datasets == 1:
        axs = [axs]

    utr5color = annot2color_dic["5'UTR"]
    cds_color = annot2color_dic["CDS"]
    utr3color = annot2color_dic["3'UTR"]

    # Plot each dataset. label: rbp_id/motif_id, data: positional counts list for each mRNA region.
    for ax, (label, data) in zip(axs, datasets.items()):

        # Concatenate data for plotting
        all_counts = data["5'UTR"] + data["CDS"] + data["3'UTR"]
        x_positions = np.arange(len(all_counts))
        
        # Fill the area under the line for each region with different colors
        ax.fill_between(x_positions[:len(data["5'UTR"])+1], all_counts[:len(data["5'UTR"])+1], color=utr5color, alpha=1, zorder=3)
        ax.fill_between(x_positions[len(data["5'UTR"]):len(data["5'UTR"]) + len(data["CDS"])+1], all_counts[len(data["5'UTR"]):len(data["5'UTR"]) + len(data["CDS"])+1], color=cds_color, alpha=1, zorder=3)
        ax.fill_between(x_positions[-len(data["3'UTR"]):], all_counts[-len(data["3'UTR"]):], color=utr3color, alpha=1, zorder=3)
        
        if same_y_scale:
            if max_pos_count > 0:
                ax.set_ylim(0, max_pos_count)

        # Use dataset ID as y-axis label
        # label_y = label + ' motif coverage'
        ax.set_ylabel(label)

        # Remove x-axis ticks for each subplot
        ax.set_xticks([])
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['bottom'].set_visible(False)

        # Enable horizontal grid lines in light gray
        ax.grid(axis='y', color='lightgray', linestyle='-', linewidth=0.7, zorder=1)

    # Set common x-axis label
    axs[-1].set_xlabel('mRNA regions')

    # Create a single legend for all plots, positioned more optimally
    fig.legend(['5\'UTR', 'CDS', '3\'UTR'], loc='upper right', bbox_to_anchor=(0.975, 0.975))

    # Reduce padding and adjust subplots to fit the layout, minimizing vertical space
    plt.subplots_adjust(left=0.08, right=0.9, top=0.95, bottom=0.05, hspace=0.1)  # Adjusted hspace here

    # # Show plot
    # plt.show()

    if plot_pdf and plot_out.endswith('.png'):
        pdf_out = plot_out[:-4] + '.pdf'
        plt.savefig(pdf_out, dpi=110, bbox_inches='tight')

    plt.savefig(plot_out, dpi=110, bbox_inches='tight')  # 110
    plt.close()


################################################################################

def create_annotation_stacked_bars_plot(rbp_id, rbp2motif2annot2c_dic, 
                                        annot2color_dic, plot_out,
                                        x_label="Annotation overlap",
                                        no_x_labels=False,
                                        plot_pdf=False,
                                        y_label=""):
    """
    Do the motif hit genomic region annotations stacked bars plot.
    Plot all bars in one plot, one for all rbp_id hits, and one for each motif_id hits. 
    
    rbp2motif2annot2c_dic:
        rbp_id -> motif_id -> annot -> annot_c.
    annot2color_dic:
        annot -> hex color for plotting.

    """
    # Determine height of plot.
    c_height = 1  # rbp_id plot is set.
    motifs_ids_with_hits = []
    for motif_id in rbp2motif2annot2c_dic[rbp_id]:
        if rbp2motif2annot2c_dic[rbp_id][motif_id]:
            motifs_ids_with_hits.append(motif_id)
            c_height += 1

    # If only one motif, just plot the RBP.
    if len(motifs_ids_with_hits) == 1:
        c_height = 1

    fheight = 0.8*c_height
    fwidth = 10

    data_dic = {}
    # First sum up all annotations for rbp_id.
    for motif_id in motifs_ids_with_hits:
        for annot in rbp2motif2annot2c_dic[rbp_id][motif_id]:
            annot_c = rbp2motif2annot2c_dic[rbp_id][motif_id][annot]
            if annot not in data_dic:
                data_dic[annot] = [annot_c]
            else:
                data_dic[annot][0] += annot_c

    # Now for the single motif IDs.
    if len(motifs_ids_with_hits) > 1:
        for motif_id in motifs_ids_with_hits:
            for annot in data_dic:
                if annot in rbp2motif2annot2c_dic[rbp_id][motif_id]:
                    annot_c = rbp2motif2annot2c_dic[rbp_id][motif_id][annot]
                    data_dic[annot].append(annot_c)
                else:
                    data_dic[annot].append(0)

        data_dic["rbp_id"] = [rbp_id]
        for motif_id in motifs_ids_with_hits:
            data_dic["rbp_id"].append(motif_id)
    else:
        data_dic["rbp_id"] = [motifs_ids_with_hits[0]]

    df = pd.DataFrame(data_dic)
    # Reverse plotting order for bars.
    df = df.sort_values('rbp_id', ascending=False)

    ax = df.set_index('rbp_id').plot(kind='barh', stacked=True, legend=False, color=annot2color_dic, edgecolor="none", figsize=(fwidth, fheight))

    plt.xlabel(x_label)
    ax.set_ylabel(y_label)
    ax.yaxis.grid(False)
    ax.xaxis.grid(True)
    ax.set_axisbelow(True)

    if no_x_labels:
        ax.set_xticklabels([])
        ax.tick_params(axis='x', which='both', bottom=False, top=False)

    # Remove border lines.
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)

    # plt.legend(loc=(1.01, 0.4), fontsize=12, framealpha=0)
    # plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.5)
    plt.legend(bbox_to_anchor=(1.01, 1), loc='upper left')

    if plot_pdf and plot_out.endswith('.png'):
        pdf_out = plot_out[:-4] + '.pdf'
        plt.savefig(pdf_out, dpi=110, bbox_inches='tight')

    plt.savefig(plot_out, dpi=110, bbox_inches='tight')
    plt.close()


################################################################################

def create_batch_annotation_stacked_bars_plot(internal_id, id2infos_dic, id2reg_annot_dic, id2hit_reg_annot_dic,
                                              annot2color_dic, plot_out,
                                              plot_pdf=False):
    """
    Created stacked genomic region annotation plot.

    """
    fheight = 0.8*2
    fwidth = 10

    rbp_id = id2infos_dic[internal_id][0]

    data_dic = {}
    for annot in sorted(id2reg_annot_dic[internal_id], reverse=True):
        data_dic[annot] = []
        data_dic[annot].append(id2reg_annot_dic[internal_id][annot])
        if annot in id2hit_reg_annot_dic[internal_id]:
            data_dic[annot].append(id2hit_reg_annot_dic[internal_id][annot])
        else:
            data_dic[annot].append(0)
    data_dic["rbp_id"] = ["All", rbp_id]

    df = pd.DataFrame(data_dic)

    # print(internal_id)
    # print(df)

    ax = df.set_index('rbp_id').plot(kind='barh', stacked=True, legend=False, color=annot2color_dic, edgecolor="none", figsize=(fwidth, fheight))

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

    plt.savefig(plot_out, dpi=110, bbox_inches='tight')

    if plot_pdf and plot_out.endswith('.png'):
        pdf_out = plot_out[:-4] + '.pdf'
        plt.savefig(pdf_out, dpi=110, bbox_inches='tight')

    plt.close()


################################################################################

def filter_rbp2regidx_dic(rbp2regidx_dic,
                          min_rbp_count=0,
                          max_rbp_rank=None):
    """
    Filter RBP to region index dictionary by count and rank.

    """
    # Go through RBPs, and count number of sites for each RBP.
    max_count = 0
    rbp_count = 0
    rbp2count_dic = {}
    for rbp_id in rbp2regidx_dic:
        rbp_count += 1
        count = len(rbp2regidx_dic[rbp_id])
        rbp2count_dic[rbp_id] = count
        if count > max_count:
            max_count = count

    # Sort RBPs by count, and keep the max_rank top counting RBP IDs in list.
    sorted_rbp2count_list = sorted(rbp2count_dic.items(), key=lambda x: x[1], reverse=True)
    max_rank = rbp_count
    if max_rbp_rank is not None and max_rbp_rank < max_rank:
        max_rank = max_rbp_rank
    top_rbp2count_list = sorted_rbp2count_list[:max_rbp_rank]

    # Remove RBPs with count less than min_count.
    # if min_rbp_count > max_count:
    #     min_rbp_count = max_count
    top_rbp2count_list = [x for x in top_rbp2count_list if x[1] >= min_rbp_count]

    # print("top_rbp2count_list after removing RBPs with count less than %d:" %(min_rbp_count))
    # print(top_rbp2count_list)

    # Remove entries from rbp2regidx_dic for RBPs not in top_rbp2count_list.
    top_rbp2count_dic = dict(top_rbp2count_list)
    rbp2regidx_dic = {k: v for k, v in rbp2regidx_dic.items() if k in top_rbp2count_dic}

    # print("rbp2regidx_dic after removing RBPs not in top_rbp2count_list:")
    # print(rbp2regidx_dic)

    return rbp2regidx_dic


################################################################################

def create_seq_var_violin_plot_plotly(kmer2stats_dic, single_cv_dic, avg_cv, plot_out, 
                                      kmer_size=3,
                                      remove_zero_val=True,
                                      color_mode=1,
                                      include_plotlyjs="cdn",
                                      full_html=False):
    """
    Create sequence k-mer variation violin plot.

    color_mode:
        If 1, color by score correlation.
        If 2, color by k-mer percentage.


    """
    assert kmer2stats_dic, "kmer2stats_dic empty"
    assert single_cv_dic, "single_cv_dic empty"

    kmer_list = []
    for kmer in sorted(single_cv_dic):
        if remove_zero_val and single_cv_dic[kmer] == 0:
            continue
        kmer_list.append(kmer)
 
    # in_df = pd.DataFrame({
    #     'k-mer': kmer_list,
    #     'k-mer CV': [round(single_cv_dic[kmer], 5) for kmer in kmer_list],
    #     'k-mer count': [kmer2stats_dic[kmer][0] for kmer in kmer_list],
    #     'k-mer percentage': [round(kmer2stats_dic[kmer][1] * 100, 5) for kmer in kmer_list],
    #     'Site percentage' : [round(kmer2stats_dic[kmer][2] * 100, 3) for kmer in kmer_list],
    #     'Score correlation' : [kmer2stats_dic[kmer][3] for kmer in kmer_list]
    # })

    # fig = px.violin(
    #     in_df, 
    #     y='k-mer CV', 
    #     box=True, 
    #     points='all', 
    #     color='Score correlation', 
    # )
    # fig.update_traces(
    #     customdata=in_df[['k-mer', 'k-mer CV', 'k-mer count', 'k-mer percentage', 'Site percentage', 'Score correlation']].values,
    #     hovertemplate=(
    #         # '<span style="font-family: \'Courier New\', monospace;">>%{customdata[0]}</span><br>'
    #         # '<span style="font-family: \'Courier New\', monospace;">%{customdata[1]}</span><br>'
    #         # '<b>%{customdata[0]}</b><br>'
    #         '<span style="font-family: \'Courier New\', monospace;">'
    #         'k-mer:       <b>%{customdata[0]}</b><br>'
    #         'Variation:    %{customdata[1]}<br>'
    #         'k-mer #:     %{customdata[2]}<br>'
    #         'k-mer %:     %{customdata[3]}<br>'
    #         'Site %:      %{customdata[4]}<br>'
    #         'Correlation: %{customdata[5]}<br>'
    #         '</span><extra></extra>'
    #     ),
    #     line_color='#2b7bba',  # Set the color of the line
    #     fillcolor='rgba(43, 123, 186, 0.5)',  # Set the color of the filled area with transparency
    #     marker=dict(color='#2b7bba')  # Set the color of the points
    # )
    # fig.update_layout(xaxis_title='Density', yaxis_title='k-mer variation', violinmode='overlay')

    color_scale = ['#c6dbef', '#9ecae1', '#6baed6', '#4292c6', '#2171b5', '#08519c', '#08306b']
    hover_data = ['k-mer', 'k-mer CV', 'k-mer count', 'k-mer %', 'Site %', 'Correlation']
    color = 'Correlation'
    if color_mode == 2:
        color = 'k-mer %'

    # Create DataFrame
    in_df = pd.DataFrame({
        'k-mer': kmer_list,
        'k-mer CV': [round(single_cv_dic[kmer], 5) for kmer in kmer_list],
        'k-mer count': [kmer2stats_dic[kmer][0] for kmer in kmer_list],
        'k-mer %': [round(kmer2stats_dic[kmer][1] * 100, 5) for kmer in kmer_list],
        'Site %': [round(kmer2stats_dic[kmer][2] * 100, 3) for kmer in kmer_list],
        'Correlation': [kmer2stats_dic[kmer][3] for kmer in kmer_list]
    })

    # Create the scatter plot
    fig = px.scatter(
        in_df,
        x='k-mer CV',
        y='Site %',
        color=color,
        color_continuous_scale=color_scale,  # Set a color scale
        # hover_data=hover_data
    )

    fig.update_traces(
        customdata=in_df[hover_data].values,
        hovertemplate=(
            # '<span style="font-family: \'Courier New\', monospace;">>%{customdata[0]}</span><br>'
            # '<span style="font-family: \'Courier New\', monospace;">%{customdata[1]}</span><br>'
            # '<b>%{customdata[0]}</b><br>'
            '<span style="font-family: \'Courier New\', monospace;">'
            'k-mer:       <b>%{customdata[0]}</b><br>'
            'Variation:   %{customdata[1]}<br>'
            'k-mer #:     %{customdata[2]}<br>'
            'k-mer %:     %{customdata[3]}<br>'
            'Site %:      %{customdata[4]}<br>'
            'Correlation: %{customdata[5]}<br>'
            '</span><extra></extra>'
        ),
        # line_color='#2b7bba',  # Set the color of the line
        # fillcolor='rgba(43, 123, 186, 0.5)',  # Set the color of the filled area with transparency
        # marker=dict(color='#2b7bba')  # Set the color of the points
    )

    # fig.update_traces(
    #     hovertemplate=(
    #         '<span style="font-family: \'Courier New\', monospace;">'
    #         'k-mer:       <b>%{customdata[0]}</b><br>'
    #         'Variation:    %{customdata[1]}<br>'
    #         'k-mer #:     %{customdata[2]}<br>'
    #         'k-mer %:     %{customdata[3]}<br>'
    #         'Site %:      %{customdata[4]}<br>'
    #         'Correlation: %{customdata[5]}<br>'
    #         '</span><extra></extra>'
    #     )
    # )


    #     fig = px.scatter(
    #         df,
    #         x='PC1',
    #         y='PC2',
    #         color='Annotation',
    #         color_discrete_map=annot2color_dic,
    #         labels={
    #             'PC1': f'PC1 ({explained_variance[0]:.2f}% variance)',
    #             'PC2': f'PC2 ({explained_variance[1]:.2f}% variance)'
    #         },
    #         # hover_name='Region ID',
    #         hover_data=hover_data
    #     )


    # fig.update_traces(
    #     hovertemplate=(
    #         '<span style="font-family: \'Courier New\', monospace;">'
    #         'k-mer:       <b>%{customdata[0]}</b><br>'
    #         'Variation:    %{customdata[1]}<br>'
    #         'k-mer #:     %{customdata[2]}<br>'
    #         'k-mer %:     %{customdata[3]}<br>'
    #         'Site %:      %{customdata[4]}<br>'
    #         'Correlation: %{customdata[5]}<br>'
    #         '</span><extra></extra>'
    #     )
    # )

    # Adding '<extra></extra>' removes si
    # Update layout for titles and axis labels
    fig.update_layout(
        xaxis_title="k-mer variation",
        yaxis_title="Site %"
        # coloraxis_colorbar=dict(title="Score Correlation")  # Add color bar title
    )

    fig.write_html(plot_out,
                   full_html=full_html,
                   include_plotlyjs=include_plotlyjs)


################################################################################

def create_len_distr_violin_plot_plotly(in_df, plot_out,
                                        include_plotlyjs="cdn",
                                        full_html=False):
    """
    Create sequence lengths violin plot, including various infos in hover box,
    such as sequences and motif hits.
    
    in_df format created from seqs_dic:

    # Sequences dataframe for plotting sequence lengths violin plot.
    sequences = []
    seq_ids = []
    for seq_id in out_seqs_dic:
        seq_ids.append(seq_id)
        sequences.append(seqs_dic[seq_id])

    motif_hits = []
    for seq_id in seq_ids:
        # region_rbp_motif_pos_dic[seq_id].sort()
        motif_hits.append(";".join(region_rbp_motif_pos_dic[seq_id]))

    seq_len_df = pd.DataFrame({
        'Sequence ID': seq_ids,
        'Sequence Length': [len(seq) for seq in sequences],
        'Sequence': [benchlib.insert_line_breaks(seq, line_len=60) for seq in sequences],
        'Motif hits': motif_hits
    })

    """

    fig = px.violin(in_df, y='Sequence Length', box=True, points='all')
    # Set customdata and hovertemplate
    # fig.update_traces(customdata=in_df[['Sequence ID', 'Sequence', 'Motif hits']].values, hovertemplate='>%{customdata[0]}<br>%{customdata[1]}<br>Sequence Length: %{y}<br>Motif hits:<br>%{customdata[2]}')
    fig.update_traces(
        customdata=in_df[['Sequence ID', 'Sequence', 'Motif hits']].values,
        hovertemplate=(
            '<span style="font-family: \'Courier New\', monospace;">>%{customdata[0]}</span><br>'
            '<span style="font-family: \'Courier New\', monospace;">%{customdata[1]}</span><br>'
            # '>%{customdata[0]}<br>'
            # '%{customdata[1]}<br>'
            'Sequence Length: %{y}<br>'
            'Motif hits:<br>'
            '%{customdata[2]}'
        ),
        line_color='#2b7bba',  # Set the color of the line
        fillcolor='rgba(43, 123, 186, 0.5)',  # Set the color of the filled area with transparency
        marker=dict(color='#2b7bba')  # Set the color of the points
    )
    fig.update_layout(xaxis_title='Density', yaxis_title='Sequence Length', violinmode='overlay')
    fig.write_html(plot_out,
                   full_html=full_html,
                   include_plotlyjs=include_plotlyjs)


################################################################################

def get_sequence_length_statistics(seq_dic,
                                   unstranded=False):
    """
    Given a dictionary of sequences (key: sequence ID, value: sequence),
    calculate mean, median, q1 q3, min max for sequence lengths.

    unstranded:
        If True, seqs_dic contains both strands of a region, but should be counted 
        as one region. Format of sequence ID:  chr1:100-200(+), chr1:100-200(-)

    Return list of statistics.

    """
    
    assert seq_dic, "no sequences given for length statistics calculation"

    seen_ids_dic = {}

    seq_len_list = []
    for seq_id in seq_dic:
        if unstranded:
            core_id = reg_get_core_id(seq_id)
            if core_id in seen_ids_dic:
                continue
            seen_ids_dic[core_id] = 1
            seq_len_list.append(len(seq_dic[seq_id]))
        else:
            seq_len_list.append(len(seq_dic[seq_id]))

    nr_seqs = len(seq_len_list)
    seq_len_arr = np.array(seq_len_list)
    seq_len_mean = np.mean(seq_len_arr)
    seq_len_median = np.median(seq_len_arr)
    seq_len_q1 = np.percentile(seq_len_arr, 25)
    seq_len_q3 = np.percentile(seq_len_arr, 75)
    seq_len_min = np.min(seq_len_arr)
    seq_len_max = np.max(seq_len_arr)

    return [nr_seqs, seq_len_mean, seq_len_median, seq_len_q1, seq_len_q3, seq_len_min, seq_len_max]
     
################################################################################

def get_motif_similarites_ll(motif_ids_list, motif_pair2sim_dic,
                             min_max_norm=False):
    """
    Get motif similarities list of lists, in order of sorted motif IDs list.
    So e.g. 3 motifs: A, B, C, then list of lists with entry format A-A: similarity of A with A, ...
    [[A-A, A-B, A-C], [B-A, B-B, B-C], [C-A, C-B, C-C]]

    motif_ids_list: list of motif IDs.
    motif_pair2sim_dic format:
    motif_pair2sim_dic[motif_id1 + "," + motif_id2] = sim_score

    >>> motif_ids_list = ["A", "B"]
    >>> motif_pair2sim_dic = {"A,A": 1.0, "A,B": 0.6, "B,A": 0.6, "B,B": 1.0}
    >>> get_motif_similarites_ll(motif_ids_list, motif_pair2sim_dic)
    [[1.0, 0.6], [0.6, 1.0]]

    """
    assert motif_ids_list, "motif_ids_list empty"
    assert motif_pair2sim_dic, "motif_pair2sim_dic empty"

    motif_sim_ll = []

    # Sort motif IDs list.
    motif_ids_list.sort()
    for m1 in motif_ids_list:
        motif_sim_list = []
        for m2 in motif_ids_list:
            m1m2_id = m1 + "," + m2
            assert m1m2_id in motif_pair2sim_dic, "no similarity score stored for motif pair %s" %(m1m2_id)
            m1m2_sim = motif_pair2sim_dic[m1m2_id]
            motif_sim_list.append(m1m2_sim)
        motif_sim_ll.append(motif_sim_list)

    if min_max_norm:
        motif_sim_ll = min_max_normalize_list_of_lists(motif_sim_ll)

    return motif_sim_ll


################################################################################

def min_max_normalize_list_of_lists(in_ll):
    """
    Min-max normalize list of lists (i.e. 2d list of numeric values).

    >>> in_ll = [[0, 1, 2], [3, 4, 5], [6, 7, 10]]
    >>> min_max_normalize_list_of_lists(in_ll)
    [[0.0, 0.1, 0.2], [0.3, 0.4, 0.5], [0.6, 0.7, 1.0]]

    """
    assert in_ll, "input list of lists empty"

    min_val = 1000000
    max_val = -1000000
    for l in in_ll:
        for v in l:
            if v < min_val:
                min_val = v
            if v > max_val:
                max_val = v

    out_ll = []
    for l in in_ll:
        out_l = []
        for v in l:
            out_v = (v - min_val) / (max_val - min_val)
            out_l.append(out_v)
        out_ll.append(out_l)
        
    return out_ll


################################################################################

def enmo_generate_html_report(args,
                              motif_enrich_stats_dic,
                              seq_motif_blocks_dic,
                              benchlib_path,
                              df_pval=False, 
                              pval_cont_lll=False,
                              motif_pair2sim_dic=False,
                              pos_seqs_dic=False,
                              neg_seqs_dic=False,
                              pos_reg2annot_dic=False,
                              neg_reg2annot_dic=False,
                              annot2color_dic=False,
                              plotly_full_html=False,
                              plotly_embed_style=1,
                              rbpbench_mode="enmo",
                              html_report_out="report.rbpbench_enmo.html",
                              plots_subfolder="html_report_plots"):

    """
    Create motif enrichment statistics / plots (enmo mode).

    """

    # Use absolute paths?
    out_folder = args.out_folder
    if args.plot_abs_paths:
        out_folder = os.path.abspath(out_folder)

    # Version string.
    version_str = "v" + args.version

    plots_folder = plots_subfolder
    plots_out_folder = out_folder + "/" + plots_folder
    if args.plot_abs_paths:
        plots_folder = plots_out_folder

    # Delete folder if already present.
    if os.path.exists(plots_out_folder):
        shutil.rmtree(plots_out_folder)
    os.makedirs(plots_out_folder)

    html_out = out_folder + "/" + "report.rbpbench_enmo.html"
    md_out = out_folder + "/" + "report.rbpbench_enmo.md"
    if html_report_out:
        html_out = html_report_out

    # Number of input + background regions/sites.
    c_input_sites = args.c_input_sites
    c_bg_sites = args.c_bg_sites
    
    site_type_uc = "Genomic"
    site_type = "genomic"
    if not args.genomic_sites_input:
        site_type_uc = "Transcript"
        site_type = "transcript"

    regex_motif_info = ""
    # if args.regex:
    #     regex_motif_info = "Used regex motif: '%s'." %(args.regex)

    """
    Setup plotly .js to support plotly plots.

    """

    include_plotlyjs = "cdn"
    # plotly_full_html = False
    plotly_js_html = ""
    plotly_js_path = benchlib_path + "/content/plotly-2.20.0.min.js"
    assert os.path.exists(plotly_js_path), "plotly .js %s not found" %(plotly_js_path)
    if args.plotly_js_mode == 2:
        include_plotlyjs = plotly_js_path
    elif args.plotly_js_mode == 3:
        shutil.copy(plotly_js_path, plots_out_folder)
        include_plotlyjs = "plotly-2.20.0.min.js" # Or plots_folder + "/plotly-2.20.0.min.js" ?
    elif args.plotly_js_mode == 4:
        include_plotlyjs = True
        # plotly_full_html = False # Don't really need full html (head body ..) in plotly html.
    elif args.plotly_js_mode == 5:
        plotly_js_web = "https://cdn.plot.ly/plotly-2.25.2.min.js"
        plotly_js_html = '<script src="' + plotly_js_web + '"></script>' + "\n"
        include_plotlyjs = False
        # plotly_full_html = True
    elif args.plotly_js_mode == 6:
        shutil.copy(plotly_js_path, plots_out_folder)
        plotly_js = plots_folder + "/plotly-2.20.0.min.js"
        plotly_js_html = '<script src="' + plotly_js + '"></script>' + "\n"
        include_plotlyjs = False
    elif args.plotly_js_mode == 7:
        js_code = read_file_content_into_str_var(plotly_js_path)
        plotly_js_html = "<script>\n" + js_code + "\n</script>\n"
        include_plotlyjs = False
        # plotly_full_html = True

    """
    Setup sorttable.js to make tables in HTML sortable.

    """
    sorttable_js_path = benchlib_path + "/content/sorttable.js"
    assert os.path.exists(sorttable_js_path), "sorttable.js not at %s" %(sorttable_js_path)
    sorttable_js_html = '<script src="' + sorttable_js_path + '" type="text/javascript"></script>'
    if args.sort_js_mode == 2:
        shutil.copy(sorttable_js_path, plots_out_folder)
        sorttable_js_path = plots_folder + "/sorttable.js"
        sorttable_js_html = '<script src="' + sorttable_js_path + '" type="text/javascript"></script>'
    elif args.sort_js_mode == 3:
        js_code = read_file_content_into_str_var(sorttable_js_path)
        sorttable_js_html = "<script>\n" + js_code + "\n</script>\n"

    # Logo path.
    logo_path_html = plots_folder + "/logo.png"
    if not os.path.exists(logo_path_html):
        logo_path = benchlib_path + "/content/logo.png"
        assert os.path.exists(logo_path), "logo.png not found in %s" %(logo_path)
        shutil.copy(logo_path, plots_out_folder)

    # HTML head section.
    html_head = """<!DOCTYPE html>
<html>
<head>
<title>RBPBench - Motif Enrichment Report</title>
%s

<style>
    th, td {
        border: 1px solid black;
        padding: 8px;
        text-align: center;
    }
    th {
        background-color: #f2f2f2;
    }
    .page {
        page-break-after: always;
    }
    .title-container {
        display: flex;
        align-items: center;
    }
    .title-container img {
        margin-right: 10px; /* Adjust the spacing as needed */
    }
</style>

</head>

<div class="title-container">
    <img src="%s" alt="Logo" width="175">
    <h1>Motif Enrichment Report</h1>
</div>

<body>
""" %(plotly_js_html, logo_path_html)

    # HTML tail section.
    html_tail = """
%s
</body>
</html>
""" %(sorttable_js_html)

    # Markdown part.
    mdtext = """

List of available statistics and plots generated
by RBPBench (%s, rbpbench %s):

- [Motif enrichment statistics](#enmo-stats)
- [Motif co-occurrences heat map](#cooc-heat-map)
- [Sequence motif similarity vs significance PCA plot](#motif-sim-sig-plot)""" %(version_str, rbpbench_mode)

    mdtext += "\n"

    if pos_reg2annot_dic or neg_reg2annot_dic:
        mdtext += "- [Genomic region annotations](#reg-annot)\n"
    if pos_seqs_dic and neg_seqs_dic:
        mdtext += "- [k-mer distributions](#kmer-plotly)\n"

    add_head_info = ""
    if args.bed_sc_thr is not None:
        add_head_info = " BED score threshold (--bed-sc-thr) = %s" %(str(args.bed_sc_thr))
        if args.bed_sc_thr_rev_filter:
            add_head_info += " (reverse filtering applied, i.e., the lower the better)."
        else:
            add_head_info += "."

    mdtext += "\nFIMO p-value threshold (--fimo-pval) = %s.%s # of considered input regions = %i. Region extension (upstream, downstream) = (%i, %i).\n" %(str(args.fimo_pval), add_head_info, c_input_sites, args.ext_up, args.ext_down)
    mdtext += "\n&nbsp;\n"


    """
    Motif enrichment statistics.

    """

    motif_add_info = "Enriched"
    fisher_mode_info = "Fisher exact test alternative hypothesis is set to 'greater', i.e., significantly overrepresented motifs are reported."
    if args.fisher_mode == 2:
        fisher_mode_info = "Fisher exact test alternative hypothesis is set to 'two-sided', i.e., significantly over- and underrepresented motifs are reported."
        motif_add_info = "Enriched and depleted"
    elif args.fisher_mode == 3:
        fisher_mode_info = "Fisher exact test alternative hypothesis is set to 'less', i.e., significantly underrepresented motifs are reported."
        motif_add_info = "Depleted"

    p_val_info = "P-values below %s are considered significant." %(str(args.enmo_pval_thr))
    if args.enmo_pval_mode == 1:
        p_val_info = "%s motifs with Benjamini-Hochberg multiple testing corrected p-values below %s are considered significant." %(motif_add_info, str(args.enmo_pval_thr))
    elif args.enmo_pval_mode == 2:
        p_val_info = "%s motifs with p-values below %s (p-value threshold Bonferroni multiple testing corrected) are considered significant." %(motif_add_info, str(args.enmo_pval_thr))
    elif args.enmo_pval_mode == 3:
        p_val_info = "%s motifs with p-values below %s are considered significant." %(motif_add_info, str(args.enmo_pval_thr))
    else:
        assert False, "Invalid motif enrichment p-value mode (--enmo-pval-mode) set: %i" %(args.enmo_pval_mode)

    pval_dic = {}
    c_sig_motifs = 0
    sig_seq_motif_ids_list = []

    for motif_id in motif_enrich_stats_dic:
        pval = motif_enrich_stats_dic[motif_id].fisher_pval_corr
        pval_dic[motif_id] = pval
        if pval <= args.enmo_pval_thr:
            c_sig_motifs += 1
            if motif_enrich_stats_dic[motif_id].motif_type == "meme_xml":
                sig_seq_motif_ids_list.append(motif_id)

    sig_seq_motif_ids_list.sort()

    mdtext += """
## Motif enrichment statistics ### {#enmo-stats}

**Table:** RBP binding motif enrichment statistics. # of significant motifs = %i. Enrichment is calculated by comparing motif occurrences in input and background dataset. 
Based on the numbers of input and background sites with and without motif hits, 
Fisher's exact test is used to assess the significance of motif enrichment.
%s
%s
For full motif results list regardless of significance, see *motif_enrichment_stats.tsv* output table.
%s

""" %(c_sig_motifs, p_val_info, fisher_mode_info, regex_motif_info)

    mdtext += '<table style="max-width: 1200px; width: 100%; border-collapse: collapse; line-height: 0.8;">' + "\n"
    mdtext += "<thead>\n"
    mdtext += "<tr>\n"
    mdtext += "<th>RBP ID</th>\n"
    mdtext += "<th>Motif ID</th>\n"
    mdtext += "<th>Motif plot</th>\n"
    mdtext += "<th># in hits</th>\n"
    mdtext += "<th># bg hits</th>\n"
    mdtext += "<th># in</th>\n"
    mdtext += "<th># not in</th>\n"
    mdtext += "<th># bg</th>\n"
    mdtext += "<th># not bg</th>\n"
    mdtext += "<th>p-value</th>\n"
    mdtext += "</tr>\n"
    mdtext += "</thead>\n"
    mdtext += "<tbody>\n"

    for motif_id, fisher_pval in sorted(pval_dic.items(), key=lambda item: item[1], reverse=False):

        if fisher_pval > args.enmo_pval_thr:
            break

        rbp_id = motif_enrich_stats_dic[motif_id].rbp_id
        c_pos_hits = motif_enrich_stats_dic[motif_id].c_pos_hits
        c_neg_hits = motif_enrich_stats_dic[motif_id].c_neg_hits
        c_pos_hit_regions = motif_enrich_stats_dic[motif_id].c_pos_hit_regions
        c_neg_hit_regions = motif_enrich_stats_dic[motif_id].c_neg_hit_regions
        c_pos_regions = motif_enrich_stats_dic[motif_id].c_pos_regions
        c_neg_regions = motif_enrich_stats_dic[motif_id].c_neg_regions
        a_con = c_pos_hit_regions
        b_con = c_pos_regions - c_pos_hit_regions
        c_con = c_neg_hit_regions
        d_con = c_neg_regions - c_neg_hit_regions

        plot_str = "-"

        if motif_id in seq_motif_blocks_dic:

            motif_plot = "%s.%s.png" %(rbp_id, motif_id)
            motif_plot_out = plots_out_folder + "/" + motif_plot
            plot_path = plots_folder + "/" + motif_plot

            # Check if motif in motif plots folder.
            motif_path = benchlib_path + "/content/motif_plots/%s" %(motif_plot)
            if os.path.exists(motif_path):
                shutil.copy(motif_path, motif_plot_out)
                if args.plot_pdf:
                    create_motif_plot(motif_id, seq_motif_blocks_dic,
                                      motif_plot_out,
                                      plot_pdf=True,
                                      plot_png=False)

            if not os.path.exists(motif_plot_out):
                create_motif_plot(motif_id, seq_motif_blocks_dic,
                                  motif_plot_out,
                                  plot_pdf=args.plot_pdf,
                                  plot_png=True)

            plot_str = '<image src = "' + plot_path + '" width="300px"></image>'

        mdtext += '<tr>' + "\n"
        mdtext += "<td>" + rbp_id + "</td>\n"
        mdtext += "<td>" + motif_id + "</td>\n"
        mdtext += "<td>" + plot_str + "</td>\n"
        mdtext += "<td>" + str(c_pos_hits) + "</td>\n"
        mdtext += "<td>" + str(c_neg_hits) + "</td>\n"
        mdtext += "<td>" + str(a_con) + "</td>\n"
        mdtext += "<td>" + str(b_con) + "</td>\n"
        mdtext += "<td>" + str(c_con) + "</td>\n"
        mdtext += "<td>" + str(d_con) + "</td>\n"
        mdtext += "<td>" + str(fisher_pval) + "</td>\n"
        mdtext += '</tr>' + "\n"

    mdtext += '</tbody>' + "\n"
    mdtext += '</table>' + "\n"

    mdtext += "\n&nbsp;\n&nbsp;\n"
    mdtext += "\nColumn IDs have the following meanings: "
    mdtext += "**RBP ID** -> RBP ID belonging to motif ID, "
    mdtext += "**Motif ID** -> motif ID, "
    mdtext += "**Motif plot** -> visualization of motif (sequence logo for sequence motifs, otherwise -), "
    mdtext += '**# in hits** -> number of motif hits in input sites, '
    mdtext += '**# bg hits** -> number of motif hits in background sites, '
    mdtext += '**# in** -> number of input sites with motif hits, '
    mdtext += '**# not in** -> number of input sites without motif hits, '
    mdtext += '**# bg** -> number of background sites with motif hits, '
    mdtext += '**# not bg** -> number of background sites without motif hits, '
    mdtext += '**p-value** -> Fisher exact test p-value (corrected).' + "\n"
    mdtext += "\n&nbsp;\n"


    """
    Motif co-occurrences heat map.

    """

    mdtext += """
## Motif co-occurrences heat map ### {#cooc-heat-map}

"""

    if pval_cont_lll:

        cooc_plot_plotly =  "co-occurrence_plot.plotly.html"
        cooc_plot_plotly_out = plots_out_folder + "/" + cooc_plot_plotly

        create_cooc_plot_plotly(df_pval, pval_cont_lll, cooc_plot_plotly_out,
                                max_motif_dist=args.max_motif_dist,
                                min_motif_dist=args.min_motif_dist,
                                id1="Motif1",
                                id2="Motif2",
                                ids="Motifs",
                                include_plotlyjs=include_plotlyjs,
                                full_html=plotly_full_html)

        plot_path = plots_folder + "/" + cooc_plot_plotly

        if args.plotly_js_mode in [5, 6, 7]:
            # Read in plotly code.
            # mdtext += '<div style="width: 1200px; height: 1200px; align-items: center;">' + "\n"
            js_code = read_file_content_into_str_var(cooc_plot_plotly_out)
            js_code = js_code.replace("height:100%; width:100%;", "height:1200px; width:1200px;")
            mdtext += js_code + "\n"
            # mdtext += '</div>'
        else:
            if plotly_embed_style == 1:
                # mdtext += '<div class="container-fluid" style="margin-top:40px">' + "\n"
                mdtext += "<div>\n"
                mdtext += '<iframe src="' + plot_path + '" width="1200" height="1200"></iframe>' + "\n"
                mdtext += '</div>'
            elif plotly_embed_style == 2:
                mdtext += '<object data="' + plot_path + '" width="1200" height="1200"> </object>' + "\n"

        p_val_info = "P-values below %s are considered significant." %(str(args.cooc_pval_thr))

        min_motif_dist_info = ""
        if args.min_motif_dist > 0:
            min_motif_dist_info = " + a mean minimum motif distance >= %i nt " %(args.min_motif_dist)

        if args.cooc_pval_mode == 1:
            p_val_info = "Motif co-occurrences with Benjamini-Hochberg multiple testing corrected p-values below %s %sare considered significant." %(str(args.cooc_pval_thr), min_motif_dist_info)
        elif args.cooc_pval_mode == 2:
            p_val_info = "Motif co-occurrences with p-values below %s (p-value threshold Bonferroni multiple testing corrected) %sare considered significant." %(str(args.cooc_pval_thr), min_motif_dist_info)
        elif args.cooc_pval_mode == 3:
            p_val_info = "Motif co-occurrences with p-values below %s %sare considered significant." %(str(args.cooc_pval_thr), min_motif_dist_info)
        else:
            assert False, "Invalid co-occurrence p-value mode (--cooc-pval-mode) set: %i" %(args.cooc_pval_mode)
        
        p_val_info += " # of motif co-occurrence comparisons: %i. # of significant co-occurrences: %i (%.2f%%)." %(args.c_all_fisher_pval, args.c_sig_fisher_pval, args.perc_sig_fisher_pval)


        # Inform about set alterntive hypothesis for Fisher exact test on significant motif co-occurrences.
        motif_add_info = "enriched"
        fisher_mode_info = "Fisher exact test alternative hypothesis is set to 'greater', i.e., significantly overrepresented motif co-occurrences are reported."
        if args.fisher_mode == 2:
            fisher_mode_info = "Fisher exact test alternative hypothesis is set to 'two-sided', i.e., significantly over- and underrepresented motif co-occurrences are reported."
            motif_add_info = "enriched and depleted"
        elif args.fisher_mode == 3:
            fisher_mode_info = "Fisher exact test alternative hypothesis is set to 'less', i.e., significantly underrepresented motif co-occurrences are reported."
            motif_add_info = "depleted"

        mdtext += """

**Figure:** Heat map of co-occurrences (Fisher's exact test p-values) between motifs. 
Only significantly %s motifs (listed in upper table) are used in checking for siginificant co-occurrences.
Motif hit co-occurrences that are not significant are colored white, while
significant co-occurrences are colored according to their -log10 p-value (used as legend color, i.e., the higher the more significant).
%s
%s
Hover box: 
**1)** Motif1 in pair.
**2)** Motif2 in pair.
**3)** p-value: Fisher's exact test p-value (calculated based on contingency table (6) between Motif1 and Motif2). 
**4)** p-value after filtering: p-value after filtering, i.e., p-value is kept if significant (< %s), otherwise it is set to 1.0.
**5)** Motifs compaired.
**6)** Counts[]: contingency table of co-occurrence counts (i.e., number of %s regions with/without shared motif hits) between compaired motifs, 
with format [[A, B], [C, D]], where 
A: Motif1 AND Motif2, 
B: NOT Motif1 AND Motif2,
C: Motif1 AND NOT Motif2,
D: NOT Motif1 AND NOT Motif2.
**7)** Mean minimum distance of Motif1 and Motif2 hits (mean over all regions containing Motif1 + Motif2 motif hits). Distances measured from motif center positions.
**8)** Over all regions containing Motif1 and Motif2 pairs, percentage of regions where Motif1 + Motif2 motifs are within %i nt distance (set via --max-motif-dist).
**9)** Correlation: Pearson correlation coefficient between Motif1 and Motif2.
%s regions are labelled 1 or 0 (motif present or not), resulting in a vector of 1s and 0s for each motif.
Correlations are then calculated by comparing vectors for every pair of motifs.
**10)** -log10 of p-value after filtering, used for legend coloring. Using p-value after filtering, all non-significant p-values become 0 
for easier distinction between significant and non-significant co-occurrences.
%s

&nbsp;

""" %(motif_add_info, p_val_info, fisher_mode_info, str(args.cooc_pval_thr), site_type, args.max_motif_dist, site_type_uc, regex_motif_info)


    else:

        mdtext += """

No co-occurrences calculated as there are no significant motifs (see upper table).
        
&nbsp;

"""


    """
    Motif similarity (only for MEME formatted sequence motifs) vs significance PCA plot.
    
    """

    mdtext += """
## Sequence motif similarity vs significance PCA plot ### {#motif-sim-sig-plot}

"""

    # Get motif similarities of significant motifs.
    motif_sim_ll = False
    motif_sim_stats_dic = {}
    if motif_pair2sim_dic and len(sig_seq_motif_ids_list) > 2:  # at least 3 significant sequence motifs needed.
        motif_sim_ll = get_motif_similarites_ll(sig_seq_motif_ids_list, motif_pair2sim_dic,
                                                min_max_norm=args.motif_sim_norm)

        for motif_id in sig_seq_motif_ids_list:
            conseq = motif_enrich_stats_dic[motif_id].consensus_seq
            pval = motif_enrich_stats_dic[motif_id].fisher_pval_corr
            con_table_str = motif_enrich_stats_dic[motif_id].con_table
            pval = round_to_n_significant_digits_v2(pval, 4,
                                                    min_val=1e-304)
            log_pval = log_tf_pval(pval)
            log_pval = round_to_n_significant_digits_v2(log_pval, 4,
                                                        min_val=0)
            motif_sim_stats_dic[motif_id] = [conseq, con_table_str, pval, log_pval]

        # for idx, motif_id in enumerate(sig_seq_motif_ids_list):
        #     print(motif_id, motif_sim_stats_dic[motif_id][0], motif_sim_ll[idx])

        motif_sim_plot_plotly =  "motif_sim_sig_pca_plot.plotly.html"
        motif_sim_plot_plotly_out = plots_out_folder + "/" + motif_sim_plot_plotly

        create_pca_motif_sim_sig_plot_plotly(sig_seq_motif_ids_list, motif_sim_ll,
                                         motif_sim_stats_dic, motif_sim_plot_plotly_out,
                                         include_plotlyjs=include_plotlyjs,
                                         full_html=plotly_full_html)

        plot_path = plots_folder + "/" + motif_sim_plot_plotly

        if args.plotly_js_mode in [5, 6, 7]:
            js_code = read_file_content_into_str_var(motif_sim_plot_plotly_out)
            js_code = js_code.replace("height:100%; width:100%;", "height:900px; width:1000px;")
            mdtext += js_code + "\n"
        else:
            if plotly_embed_style == 1:
                mdtext += "<div>\n"
                mdtext += '<iframe src="' + plot_path + '" width="1000" height="900"></iframe>' + "\n"
                mdtext += '</div>'
            elif plotly_embed_style == 2:
                mdtext += '<object data="' + plot_path + '" width="1000" height="900"> </object>' + "\n"


        mdtext += """

**Figure:** Sequence motif similarity vs significance PCA plot. Motifs are arranged by their similarity and colored by their significance, i.e., their -log10 p-value (used as legend color, i.e., the higher the more significant).
Motifs closer together in the plot translates to higher motif similarity.
Motif similarity is measured using TOMTOM's euclidean distance measure between motif position weight matrices (PWMs), and motif similarity vectors are used for 2-dimensional PCA.
Only motifs that are sequence motifs and that are significantly %s (from upper table, sequence logos shown as well) are used for the comparison.
Hover box: 
**Motif** -> Motif ID.
**Consensus** -> Consensus sequence derived from PWM (PWM sequence logo can be found in upper motif enrichment statistics table).
**Counts** -> contingency table of motif occurrence counts (i.e., number of input and background regions with/without motif hits), 
with format [[A, B], [C, D]], where 
A: # input regions with motif hits, 
B: # input regions without motif hits,
C: # background regions with motif hits,
D: # background regions without motif hits.
**p-value** -> Fisher exact test p-value (corrected) derived from contingency table.
**-log10(p-value)** -> -log10 p-value of Fisher exact test p-value, used for coloring of motifs.

&nbsp;

""" %(motif_add_info)

    else:

        mdtext += """

No motif similarity vs significance plot generated since there are < 3 significant sequence motifs.
        
&nbsp;

"""



    """
    Genomic region annotations.

    """

    if pos_reg2annot_dic or neg_reg2annot_dic:

        mdtext += """
## Genomic region annotations ### {#reg-annot}

"""

    if pos_reg2annot_dic:

        annot_bar_plot =  "gene_region_annotation_bar_plot.input.png"
        annot_bar_plot_out = plots_out_folder + "/" + annot_bar_plot

        create_enmo_annotation_bar_plot(pos_reg2annot_dic, 
                                        annot2color_dic=annot2color_dic,
                                        data_id="",
                                        plot_pdf=args.plot_pdf,
                                        plot_out=annot_bar_plot_out)

        plot_path = plots_folder + "/" + annot_bar_plot

        mdtext += '<img src="' + plot_path + '" alt="Annotation bar plot input"' + "\n"
        mdtext += 'title="Annotation bar plot input" />' + "\n"
        mdtext += """
**Figure:** Gene region annotations for input sites (# sites = %i).

&nbsp;

""" %(c_input_sites)

    if neg_reg2annot_dic:

        annot_bar_plot =  "gene_region_annotation_bar_plot.background.png"
        annot_bar_plot_out = plots_out_folder + "/" + annot_bar_plot

        create_enmo_annotation_bar_plot(neg_reg2annot_dic,
                                        annot2color_dic=annot2color_dic,
                                        data_id="",
                                        plot_pdf=args.plot_pdf,
                                        plot_out=annot_bar_plot_out)

        plot_path = plots_folder + "/" + annot_bar_plot

        mdtext += '<img src="' + plot_path + '" alt="Annotation bar plot background"' + "\n"
        mdtext += 'title="Annotation bar plot background" />' + "\n"
        mdtext += """
**Figure:** Gene region annotations for background sites (# sites = %i).

&nbsp;

""" %(c_bg_sites)


    """
    k-mer distributions.
    
    """

    if pos_seqs_dic and neg_seqs_dic:


        # Get 3-mer percentages.
        pos_3mer_dic = seqs_dic_count_kmer_freqs(pos_seqs_dic, 3, rna=False,
                                                 return_ratios=True,
                                                 perc=True,
                                                 report_key_error=False,
                                                 skip_non_dic_keys=True,
                                                 convert_to_uc=True)
        neg_3mer_dic = seqs_dic_count_kmer_freqs(neg_seqs_dic, 3, rna=False,
                                                 return_ratios=True,
                                                 perc=True,
                                                 report_key_error=False,
                                                 skip_non_dic_keys=True,
                                                 convert_to_uc=True)
        # Get 4-mer percentages.
        pos_4mer_dic = seqs_dic_count_kmer_freqs(pos_seqs_dic, 4, rna=False,
                                                 return_ratios=True,
                                                 perc=True,
                                                 report_key_error=False,
                                                 skip_non_dic_keys=True,
                                                 convert_to_uc=True)
        neg_4mer_dic = seqs_dic_count_kmer_freqs(neg_seqs_dic, 4, rna=False,
                                                 return_ratios=True,
                                                 perc=True,
                                                 report_key_error=False,
                                                 skip_non_dic_keys=True,
                                                 convert_to_uc=True)
        # Get 5-mer percentages.
        pos_5mer_dic = seqs_dic_count_kmer_freqs(pos_seqs_dic, 5, rna=False,
                                                 return_ratios=True,
                                                 perc=True,
                                                 report_key_error=False,
                                                 skip_non_dic_keys=True,
                                                 convert_to_uc=True)
        neg_5mer_dic = seqs_dic_count_kmer_freqs(neg_seqs_dic, 5, rna=False,
                                                 return_ratios=True,
                                                 perc=True,
                                                 report_key_error=False,
                                                 skip_non_dic_keys=True,
                                                 convert_to_uc=True)

        mdtext += """
## k-mer distributions ### {#kmer-plotly}

Frequency distributions of k-mers (in percent) for the input and background dataset.

"""

        plotly_3mer_plot = "plotly_scatter_3mer.html"
        plotly_4mer_plot = "plotly_scatter_4mer.html"
        plotly_5mer_plot = "plotly_scatter_5mer.html"
        plotly_3mer_plot_out = plots_out_folder + "/" + plotly_3mer_plot
        plotly_4mer_plot_out = plots_out_folder + "/" + plotly_4mer_plot
        plotly_5mer_plot_out = plots_out_folder + "/" + plotly_5mer_plot

        # Create 3-mer plotly scatter plot.
        create_kmer_sc_plotly_scatter_plot(pos_3mer_dic, neg_3mer_dic, 3,
                                        plotly_3mer_plot_out,
                                        plotly_js_path)
        # Create 4-mer plotly scatter plot.
        create_kmer_sc_plotly_scatter_plot(pos_4mer_dic, neg_4mer_dic, 4,
                                        plotly_4mer_plot_out,
                                        plotly_js_path)
        # Create 5-mer plotly scatter plot.
        create_kmer_sc_plotly_scatter_plot(pos_5mer_dic, neg_5mer_dic, 5,
                                        plotly_5mer_plot_out,
                                        plotly_js_path)
        # Plot paths inside html report.
        plotly_3mer_plot_path = plots_folder + "/" + plotly_3mer_plot
        plotly_4mer_plot_path = plots_folder + "/" + plotly_4mer_plot
        plotly_5mer_plot_path = plots_folder + "/" + plotly_5mer_plot

        # R2 scores.
        r2_3mer = calc_r2_corr_measure(pos_3mer_dic, neg_3mer_dic,
                                    is_dic=True)
        r2_4mer = calc_r2_corr_measure(pos_4mer_dic, neg_4mer_dic,
                                    is_dic=True)
        r2_5mer = calc_r2_corr_measure(pos_5mer_dic, neg_5mer_dic,
                                    is_dic=True)

        if args.plotly_js_mode in [5, 6, 7]:
            js_code = read_file_content_into_str_var(plotly_3mer_plot_out)
            js_code = js_code.replace("height:100%; width:100%;", "height:800px; width:800px;")
            mdtext += js_code + "\n"
        else:
            mdtext += "<div>\n"
            mdtext += '<iframe src="' + plotly_3mer_plot_path + '" width="800" height="800"></iframe>' + "\n"
            mdtext += '</div>'

        # Old rplib style.
        # mdtext += '<div class=class="container-fluid" style="margin-top:40px">' + "\n"
        # mdtext += '<iframe src="' + plotly_3mer_plot_path + '" width="500" height="500"></iframe>' + "\n"
        # mdtext += '</div>'

        mdtext += """

**Figure:** 3-mer percentages in the input and background dataset. In case of
a uniform distribution with all 3-mers present, each 3-mer would have a
percentage = 1.5625. R2 = %.6f.

&nbsp;

""" %(r2_3mer)
    
        if args.plotly_js_mode in [5, 6, 7]:
            js_code = read_file_content_into_str_var(plotly_4mer_plot_out)
            js_code = js_code.replace("height:100%; width:100%;", "height:800px; width:800px;")
            mdtext += js_code + "\n"
        else:
            mdtext += "<div>\n"
            mdtext += '<iframe src="' + plotly_4mer_plot_path + '" width="800" height="800"></iframe>' + "\n"
            mdtext += '</div>'

        # mdtext += '<div class="container-fluid" style="margin-top:40px">' + "\n"
        # mdtext += '<iframe src="' + plotly_4mer_plot_path + '" width="600" height="600"></iframe>' + "\n"
        # mdtext += '</div>'

        mdtext += """

**Figure:** 4-mer percentages in the input and background dataset. In case of
a uniform distribution with all 4-mers present, each 4-mer would have a
percentage = 0.390625. R2 = %.6f.

&nbsp;

""" %(r2_4mer)

        if args.plotly_js_mode in [5, 6, 7]:
            js_code = read_file_content_into_str_var(plotly_5mer_plot_out)
            js_code = js_code.replace("height:100%; width:100%;", "height:800px; width:800px;")
            mdtext += js_code + "\n"
        else:
            mdtext += "<div>\n"
            mdtext += '<iframe src="' + plotly_5mer_plot_path + '" width="800" height="800"></iframe>' + "\n"
            mdtext += '</div>'

        # mdtext += '<div class="container-fluid" style="margin-top:40px">' + "\n"
        # mdtext += '<iframe src="' + plotly_5mer_plot_path + '" width="700" height="700"></iframe>' + "\n"
        # mdtext += '</div>'

        mdtext += """

**Figure:** 5-mer percentages in the input and background dataset. In case of
a uniform distribution with all 5-mers present, each 5-mer would have a
percentage = 0.09765625. R2 = %.6f.

&nbsp;

""" %(r2_5mer)

    # Convert mdtext to html.
    md2html = markdown(mdtext, extensions=['attr_list', 'tables'])

    # OUTMD = open(md_out,"w")
    # OUTMD.write("%s\n" %(mdtext))
    # OUTMD.close()

    html_content = html_head + md2html + html_tail

    OUTHTML = open(html_out,"w")
    OUTHTML.write("%s\n" %(html_content))
    OUTHTML.close()


################################################################################

def nemo_generate_html_report(args,
                              motif_enrich_stats_dic,
                              seq_motif_blocks_dic,
                              benchlib_path,
                              df_pval=False, 
                              pval_cont_lll=False,
                              motif_pair2sim_dic=False,
                              pos_seqs_dic=False,
                              neg_seqs_dic=False,
                              pos_reg2annot_dic=False,
                              neg_reg2annot_dic=False,
                              annot2color_dic=False,
                              plotly_full_html=False,
                              plotly_embed_style=1,
                              rbpbench_mode="nemo",
                              html_report_out="report.rbpbench_nemo.html",
                              plots_subfolder="html_report_plots"):

    """
    Create neighboring motif enrichment statistics / plots (nemo mode).

    """

    # Use absolute paths?
    out_folder = args.out_folder
    if args.plot_abs_paths:
        out_folder = os.path.abspath(out_folder)

    # Version string.
    version_str = "v" + args.version

    plots_folder = plots_subfolder
    plots_out_folder = out_folder + "/" + plots_folder
    if args.plot_abs_paths:
        plots_folder = plots_out_folder

    # Delete folder if already present.
    if os.path.exists(plots_out_folder):
        shutil.rmtree(plots_out_folder)
    os.makedirs(plots_out_folder)

    html_out = out_folder + "/" + "report.rbpbench_nemo.html"
    md_out = out_folder + "/" + "report.rbpbench_nemo.md"
    if html_report_out:
        html_out = html_report_out

    # Number of input + background regions/sites.
    c_input_sites = args.c_input_sites
    c_bg_sites = args.c_bg_sites
    
    site_type_uc = "Genomic"
    site_type = "genomic"
    if not args.genomic_sites_input:
        site_type_uc = "Transcript"
        site_type = "transcript"

    regex_motif_info = ""
    # if args.regex:
    #     regex_motif_info = "Used regex motif: '%s'." %(args.regex)

    """
    Setup plotly .js to support plotly plots.

    """

    include_plotlyjs = "cdn"
    # plotly_full_html = False
    plotly_js_html = ""
    plotly_js_path = benchlib_path + "/content/plotly-2.20.0.min.js"
    assert os.path.exists(plotly_js_path), "plotly .js %s not found" %(plotly_js_path)
    if args.plotly_js_mode == 2:
        include_plotlyjs = plotly_js_path
    elif args.plotly_js_mode == 3:
        shutil.copy(plotly_js_path, plots_out_folder)
        include_plotlyjs = "plotly-2.20.0.min.js" # Or plots_folder + "/plotly-2.20.0.min.js" ?
    elif args.plotly_js_mode == 4:
        include_plotlyjs = True
        # plotly_full_html = False # Don't really need full html (head body ..) in plotly html.
    elif args.plotly_js_mode == 5:
        plotly_js_web = "https://cdn.plot.ly/plotly-2.25.2.min.js"
        plotly_js_html = '<script src="' + plotly_js_web + '"></script>' + "\n"
        include_plotlyjs = False
        # plotly_full_html = True
    elif args.plotly_js_mode == 6:
        shutil.copy(plotly_js_path, plots_out_folder)
        plotly_js = plots_folder + "/plotly-2.20.0.min.js"
        plotly_js_html = '<script src="' + plotly_js + '"></script>' + "\n"
        include_plotlyjs = False
    elif args.plotly_js_mode == 7:
        js_code = read_file_content_into_str_var(plotly_js_path)
        plotly_js_html = "<script>\n" + js_code + "\n</script>\n"
        include_plotlyjs = False
        # plotly_full_html = True

    """
    Setup sorttable.js to make tables in HTML sortable.

    """
    sorttable_js_path = benchlib_path + "/content/sorttable.js"
    assert os.path.exists(sorttable_js_path), "sorttable.js not at %s" %(sorttable_js_path)
    sorttable_js_html = '<script src="' + sorttable_js_path + '" type="text/javascript"></script>'
    if args.sort_js_mode == 2:
        shutil.copy(sorttable_js_path, plots_out_folder)
        sorttable_js_path = plots_folder + "/sorttable.js"
        sorttable_js_html = '<script src="' + sorttable_js_path + '" type="text/javascript"></script>'
    elif args.sort_js_mode == 3:
        js_code = read_file_content_into_str_var(sorttable_js_path)
        sorttable_js_html = "<script>\n" + js_code + "\n</script>\n"

    # Logo path.
    logo_path_html = plots_folder + "/logo.png"
    if not os.path.exists(logo_path_html):
        logo_path = benchlib_path + "/content/logo.png"
        assert os.path.exists(logo_path), "logo.png not found in %s" %(logo_path)
        shutil.copy(logo_path, plots_out_folder)

    # HTML head section.
    html_head = """<!DOCTYPE html>
<html>
<head>
<title>RBPBench - Neighboring Motif Enrichment Report</title>
%s

<style>
    th, td {
        border: 1px solid black;
        padding: 8px;
        text-align: center;
    }
    th {
        background-color: #f2f2f2;
    }
    .page {
        page-break-after: always;
    }
    .title-container {
        display: flex;
        align-items: center;
    }
    .title-container img {
        margin-right: 10px; /* Adjust the spacing as needed */
    }
</style>

</head>

<div class="title-container">
    <img src="%s" alt="Logo" width="175">
    <h1>Neighboring Motif Enrichment Report</h1>
</div>

<body>
""" %(plotly_js_html, logo_path_html)

    # HTML tail section.
    html_tail = """
%s
</body>
</html>
""" %(sorttable_js_html)

    # Markdown part.
    mdtext = """

List of available statistics and plots generated
by RBPBench (%s, rbpbench %s):

- [Neighboring motif enrichment statistics](#nemo-stats)
- [Motif co-occurrences heat map](#cooc-heat-map)
- [Sequence motif similarity vs significance PCA plot](#motif-sim-sig-plot)
- [Sequence motif similarity vs direction PCA plot](#motif-sim-dir-plot)""" %(version_str, rbpbench_mode)

    mdtext += "\n"

    if pos_reg2annot_dic or neg_reg2annot_dic:
        mdtext += "- [Genomic region annotations](#reg-annot)\n"
    if pos_seqs_dic and neg_seqs_dic:
        mdtext += "- [k-mer distributions](#kmer-plotly)\n"

    add_head_info = ""
    if args.bed_sc_thr is not None:
        add_head_info = " BED score threshold (--bed-sc-thr) = %s" %(str(args.bed_sc_thr))
        if args.bed_sc_thr_rev_filter:
            add_head_info += " (reverse filtering applied, i.e., the lower the better)."
        else:
            add_head_info += "."

    mdtext += "\nFIMO p-value threshold (--fimo-pval) = %s.%s # of considered input regions = %i. Region extension (upstream, downstream) = (%i, %i).\n" %(str(args.fimo_pval), add_head_info, c_input_sites, args.ext_up, args.ext_down)
    mdtext += "\n&nbsp;\n"


    """
    Motif enrichment statistics.

    """

    motif_add_info = "Enriched"
    fisher_mode_info = "Fisher exact test alternative hypothesis is set to 'greater', i.e., significantly overrepresented motifs are reported."
    if args.fisher_mode == 2:
        fisher_mode_info = "Fisher exact test alternative hypothesis is set to 'two-sided', i.e., significantly over- and underrepresented motifs are reported."
        motif_add_info = "Enriched and depleted"
    elif args.fisher_mode == 3:
        fisher_mode_info = "Fisher exact test alternative hypothesis is set to 'less', i.e., significantly underrepresented motifs are reported."
        motif_add_info = "Depleted"

    p_val_info = "P-values (p-value column) below %s are considered significant." %(str(args.nemo_pval_thr))
    if args.nemo_pval_mode == 1:
        p_val_info = "%s motifs with Benjamini-Hochberg multiple testing corrected p-values (p-value column) below %s are considered significant." %(motif_add_info, str(args.nemo_pval_thr))
    elif args.nemo_pval_mode == 2:
        p_val_info = "%s motifs with p-values (p-value column) below %s (p-value threshold Bonferroni multiple testing corrected) are considered significant." %(motif_add_info, str(args.nemo_pval_thr))
    elif args.nemo_pval_mode == 3:
        p_val_info = "%s motifs with p-values (p-value column) below %s are considered significant." %(motif_add_info, str(args.nemo_pval_thr))
    else:
        assert False, "Invalid motif enrichment p-value mode (--nemo-pval-mode) set: %i" %(args.nemo_pval_mode)

    # Inform about set alterntive hypothesis for Wilcoxon rank-sum test.
    # wrs_mode_info = ""
    # if args.wrs_mode == 1:
    #     wrs_mode_info = "Wilcoxon rank sum test alternative hypothesis is set to 'two-sided', i.e., low WRS p-values (WRS p-value column) mean either up- or downstream context regions have significantly higher motif hit counts."
    # elif args.wrs_mode == 2:
    #     wrs_mode_info = "Wilcoxon rank sum test alternative hypothesis is set to 'greater', i.e., low WRS p-values (WRS p-value column) mean upstream context regions have significantly higher motif hit counts."
    # elif args.wrs_mode == 3:
    #     wrs_mode_info = "Wilcoxon rank sum test alternative hypothesis is set to 'less', i.e., low WRS p-values (WRS p-value column) mean downstream context regions have significantly higher motif hit counts."
    # else:
    #     assert False, "Invalid Wilcoxon rank sum test mode (--wrs-mode) set: %i" %(args.wrs_mode)
    wrs_mode_info = "Wilcoxon rank-sum test alternative hypothesis is set to 'two-sided', i.e., low WRS p-values (WRS p-value column) mean either up- or downstream context regions have significantly higher motif hit counts."

    ol_info = "Motif hits that overlap with the actual input sites are not counted."
    if args.allow_overlaps:
        ol_info = "Motif hits that overlap with the actual input sites are counted as well (--allow-overlaps enabled)."

    pval_dic = {}
    c_sig_motifs = 0
    sig_seq_motif_ids_list = []

    for motif_id in motif_enrich_stats_dic:
        pval = motif_enrich_stats_dic[motif_id].fisher_pval_corr
        pval_dic[motif_id] = pval
        if pval <= args.nemo_pval_thr:
            c_sig_motifs += 1
            if motif_enrich_stats_dic[motif_id].motif_type == "meme_xml":
                sig_seq_motif_ids_list.append(motif_id)

    sig_seq_motif_ids_list.sort()

    mdtext += """
## Neighboring motif enrichment statistics ### {#nemo-stats}

**Table:** Neighboring RBP binding motif enrichment statistics. # of significant motifs = %i. Enrichment is calculated by comparing motif occurrences in the context regions 
surrounding given input sites (up- and downstream context region size specified via --ext), effectively comparing the input with the background context regions.
%s
Based on the numbers of input and background context regions with and without motif hits, 
Fisher's exact test is used to assess the significance of motif enrichment.
%s
%s
For full motif results list regardless of significance, see *motif_enrichment_stats.tsv* output table.
To test whether up- or downstream regions have significantly higher motif hit counts, 
Wilcoxon rank-sum (WRS) test is applied.
%s
%s

""" %(c_sig_motifs, ol_info, p_val_info, fisher_mode_info, wrs_mode_info, regex_motif_info)

    mdtext += '<table style="max-width: 1400px; width: 100%; border-collapse: collapse; line-height: 0.9;">' + "\n"
    mdtext += "<thead>\n"
    mdtext += "<tr>\n"
    mdtext += "<th>RBP ID</th>\n"
    mdtext += "<th>Motif ID</th>\n"
    mdtext += "<th>Motif plot</th>\n"
    mdtext += "<th># in hits</th>\n"
    mdtext += "<th># bg hits</th>\n"
    mdtext += "<th># in</th>\n"
    mdtext += "<th># not in</th>\n"
    mdtext += "<th># bg</th>\n"
    mdtext += "<th># not bg</th>\n"
    mdtext += "<th>avg in dist</th>\n"
    # mdtext += "<th>max in dist</th>\n"
    # mdtext += "<th># max in dist</th>\n"
    mdtext += "<th>avg bg dist</th>\n"
    # mdtext += "<th>max bg dist</th>\n"
    # mdtext += "<th># max bg dist</th>\n"
    mdtext += "<th>WRS p-value</th>\n"
    mdtext += "<th>Motif distance plot</th>\n"
    mdtext += "<th>p-value</th>\n"
    mdtext += "</tr>\n"
    mdtext += "</thead>\n"
    mdtext += "<tbody>\n"

    for motif_id, fisher_pval in sorted(pval_dic.items(), key=lambda item: item[1], reverse=False):

        if fisher_pval > args.nemo_pval_thr:
            break

        rbp_id = motif_enrich_stats_dic[motif_id].rbp_id
        c_pos_hits = motif_enrich_stats_dic[motif_id].c_pos_hits
        c_neg_hits = motif_enrich_stats_dic[motif_id].c_neg_hits
        c_pos_hit_regions = motif_enrich_stats_dic[motif_id].c_pos_hit_regions
        c_neg_hit_regions = motif_enrich_stats_dic[motif_id].c_neg_hit_regions
        c_pos_regions = motif_enrich_stats_dic[motif_id].c_pos_regions
        c_neg_regions = motif_enrich_stats_dic[motif_id].c_neg_regions
        pos_avg_center_dist = round(motif_enrich_stats_dic[motif_id].pos_set_avg_center_dist, 1)
        neg_avg_center_dist = round(motif_enrich_stats_dic[motif_id].neg_set_avg_center_dist, 1)
        # pos_max_center_dist = motif_enrich_stats_dic[motif_id].pos_set_max_center_dist
        # pos_max_center_dist_c = motif_enrich_stats_dic[motif_id].pos_set_max_center_dist_c
        # neg_max_center_dist = motif_enrich_stats_dic[motif_id].neg_set_max_center_dist
        # neg_max_center_dist_c = motif_enrich_stats_dic[motif_id].neg_set_max_center_dist_c

        # dist_plot_counts_dic format: {pos: count, ...}, e.g. from -5 to 5: {5: 1, 4: 2, 3: 5, 2: 4, 1: 3, 0: 0, -1: 3, -2: 4, -3: 5, -4: 2, -5: 1}
        dist_plot_counts_dic = motif_enrich_stats_dic[motif_id].dist_plot_counts_dic
        wrs_pval = motif_enrich_stats_dic[motif_id].wrs_pval_two_sided

        a_con = c_pos_hit_regions
        b_con = c_pos_regions - c_pos_hit_regions
        c_con = c_neg_hit_regions
        d_con = c_neg_regions - c_neg_hit_regions

        motif_plot_str = "-"

        if motif_id in seq_motif_blocks_dic:

            motif_plot = "%s.%s.png" %(rbp_id, motif_id)
            motif_plot_out = plots_out_folder + "/" + motif_plot
            plot_path = plots_folder + "/" + motif_plot

            # Check if motif in motif plots folder.
            motif_path = benchlib_path + "/content/motif_plots/%s" %(motif_plot)
            if os.path.exists(motif_path):
                shutil.copy(motif_path, motif_plot_out)
                if args.plot_pdf:
                    create_motif_plot(motif_id, seq_motif_blocks_dic,
                                      motif_plot_out,
                                      plot_pdf=True,
                                      plot_png=False)

            if not os.path.exists(motif_plot_out):
                create_motif_plot(motif_id, seq_motif_blocks_dic,
                                  motif_plot_out,
                                  plot_pdf=args.plot_pdf,
                                  plot_png=True)

            motif_plot_str = '<image src = "' + plot_path + '" width="300px"></image>'

        # Plot motif distances to center as positions on x-axis and counts on y-axis.
        # positions = list(dist_plot_counts_dic.keys())
        # counts = list(dist_plot_counts_dic.values())
        positions = sorted(dist_plot_counts_dic.keys(), reverse=True)
        counts = [dist_plot_counts_dic[pos] for pos in positions]

        plt.figure(figsize=(16, 2))
        plt.plot(positions, counts, marker='o', linestyle='-', color='b', markersize=3)
        plt.xlabel('Position')
        plt.ylabel('Count')

        plt.xlim(min(positions) - 1, max(positions) + 1)

        plt.grid(False)
        plt.gca().spines['top'].set_visible(False)
        plt.gca().spines['right'].set_visible(False)
        plt.axvline(x=0, color='r', linestyle='--')

        plt.subplots_adjust(left=0.05, right=0.95, top=0.9, bottom=0.25)

        dist_plot = "%s.%s.motif_hit_dist_center.png" %(rbp_id, motif_id)
        dist_plot_out = plots_out_folder + "/" + dist_plot
        plot_path = plots_folder + "/" + dist_plot

        plt.savefig(dist_plot_out)

        if args.plot_pdf and dist_plot_out.endswith('.png'):
            pdf_out = dist_plot_out[:-4] + '.pdf'
            plt.savefig(pdf_out)

        plt.close()
        dist_plot_str = '<image src = "' + plot_path + '" width="300px"></image>'

        mdtext += '<tr>' + "\n"
        mdtext += "<td>" + rbp_id + "</td>\n"
        mdtext += "<td>" + motif_id + "</td>\n"
        mdtext += "<td>" + motif_plot_str + "</td>\n"
        mdtext += "<td>" + str(c_pos_hits) + "</td>\n"
        mdtext += "<td>" + str(c_neg_hits) + "</td>\n"
        mdtext += "<td>" + str(a_con) + "</td>\n"
        mdtext += "<td>" + str(b_con) + "</td>\n"
        mdtext += "<td>" + str(c_con) + "</td>\n"
        mdtext += "<td>" + str(d_con) + "</td>\n"
        mdtext += "<td>" + str(pos_avg_center_dist) + "</td>\n"
        # mdtext += "<td>" + str(pos_max_center_dist) + "</td>\n"
        # mdtext += "<td>" + str(pos_max_center_dist_c) + "</td>\n"
        mdtext += "<td>" + str(neg_avg_center_dist) + "</td>\n"
        # mdtext += "<td>" + str(neg_max_center_dist) + "</td>\n"
        # mdtext += "<td>" + str(neg_max_center_dist_c) + "</td>\n"
        mdtext += "<td>" + str(wrs_pval) + "</td>\n"
        mdtext += "<td>" + dist_plot_str + "</td>\n"
        mdtext += "<td>" + str(fisher_pval) + "</td>\n"
        mdtext += '</tr>' + "\n"

    mdtext += '</tbody>' + "\n"
    mdtext += '</table>' + "\n"

    mdtext += "\n&nbsp;\n&nbsp;\n"
    mdtext += "\nColumn IDs have the following meanings: "
    mdtext += "**RBP ID** -> RBP ID belonging to motif ID, "
    mdtext += "**Motif ID** -> motif ID, "
    mdtext += "**Motif plot** -> visualization of motif (sequence logo for sequence motifs, otherwise -), "
    mdtext += '**# in hits** -> number of motif hits in input sites, '
    mdtext += '**# bg hits** -> number of motif hits in background sites, '
    mdtext += '**# in** -> number of input sites with motif hits, '
    mdtext += '**# not in** -> number of input sites without motif hits, '
    mdtext += '**# bg** -> number of background sites with motif hits, '
    mdtext += '**# not bg** -> number of background sites without motif hits, '
    mdtext += '**avg in dist** -> average distance of motif hits to center of input sites (positive value indicates motifs tend to be located upstream of input sites, whereas negative value indicates downstream), '
    # mdtext += '**max in dist** -> distance position with maximum count (i.e., where most motif hit centers lie relative to input site centers), '
    # mdtext += '**# max in dist** -> number of motif hits at distance position with maximum count for input sites, '
    mdtext += '**avg bg dist** -> average distance of motif hits to center of background sites (positive value indicates motifs tend to be located upstream of input sites, whereas negative value indicates downstream), '
    # mdtext += '**max bg dist** -> distance position with maximum count (i.e., where most motif hit centers lie relative to background site centers), '
    # mdtext += '**# max bg dist** -> number of motif hits at distance position with maximum count for background sites, '
    mdtext += '**WRS p-value** -> Wilcoxon rank-sum test p-value to test for significantly different counts in up- and downstream context regions.' + "\n"
    mdtext += "**Motif distance plot** -> visualization of motif distance plot (counting motif hit center occurrences relative to input site centers), "
    mdtext += '**p-value** -> Fisher exact test p-value (corrected).' + "\n"
    mdtext += "\n&nbsp;\n"


    """
    Motif co-occurrences heat map.

    """

    mdtext += """
## Motif co-occurrences heat map ### {#cooc-heat-map}

"""

    if pval_cont_lll:

        cooc_plot_plotly =  "co-occurrence_plot.plotly.html"
        cooc_plot_plotly_out = plots_out_folder + "/" + cooc_plot_plotly

        create_cooc_plot_plotly(df_pval, pval_cont_lll, cooc_plot_plotly_out,
                                max_motif_dist=args.max_motif_dist,
                                min_motif_dist=args.min_motif_dist,
                                id1="Motif1",
                                id2="Motif2",
                                ids="Motifs",
                                include_plotlyjs=include_plotlyjs,
                                full_html=plotly_full_html)

        plot_path = plots_folder + "/" + cooc_plot_plotly

        if args.plotly_js_mode in [5, 6, 7]:
            # Read in plotly code.
            # mdtext += '<div style="width: 1200px; height: 1200px; align-items: center;">' + "\n"
            js_code = read_file_content_into_str_var(cooc_plot_plotly_out)
            js_code = js_code.replace("height:100%; width:100%;", "height:1200px; width:1200px;")
            mdtext += js_code + "\n"
            # mdtext += '</div>'
        else:
            if plotly_embed_style == 1:
                # mdtext += '<div class="container-fluid" style="margin-top:40px">' + "\n"
                mdtext += "<div>\n"
                mdtext += '<iframe src="' + plot_path + '" width="1200" height="1200"></iframe>' + "\n"
                mdtext += '</div>'
            elif plotly_embed_style == 2:
                mdtext += '<object data="' + plot_path + '" width="1200" height="1200"> </object>' + "\n"

        p_val_info = "P-values below %s are considered significant." %(str(args.cooc_pval_thr))

        min_motif_dist_info = ""
        if args.min_motif_dist > 0:
            min_motif_dist_info = " + a mean minimum motif distance >= %i nt " %(args.min_motif_dist)

        if args.cooc_pval_mode == 1:
            p_val_info = "Motif co-occurrences with Benjamini-Hochberg multiple testing corrected p-values below %s %sare considered significant." %(str(args.cooc_pval_thr), min_motif_dist_info)
        elif args.cooc_pval_mode == 2:
            p_val_info = "Motif co-occurrences with p-values below %s (p-value threshold Bonferroni multiple testing corrected) %sare considered significant." %(str(args.cooc_pval_thr), min_motif_dist_info)
        elif args.cooc_pval_mode == 3:
            p_val_info = "Motif co-occurrences with p-values below %s %sare considered significant." %(str(args.cooc_pval_thr), min_motif_dist_info)
        else:
            assert False, "Invalid co-occurrence p-value mode (--cooc-pval-mode) set: %i" %(args.cooc_pval_mode)
        
        p_val_info += " # of motif co-occurrence comparisons: %i. # of significant co-occurrences: %i (%.2f%%)." %(args.c_all_fisher_pval, args.c_sig_fisher_pval, args.perc_sig_fisher_pval)

        # Inform about set alterntive hypothesis for Fisher exact test on significant motif co-occurrences.
        motif_add_info = "enriched"
        fisher_mode_info = "Fisher exact test alternative hypothesis is set to 'greater', i.e., significantly overrepresented motif co-occurrences are reported."
        if args.fisher_mode == 2:
            fisher_mode_info = "Fisher exact test alternative hypothesis is set to 'two-sided', i.e., significantly over- and underrepresented motif co-occurrences are reported."
            motif_add_info = "enriched and depleted"
        elif args.fisher_mode == 3:
            fisher_mode_info = "Fisher exact test alternative hypothesis is set to 'less', i.e., significantly underrepresented motif co-occurrences are reported."
            motif_add_info = "depleted"

        mdtext += """

**Figure:** Heat map of co-occurrences (Fisher's exact test p-values) between motifs. 
Only significantly %s context region motifs (listed in upper table) are used in checking for siginificant co-occurrences.
Motif hit co-occurrences that are not significant are colored white, while
significant co-occurrences are colored according to their -log10 p-value (used as legend color, i.e., the higher the more significant).
%s
%s
Hover box: 
**1)** Motif1 in pair.
**2)** Motif2 in pair.
**3)** p-value: Fisher's exact test p-value (calculated based on contingency table (6) between Motif1 and Motif2). 
**4)** p-value after filtering: p-value after filtering, i.e., p-value is kept if significant (< %s), otherwise it is set to 1.0.
**5)** Motifs compaired.
**6)** Counts[]: contingency table of co-occurrence counts (i.e., number of %s regions with/without shared motif hits) between compaired motifs, 
with format [[A, B], [C, D]], where 
A: Motif1 AND Motif2, 
B: NOT Motif1 AND Motif2,
C: Motif1 AND NOT Motif2,
D: NOT Motif1 AND NOT Motif2.
**7)** Mean minimum distance of Motif1 and Motif2 hits (mean over all regions containing Motif1 + Motif2 motif hits). Distances measured from motif center positions.
**8)** Over all regions containing Motif1 and Motif2 pairs, percentage of regions where Motif1 + Motif2 motifs are within %i nt distance (set via --max-motif-dist).
**9)** Correlation: Pearson correlation coefficient between Motif1 and Motif2.
%s regions are labelled 1 or 0 (motif present or not), resulting in a vector of 1s and 0s for each motif.
Correlations are then calculated by comparing vectors for every pair of motifs.
**10)** -log10 of p-value after filtering, used for legend coloring. Using p-value after filtering, all non-significant p-values become 0 
for easier distinction between significant and non-significant co-occurrences.
%s

&nbsp;

""" %(motif_add_info, p_val_info, fisher_mode_info, str(args.cooc_pval_thr), site_type, args.max_motif_dist, site_type_uc, regex_motif_info)


    else:

        mdtext += """

No co-occurrences calculated as no significant context region motifs were found (see upper table).
        
&nbsp;

"""



    """
    Motif similarity (only for MEME formatted sequence motifs) vs significance PCA plot.

    """

    mdtext += """
## Sequence motif similarity vs significance PCA plot ### {#motif-sim-sig-plot}

"""

    # Get motif similarities of significant motifs.
    motif_sim_ll = False
    motif_sim_stats_dic = {}
    if motif_pair2sim_dic and len(sig_seq_motif_ids_list) > 2:  # at least 3 significant sequence motifs needed.
        motif_sim_ll = get_motif_similarites_ll(sig_seq_motif_ids_list, motif_pair2sim_dic,
                                                min_max_norm=args.motif_sim_norm)

        for motif_id in sig_seq_motif_ids_list:
            conseq = motif_enrich_stats_dic[motif_id].consensus_seq
            pval = motif_enrich_stats_dic[motif_id].fisher_pval_corr
            con_table_str = motif_enrich_stats_dic[motif_id].con_table
            pval = round_to_n_significant_digits_v2(pval, 4,
                                                    min_val=1e-304)
            log_pval = log_tf_pval(pval)
            log_pval = round_to_n_significant_digits_v2(log_pval, 4,
                                                        min_val=0)
            motif_sim_stats_dic[motif_id] = [conseq, con_table_str, pval, log_pval]

        # for idx, motif_id in enumerate(sig_seq_motif_ids_list):
        #     print(motif_id, motif_sim_stats_dic[motif_id][0], motif_sim_ll[idx])

        motif_sim_plot_plotly =  "motif_sim_sig_pca_plot.plotly.html"
        motif_sim_plot_plotly_out = plots_out_folder + "/" + motif_sim_plot_plotly

        create_pca_motif_sim_sig_plot_plotly(sig_seq_motif_ids_list, motif_sim_ll,
                                         motif_sim_stats_dic, motif_sim_plot_plotly_out,
                                         include_plotlyjs=include_plotlyjs,
                                         full_html=plotly_full_html)

        plot_path = plots_folder + "/" + motif_sim_plot_plotly

        if args.plotly_js_mode in [5, 6, 7]:
            js_code = read_file_content_into_str_var(motif_sim_plot_plotly_out)
            js_code = js_code.replace("height:100%; width:100%;", "height:900px; width:1000px;")
            mdtext += js_code + "\n"
        else:
            if plotly_embed_style == 1:
                mdtext += "<div>\n"
                mdtext += '<iframe src="' + plot_path + '" width="1000" height="900"></iframe>' + "\n"
                mdtext += '</div>'
            elif plotly_embed_style == 2:
                mdtext += '<object data="' + plot_path + '" width="1000" height="900"> </object>' + "\n"


        mdtext += """

**Figure:** Sequence motif similarity vs significance PCA plot. Motifs are arranged by their similarity and colored by their significance, i.e., their -log10 p-value (used as legend color, i.e., the higher the more significant).
Motifs closer together in the plot translates to higher motif similarity.
Motif similarity is measured using TOMTOM's euclidean distance measure between motif position weight matrices (PWMs), and motif similarity vectors are used for 2-dimensional PCA.
Only motifs that are sequence motifs and that are significantly %s (from upper table, sequence logos shown as well) are used for the comparison.
Hover box: 
**Motif** -> Motif ID.
**Consensus** -> Consensus sequence derived from PWM (PWM sequence logo can be found in upper motif enrichment statistics table).
**Counts** -> contingency table of motif occurrence counts (i.e., number of input and background regions with/without motif hits), 
with format [[A, B], [C, D]], where 
A: # input regions with motif hits, 
B: # input regions without motif hits,
C: # background regions with motif hits,
D: # background regions without motif hits.
**p-value** -> Fisher exact test p-value (corrected) derived from contingency table.
**-log10(p-value)** -> -log10 p-value of Fisher exact test p-value, used for coloring of motifs.

&nbsp;

""" %(motif_add_info)

    else:

        mdtext += """

No motif similarity vs significance plot generated since there are < 3 significant sequence motifs.
        
&nbsp;

"""






    """
    Motif similarity (only for MEME formatted sequence motifs) vs direction / context PCA plot.

    """

    mdtext += """
## Sequence motif similarities vs direction PCA plot ### {#motif-sim-dir-plot}

"""

    # Get motif similarities of significant motifs.
    motif_sim_ll = False
    motif_sim_stats_dic = {}
    if motif_pair2sim_dic and len(sig_seq_motif_ids_list) > 2:  # at least 3 significant sequence motifs needed.
        motif_sim_ll = get_motif_similarites_ll(sig_seq_motif_ids_list, motif_pair2sim_dic,
                                                min_max_norm=args.motif_sim_norm)

        for motif_id in sig_seq_motif_ids_list:
            conseq = motif_enrich_stats_dic[motif_id].consensus_seq
            pval = motif_enrich_stats_dic[motif_id].fisher_pval_corr
            con_table_str = motif_enrich_stats_dic[motif_id].con_table
            wrs_pval_greater = motif_enrich_stats_dic[motif_id].wrs_pval_greater
            wrs_pval_less = motif_enrich_stats_dic[motif_id].wrs_pval_less

            pval = round_to_n_significant_digits_v2(pval, 4,
                                                    min_val=1e-304)
            log_pval = log_tf_pval(pval)
            log_pval = round_to_n_significant_digits_v2(log_pval, 4,
                                                        min_val=0)

            wrs_pval = 1.0
            sign = 1
            if wrs_pval_greater < wrs_pval_less: # upstream p-value lower.
                wrs_pval = wrs_pval_greater
                sign = -1
            elif wrs_pval_less < wrs_pval_greater:
                wrs_pval = wrs_pval_less
            
            wrs_pval = round_to_n_significant_digits_v2(wrs_pval, 4,
                                                    min_val=1e-304)
            log_wrs_pval = log_tf_pval(wrs_pval)
            log_wrs_pval = round_to_n_significant_digits_v2(log_wrs_pval, 4,
                                                            min_val=0)
            
            # log_wrs_pval zero if wrs_pval_greater == wrs_pval_less, negative if wrs_pval_greater greater, positive if wrs_pval_less greater.
            log_wrs_pval = sign * log_wrs_pval

            motif_sim_stats_dic[motif_id] = [conseq, con_table_str, pval, log_pval, wrs_pval_greater, wrs_pval_less, log_wrs_pval]


        motif_sim_plot_plotly =  "motif_sim_dir_pca_plot.plotly.html"
        motif_sim_plot_plotly_out = plots_out_folder + "/" + motif_sim_plot_plotly

        create_pca_motif_sim_dir_plot_plotly(sig_seq_motif_ids_list, motif_sim_ll,
                                         motif_sim_stats_dic, motif_sim_plot_plotly_out,
                                         include_plotlyjs=include_plotlyjs,
                                         full_html=plotly_full_html)

        plot_path = plots_folder + "/" + motif_sim_plot_plotly

        if args.plotly_js_mode in [5, 6, 7]:
            js_code = read_file_content_into_str_var(motif_sim_plot_plotly_out)
            js_code = js_code.replace("height:100%; width:100%;", "height:900px; width:1000px;")
            mdtext += js_code + "\n"
        else:
            if plotly_embed_style == 1:
                mdtext += "<div>\n"
                mdtext += '<iframe src="' + plot_path + '" width="1000" height="900"></iframe>' + "\n"
                mdtext += '</div>'
            elif plotly_embed_style == 2:
                mdtext += '<object data="' + plot_path + '" width="1000" height="900"> </object>' + "\n"


        mdtext += """

**Figure:** Sequence motif similarity vs direction PCA plot. Motifs are arranged by their similarity and colored by their preferred binding direction relative to the provided central sites (up- or downstream context). 
Significance of directional preference is given as the -log10 of the Wilcoxon rank sum-test p-value (used as legend color, i.e., the more negative or positive the value the more significant the preference).
A negative -log10 p-value indicates a preference for upstream context, whereas a positive value indicates a preference for downstream context.
A -log10 p-value close to 0 indicates no significant directional preference.
Motifs closer together in the plot translates to higher motif similarity.
Motif similarity is measured using TOMTOM's euclidean distance measure between motif position weight matrices (PWMs), and motif similarity vectors are used for 2-dimensional PCA.
Only motifs that are sequence motifs and that are significantly %s (from upper table, also containing their sequence logos) are used for the comparison.
Hover box: 
**Motif** -> Motif ID.
**Consensus** -> Consensus sequence derived from PWM (PWM sequence logo can be found in upper motif enrichment statistics table).
**Counts** -> contingency table of motif occurrence counts (i.e., number of input and background regions with/without motif hits), 
with format [[A, B], [C, D]], where 
A: # input regions with motif hits, 
B: # input regions without motif hits,
C: # background regions with motif hits,
D: # background regions without motif hits.
**p-value** -> Fisher exact test p-value (corrected) derived from contingency table.
**WRS p-value (upstream)** -> Wilcoxon rank-sum test p-value to test for significantly higher counts in upstream context regions.
**WRS p-value (downstream)** -> Wilcoxon rank-sum test p-value to test for significantly higher counts in downstream context regions.
**-log10(WRS p-value)** -> -log10 p-value of Wilcoxon rank-sum test p-value used for coloring of motifs. 
The smaller of the two WRS p-values is taken. A negative value indicates upstream preference, a positive value downstream preference.

&nbsp;

""" %(motif_add_info)

    else:

        mdtext += """

No motif similarity vs direction plot generated since there are < 3 significant sequence motifs.
        
&nbsp;

"""


    """
    Genomic region annotations.

    """

    if pos_reg2annot_dic or neg_reg2annot_dic:

        mdtext += """
## Genomic region annotations ### {#reg-annot}

"""

    if pos_reg2annot_dic:

        annot_bar_plot =  "gene_region_annotation_bar_plot.input.png"
        annot_bar_plot_out = plots_out_folder + "/" + annot_bar_plot

        create_enmo_annotation_bar_plot(pos_reg2annot_dic, 
                                        annot2color_dic=annot2color_dic,
                                        data_id="",
                                        plot_pdf=args.plot_pdf,
                                        plot_out=annot_bar_plot_out)

        plot_path = plots_folder + "/" + annot_bar_plot

        mdtext += '<img src="' + plot_path + '" alt="Annotation bar plot input"' + "\n"
        mdtext += 'title="Annotation bar plot input" />' + "\n"
        mdtext += """
**Figure:** Gene region annotations for input sites (# sites = %i).

&nbsp;

""" %(c_input_sites)

    if neg_reg2annot_dic:

        annot_bar_plot =  "gene_region_annotation_bar_plot.background.png"
        annot_bar_plot_out = plots_out_folder + "/" + annot_bar_plot

        create_enmo_annotation_bar_plot(neg_reg2annot_dic,
                                        annot2color_dic=annot2color_dic,
                                        data_id="",
                                        plot_pdf=args.plot_pdf,
                                        plot_out=annot_bar_plot_out)

        plot_path = plots_folder + "/" + annot_bar_plot

        mdtext += '<img src="' + plot_path + '" alt="Annotation bar plot background"' + "\n"
        mdtext += 'title="Annotation bar plot background" />' + "\n"
        mdtext += """
**Figure:** Gene region annotations for background sites (# sites = %i).

&nbsp;

""" %(c_bg_sites)




    """
    k-mer distributions.
    
    """

    if pos_seqs_dic and neg_seqs_dic:


        # Get 3-mer percentages.
        pos_3mer_dic = seqs_dic_count_kmer_freqs(pos_seqs_dic, 3, rna=False,
                                                 return_ratios=True,
                                                 perc=True,
                                                 report_key_error=False,
                                                 skip_non_dic_keys=True,
                                                 convert_to_uc=True)
        neg_3mer_dic = seqs_dic_count_kmer_freqs(neg_seqs_dic, 3, rna=False,
                                                 return_ratios=True,
                                                 perc=True,
                                                 report_key_error=False,
                                                 skip_non_dic_keys=True,
                                                 convert_to_uc=True)
        # Get 4-mer percentages.
        pos_4mer_dic = seqs_dic_count_kmer_freqs(pos_seqs_dic, 4, rna=False,
                                                 return_ratios=True,
                                                 perc=True,
                                                 report_key_error=False,
                                                 skip_non_dic_keys=True,
                                                 convert_to_uc=True)
        neg_4mer_dic = seqs_dic_count_kmer_freqs(neg_seqs_dic, 4, rna=False,
                                                 return_ratios=True,
                                                 perc=True,
                                                 report_key_error=False,
                                                 skip_non_dic_keys=True,
                                                 convert_to_uc=True)
        # Get 5-mer percentages.
        pos_5mer_dic = seqs_dic_count_kmer_freqs(pos_seqs_dic, 5, rna=False,
                                                 return_ratios=True,
                                                 perc=True,
                                                 report_key_error=False,
                                                 skip_non_dic_keys=True,
                                                 convert_to_uc=True)
        neg_5mer_dic = seqs_dic_count_kmer_freqs(neg_seqs_dic, 5, rna=False,
                                                 return_ratios=True,
                                                 perc=True,
                                                 report_key_error=False,
                                                 skip_non_dic_keys=True,
                                                 convert_to_uc=True)

        mdtext += """
## k-mer distributions ### {#kmer-plotly}

Frequency distributions of k-mers (in percent) for the input and background dataset.

"""

        plotly_3mer_plot = "plotly_scatter_3mer.html"
        plotly_4mer_plot = "plotly_scatter_4mer.html"
        plotly_5mer_plot = "plotly_scatter_5mer.html"
        plotly_3mer_plot_out = plots_out_folder + "/" + plotly_3mer_plot
        plotly_4mer_plot_out = plots_out_folder + "/" + plotly_4mer_plot
        plotly_5mer_plot_out = plots_out_folder + "/" + plotly_5mer_plot

        # Create 3-mer plotly scatter plot.
        create_kmer_sc_plotly_scatter_plot(pos_3mer_dic, neg_3mer_dic, 3,
                                        plotly_3mer_plot_out,
                                        plotly_js_path)
        # Create 4-mer plotly scatter plot.
        create_kmer_sc_plotly_scatter_plot(pos_4mer_dic, neg_4mer_dic, 4,
                                        plotly_4mer_plot_out,
                                        plotly_js_path)
        # Create 5-mer plotly scatter plot.
        create_kmer_sc_plotly_scatter_plot(pos_5mer_dic, neg_5mer_dic, 5,
                                        plotly_5mer_plot_out,
                                        plotly_js_path)
        # Plot paths inside html report.
        plotly_3mer_plot_path = plots_folder + "/" + plotly_3mer_plot
        plotly_4mer_plot_path = plots_folder + "/" + plotly_4mer_plot
        plotly_5mer_plot_path = plots_folder + "/" + plotly_5mer_plot

        # R2 scores.
        r2_3mer = calc_r2_corr_measure(pos_3mer_dic, neg_3mer_dic,
                                    is_dic=True)
        r2_4mer = calc_r2_corr_measure(pos_4mer_dic, neg_4mer_dic,
                                    is_dic=True)
        r2_5mer = calc_r2_corr_measure(pos_5mer_dic, neg_5mer_dic,
                                    is_dic=True)

        if args.plotly_js_mode in [5, 6, 7]:
            js_code = read_file_content_into_str_var(plotly_3mer_plot_out)
            js_code = js_code.replace("height:100%; width:100%;", "height:800px; width:800px;")
            mdtext += js_code + "\n"
        else:
            mdtext += "<div>\n"
            mdtext += '<iframe src="' + plotly_3mer_plot_path + '" width="800" height="800"></iframe>' + "\n"
            mdtext += '</div>'

        # Old rplib style.
        # mdtext += '<div class=class="container-fluid" style="margin-top:40px">' + "\n"
        # mdtext += '<iframe src="' + plotly_3mer_plot_path + '" width="500" height="500"></iframe>' + "\n"
        # mdtext += '</div>'

        mdtext += """

**Figure:** 3-mer percentages in the input and background dataset. In case of
a uniform distribution with all 3-mers present, each 3-mer would have a
percentage = 1.5625. R2 = %.6f.

&nbsp;

""" %(r2_3mer)
    
        if args.plotly_js_mode in [5, 6, 7]:
            js_code = read_file_content_into_str_var(plotly_4mer_plot_out)
            js_code = js_code.replace("height:100%; width:100%;", "height:800px; width:800px;")
            mdtext += js_code + "\n"
        else:
            mdtext += "<div>\n"
            mdtext += '<iframe src="' + plotly_4mer_plot_path + '" width="800" height="800"></iframe>' + "\n"
            mdtext += '</div>'

        # mdtext += '<div class="container-fluid" style="margin-top:40px">' + "\n"
        # mdtext += '<iframe src="' + plotly_4mer_plot_path + '" width="600" height="600"></iframe>' + "\n"
        # mdtext += '</div>'

        mdtext += """

**Figure:** 4-mer percentages in the input and background dataset. In case of
a uniform distribution with all 4-mers present, each 4-mer would have a
percentage = 0.390625. R2 = %.6f.

&nbsp;

""" %(r2_4mer)

        if args.plotly_js_mode in [5, 6, 7]:
            js_code = read_file_content_into_str_var(plotly_5mer_plot_out)
            js_code = js_code.replace("height:100%; width:100%;", "height:800px; width:800px;")
            mdtext += js_code + "\n"
        else:
            mdtext += "<div>\n"
            mdtext += '<iframe src="' + plotly_5mer_plot_path + '" width="800" height="800"></iframe>' + "\n"
            mdtext += '</div>'

        # mdtext += '<div class="container-fluid" style="margin-top:40px">' + "\n"
        # mdtext += '<iframe src="' + plotly_5mer_plot_path + '" width="700" height="700"></iframe>' + "\n"
        # mdtext += '</div>'

        mdtext += """

**Figure:** 5-mer percentages in the input and background dataset. In case of
a uniform distribution with all 5-mers present, each 5-mer would have a
percentage = 0.09765625. R2 = %.6f.

&nbsp;

""" %(r2_5mer)

    # Convert mdtext to html.
    md2html = markdown(mdtext, extensions=['attr_list', 'tables'])

    # OUTMD = open(md_out,"w")
    # OUTMD.write("%s\n" %(mdtext))
    # OUTMD.close()

    html_content = html_head + md2html + html_tail

    OUTHTML = open(html_out,"w")
    OUTHTML.write("%s\n" %(html_content))
    OUTHTML.close()


################################################################################

def create_kmer_sc_plotly_scatter_plot(pos_mer_dic, neg_mer_dic, k,
                                       out_html, plotly_js,
                                       pos_label="k-mer % input",
                                       neg_label="k-mer % background",
                                       kmer_label="k-mer"):
    """
    Create plotly graph plot, containing k-mer scores of positive
    and negative set, and store in .html file.

    pos_mer_dic:
        dic with k-mer percentages of positive set.
    neg_mer_dic:
        dic with k-mer percentages of negative set.
    k:
        k in k-mer.
    out_html:
        Output .html path to store interactive (!) plotly graph.
    plotly_js:
        Path to plotly js plotly-latest.min.js.

    """
    assert pos_mer_dic, "given pos_mer_dic empty"
    assert neg_mer_dic, "given neg_mer_dic empty"
    assert len(pos_mer_dic) == len(neg_mer_dic), "len(pos_mer_dic) != len(neg_mer_dic)"

    data = {pos_label : [], neg_label : [], kmer_label : []}

    max_pos_perc = 0
    max_neg_perc = 0

    for kmer in pos_mer_dic:
        pos_perc = pos_mer_dic[kmer]
        neg_perc = neg_mer_dic[kmer]
        if pos_perc != 0:
            pos_perc = round(pos_perc, k)
        if neg_perc != 0:
            neg_perc = round(neg_perc, k)
        if pos_perc > max_pos_perc:
            max_pos_perc = pos_perc
        if neg_perc > max_neg_perc:
            max_neg_perc = neg_perc
        data[pos_label].append(pos_perc)
        data[neg_label].append(neg_perc)
        data[kmer_label].append(kmer)

    # Get min and max axis values for scaling.
    min_perc = 0
    max_perc = max_pos_perc
    if max_neg_perc > max_pos_perc:
        max_perc = max_neg_perc

    # Find out how to round up max_perc.
    if re.search(r"\d+\.\d+", str(max_perc)):
        m = re.search(r"(\d+)\.(\d+)", str(max_perc))
        left = str(m.group(1))
        right = str(m.group(2))
    else:
        assert False, "no pattern match on max_perc"
    if left == "0":
        for i,c in enumerate(right):
            prec = i + 1
            if c != "0":
                # Custom decimal round up.
                max_perc = decimal_ceil(max_perc, prec)
                break
    else:
        # Round up to whole number with math.ceil.
        max_perc = ceil(max_perc)

    df = pd.DataFrame(data, columns = [pos_label, neg_label, kmer_label])

    # # 3utr blue: #1f77b4
    # dot_col = "#1f77b4"  
    # # Color of dots.
    # dot_col = "#69e9f6"
    # if theme == 2:
    #     dot_col = "blue"

    # plot = px.scatter(data_frame=df, x=pos_label, y=neg_label, hover_name=kmer_label,
    #                   color_discrete_sequence=[dot_col])
    # plot.layout.template = 'seaborn'

    dot_col = "#2b7bba"

    plot = px.scatter(data_frame=df, x=pos_label, y=neg_label, hover_name=kmer_label, color_discrete_sequence=[dot_col])

    plot.update_layout(yaxis_range=[min_perc, max_perc])
    plot.update_layout(xaxis_range=[min_perc, max_perc])

    plot.write_html(out_html,
                    full_html=False,
                    include_plotlyjs=plotly_js)


################################################################################

def decimal_ceil(a, prec):
    """
    Round up a given decimal number at a certain precision.

    >>> a = 0.002489
    >>> decimal_ceil(a, 3)
    0.003
    >>> decimal_ceil(a, 2)
    0.01
    """
    a_rounded = np.round(a + 0.5 * 10**(-prec), prec)
    return float(a_rounded)


################################################################################

def calc_r2_corr_measure(scores1, scores2,
                         is_dic=False):
    """
    Calculate R2 measure.

    is_dic:
        If scores1 + scores2 are dictionaries.

    """
    assert len(scores1) == len(scores2), "len(scores1) != len(scores2)"

    if is_dic:
        sc1 = []
        sc2 = []
        for dic_key in scores1:
            sc1.append(scores1[dic_key])
            sc2.append(scores2[dic_key])
        correlation_matrix = np.corrcoef(sc1, sc2)
    else:
        correlation_matrix = np.corrcoef(scores1, scores2)
    correlation_xy = correlation_matrix[0,1]
    return correlation_xy**2


################################################################################

def split_regions_by_sc(reg2sc_dic, top_n=False, bottom_n=False,
                        rev_sort=True):
    """
    Split region -> score dictionary by score. Return top_n and bottom_n
    regions.

    If top_n and bottom_n not set, split regions sorted by score 
    at the midpoint.
    If top_n and not bottom_n, split regions at the top_n.
    If bottom_n and not top_n, split regions at the bottom_n.
    If top_n and bottom_n, return the top_n and bottom_n regions.
    If any of top_n, bottom_n or top_n + bottom_n > len(reg2sc_dic),
    split regions at midpoint.

    rev_sort: 
        If True, the regions are sorted in descending order of scores, 
        so higher scores means better regions.

    >>> reg2sc_dic = {}
    >>> reg2sc_dic["reg1"] = 0.1
    >>> reg2sc_dic["reg2"] = 0.2
    >>> reg2sc_dic["reg3"] = 0.3
    >>> reg2sc_dic["reg4"] = 0.4
    >>> reg2sc_dic["reg5"] = 0.5
    >>> reg2sc_dic["reg6"] = 0.6
    >>> reg2sc_dic["reg7"] = 0.7
    >>> reg2sc_dic["reg8"] = 0.8
    >>> split_regions_by_sc(reg2sc_dic)
    (['reg8', 'reg7', 'reg6', 'reg5'], ['reg4', 'reg3', 'reg2', 'reg1'])
    >>> split_regions_by_sc(reg2sc_dic, top_n=2)
    (['reg8', 'reg7'], ['reg6', 'reg5', 'reg4', 'reg3', 'reg2', 'reg1'])
    >>> split_regions_by_sc(reg2sc_dic, top_n=10)
    (['reg8', 'reg7', 'reg6', 'reg5'], ['reg4', 'reg3', 'reg2', 'reg1'])
    >>> split_regions_by_sc(reg2sc_dic, bottom_n=2)
    (['reg8', 'reg7', 'reg6', 'reg5', 'reg4', 'reg3'], ['reg2', 'reg1'])
    >>> split_regions_by_sc(reg2sc_dic, bottom_n=10)
    (['reg8', 'reg7', 'reg6', 'reg5'], ['reg4', 'reg3', 'reg2', 'reg1'])
    >>> split_regions_by_sc(reg2sc_dic, top_n=2, bottom_n=2)
    (['reg8', 'reg7'], ['reg2', 'reg1'])
    >>> split_regions_by_sc(reg2sc_dic, top_n=10, bottom_n=10)
    (['reg8', 'reg7', 'reg6', 'reg5'], ['reg4', 'reg3', 'reg2', 'reg1'])
    >>> reg2sc_dic = {}
    >>> reg2sc_dic["reg1"] = 0.1
    >>> reg2sc_dic["reg2"] = 0.2
    >>> reg2sc_dic["reg3"] = 0.3
    >>> reg2sc_dic["reg4"] = 0.4
    >>> reg2sc_dic["reg5"] = 0.5
    >>> reg2sc_dic["reg6"] = 0.6
    >>> reg2sc_dic["reg7"] = 0.7
    >>> split_regions_by_sc(reg2sc_dic)
    (['reg7', 'reg6', 'reg5'], ['reg4', 'reg3', 'reg2', 'reg1'])
    >>> split_regions_by_sc(reg2sc_dic, top_n=2)
    (['reg7', 'reg6'], ['reg5', 'reg4', 'reg3', 'reg2', 'reg1'])
    >>> split_regions_by_sc(reg2sc_dic, top_n=10)
    (['reg7', 'reg6', 'reg5'], ['reg4', 'reg3', 'reg2', 'reg1'])
    >>> split_regions_by_sc(reg2sc_dic, bottom_n=2)
    (['reg7', 'reg6', 'reg5', 'reg4', 'reg3'], ['reg2', 'reg1'])
    >>> split_regions_by_sc(reg2sc_dic, bottom_n=10)
    (['reg7', 'reg6', 'reg5'], ['reg4', 'reg3', 'reg2', 'reg1'])
    >>> split_regions_by_sc(reg2sc_dic, top_n=2, bottom_n=2)
    (['reg7', 'reg6'], ['reg2', 'reg1'])
    >>> split_regions_by_sc(reg2sc_dic, top_n=10, bottom_n=10)
    (['reg7', 'reg6', 'reg5'], ['reg4', 'reg3', 'reg2', 'reg1'])
    >>> split_regions_by_sc(reg2sc_dic, rev_sort=False)
    (['reg1', 'reg2', 'reg3'], ['reg4', 'reg5', 'reg6', 'reg7'])
    >>> split_regions_by_sc(reg2sc_dic, top_n=2, bottom_n=2, rev_sort=False)
    (['reg1', 'reg2'], ['reg6', 'reg7'])
    >>> split_regions_by_sc(reg2sc_dic, bottom_n=2, rev_sort=False)
    (['reg1', 'reg2', 'reg3', 'reg4', 'reg5'], ['reg6', 'reg7'])
    
    """

    sorted_reg2sc = sorted(reg2sc_dic.items(), key=lambda item: item[1], reverse=rev_sort)

    midpoint = len(sorted_reg2sc) // 2
    top = []
    bottom = []

    if top_n and bottom_n:
        if top_n + bottom_n > len(sorted_reg2sc):
            top = sorted_reg2sc[:midpoint]
            bottom = sorted_reg2sc[midpoint:]
        else:
            top = sorted_reg2sc[:top_n]
            bottom = sorted_reg2sc[-bottom_n:]
    elif top_n and not bottom_n:
        if top_n > len(sorted_reg2sc):
            top = sorted_reg2sc[:midpoint]
            bottom = sorted_reg2sc[midpoint:]
        else:
            top = sorted_reg2sc[:top_n]
            bottom = sorted_reg2sc[top_n:]
    elif bottom_n and not top_n:
        if bottom_n > len(sorted_reg2sc):
            top = sorted_reg2sc[:midpoint]
            bottom = sorted_reg2sc[midpoint:]
        else:
            top = sorted_reg2sc[:-bottom_n]
            bottom = sorted_reg2sc[-bottom_n:]
    else:
        top = sorted_reg2sc[:midpoint]
        bottom = sorted_reg2sc[midpoint:]

    top_ids = [item[0] for item in top]
    bottom_ids = [item[0] for item in bottom]

    return top_ids, bottom_ids


################################################################################

def calc_exp_kmer_perc(kmer_k):
    """
    Calculate k-mer percentage in case of a  uniform distribution.
    E.g.
    kmer_k = 1
    4^1 = 4 possible kmers
    1/4 = 0.25 = 25 kmer percentage (return this value).
    
    >>> calc_exp_kmer_perc(1)
    25.0
    >>> calc_exp_kmer_perc(2)
    6.25
    >>> calc_exp_kmer_perc(3)
    1.5625
    >>> calc_exp_kmer_perc(4)
    0.390625
    >>> calc_exp_kmer_perc(5)
    0.09765625
    >>> calc_exp_kmer_perc(0)
    100.0

    """

    kmer_perc = 1 / (4 ** kmer_k) * 100
    return kmer_perc


################################################################################

def min_max_scale(values, new_min=0, new_max=1):
    # Find the minimum and maximum of the input values
    min_val = min(values)
    max_val = max(values)

    # Handle edge case where all values are the same
    if min_val == max_val:
        return [new_min for _ in values]

    # Apply min-max scaling
    scaled_values = [
        new_min + (val - min_val) * (new_max - new_min) / (max_val - min_val)
        for val in values
    ]
    return scaled_values


################################################################################

def search_generate_html_report(args,
                                df_pval, pval_cont_lll,
                                search_rbps_dic,
                                id2name_dic, name2ids_dic,
                                region_rbp_motif_pos_dic,
                                reg2pol_dic,
                                benchlib_path,
                                rbp2regidx_dic,
                                reg_ids_list,
                                seq_len_df=None,
                                mrna_prof_dic=False,
                                ei_ol_stats_dic=False,
                                seq_motif_blocks_dic=None,
                                reg2annot_dic=False,
                                annot2color_dic=False,
                                html_report_out="report.rbpbench_search.html",
                                goa_results_df=False,
                                goa_stats_dic=False,
                                plotly_embed_style=1,
                                plotly_full_html=False,
                                rbpbench_mode="search",
                                disable_motif_enrich_table=False,
                                disable_top_kmers_plot=False,
                                reg_seq_str="regions",
                                reg2seq_dic=False,
                                reg2sc_dic=False,
                                add_annot_stats_dic=False,
                                plots_subfolder="html_report_plots"):
    """
    Create additional hit statistics for selected RBPs, 
    e.g. correlation / co-occurrence between RBPs.

    """

    # Use absolute paths?
    out_folder = args.out_folder
    if args.plot_abs_paths:
        out_folder = os.path.abspath(out_folder)

    # Version string.
    version_str = "v" + args.version

    plots_folder = plots_subfolder
    plots_out_folder = out_folder + "/" + plots_folder
    if args.plot_abs_paths:
        plots_folder = plots_out_folder

    # Delete folder if already present.
    if os.path.exists(plots_out_folder):
        shutil.rmtree(plots_out_folder)
    os.makedirs(plots_out_folder)

    html_out = out_folder + "/" + "report.rbpbench_search.html"
    md_out = out_folder + "/" + "report.rbpbench_search.md"
    if html_report_out:
        html_out = html_report_out

    # Number of genomic regions.
    c_regions = len(reg_ids_list)

    # Check if no regions have motif hits.
    no_region_hits = True
    for rbp_id in rbp2regidx_dic:
        if rbp2regidx_dic[rbp_id]:
            no_region_hits = False
            break

    site_type_uc = "Genomic"
    site_type = "genomic"
    if rbpbench_mode == "searchrna":
        site_type_uc = "Transcript"
        site_type = "transcript"

    regex_motif_info = ""
    if args.regex:
        regex_motif_info = "Used regex motif: '%s'." %(name2ids_dic[args.regex_id][0])

    # Number of input regions.
    c_in_regions = 0
    for rbp_id in search_rbps_dic:
        c_in_regions += search_rbps_dic[rbp_id].c_hit_reg
        c_in_regions += search_rbps_dic[rbp_id].c_no_hit_reg      
        break

    """
    Setup sorttable.js to make tables in HTML sortable.

    """
    sorttable_js_path = benchlib_path + "/content/sorttable.js"
    assert os.path.exists(sorttable_js_path), "sorttable.js not at %s" %(sorttable_js_path)
    sorttable_js_html = '<script src="' + sorttable_js_path + '" type="text/javascript"></script>'
    if args.sort_js_mode == 2:
        shutil.copy(sorttable_js_path, plots_out_folder)
        sorttable_js_path = plots_folder + "/sorttable.js"
        sorttable_js_html = '<script src="' + sorttable_js_path + '" type="text/javascript"></script>'
    elif args.sort_js_mode == 3:
        js_code = read_file_content_into_str_var(sorttable_js_path)
        sorttable_js_html = "<script>\n" + js_code + "\n</script>\n"

    """
    Setup plotly .js to support plotly plots.

    https://plotly.com/javascript/getting-started/#download
    plotly-latest.min.js
    Packaged version: plotly-2.20.0.min.js
    
    """

    include_plotlyjs = "cdn"
    # plotly_full_html = False
    plotly_js_html = ""
    plotly_js_path = benchlib_path + "/content/plotly-2.20.0.min.js"
    assert os.path.exists(plotly_js_path), "plotly .js %s not found" %(plotly_js_path)
    if args.plotly_js_mode == 2:
        include_plotlyjs = plotly_js_path
    elif args.plotly_js_mode == 3:
        shutil.copy(plotly_js_path, plots_out_folder)
        include_plotlyjs = "plotly-2.20.0.min.js" # Or plots_folder + "/plotly-2.20.0.min.js" ?
    elif args.plotly_js_mode == 4:
        include_plotlyjs = True
        # plotly_full_html = False # Don't really need full html (head body ..) in plotly html.
    elif args.plotly_js_mode == 5:
        plotly_js_web = "https://cdn.plot.ly/plotly-2.25.2.min.js"
        plotly_js_html = '<script src="' + plotly_js_web + '"></script>' + "\n"
        include_plotlyjs = False
        # plotly_full_html = True
    elif args.plotly_js_mode == 6:
        shutil.copy(plotly_js_path, plots_out_folder)
        plotly_js = plots_folder + "/plotly-2.20.0.min.js"
        plotly_js_html = '<script src="' + plotly_js + '"></script>' + "\n"
        include_plotlyjs = False
    elif args.plotly_js_mode == 7:
        js_code = read_file_content_into_str_var(plotly_js_path)
        plotly_js_html = "<script>\n" + js_code + "\n</script>\n"
        include_plotlyjs = False
        # plotly_full_html = True

    # Logo path.
    logo_path_html = plots_folder + "/logo.png"
    if not os.path.exists(logo_path_html):
        logo_path = benchlib_path + "/content/logo.png"
        assert os.path.exists(logo_path), "logo.png not found in %s" %(logo_path)
        shutil.copy(logo_path, plots_out_folder)

    # HTML head section.
    html_head = """<!DOCTYPE html>
<html>
<head>
<title>RBPBench - Search Report</title>
%s

<style>
    th, td {
        border: 1px solid black;
        padding: 8px;
        text-align: center;
    }
    th {
        background-color: #f2f2f2;
    }
    .page {
        page-break-after: always;
    }

    .title-container {
        display: flex;
        align-items: center;
    }
    .title-container img {
        margin-right: 10px; /* Adjust the spacing as needed */
    }

</style>

</head>

<div class="title-container">
    <img src="%s" alt="Logo" width="175">
    <h1>Search Report</h1>
</div>


<body>
""" %(plotly_js_html, logo_path_html)

    # HTML tail section.
    html_tail = """
%s
</body>

</html>
""" %(sorttable_js_html)

    motif_enrich_info = "- [RBP region score motif enrichment statistics](#rbp-enrich-stats)"
    if disable_motif_enrich_table:
        motif_enrich_info = ""

    # Markdown part.
    mdtext = """

List of available statistics and plots generated
by RBPBench (%s, rbpbench %s):

%s
- [RBP motif co-occurrences heat map](#cooc-heat-map)""" %(version_str, rbpbench_mode, motif_enrich_info)

    mdtext += "\n"

    # Additional plot if GTF annotations given.
    if seq_len_df is not None:
        mdtext += "- [Input sequence length distribution](#seq-len-plot)\n"
    if not disable_top_kmers_plot:
        mdtext += "- [Input sequences top k-mers plot](#seqs-top-kmer-plot)\n"
    if not args.disable_kmer_tb_plot:
        mdtext += "- [Top vs bottom scoring regions k-mer distribution](#kmer-dist)\n"
    if reg2seq_dic and not args.disable_kmer_var_plot:
        mdtext += "- [Input sequences k-mer variation plot](#seq-var-plot)\n"
    if not args.disable_kmer_pca_plot:
        mdtext += "- [Input sequences by k-mer content plot](#input-seqs-kmer)\n"

    # Additional plot if GTF annotations given.
    if reg2annot_dic:
        mdtext += "- [Region annotations per RBP](#annot-rbp-plot)\n"
    # Additional region annotations statistics.
    if add_annot_stats_dic:
        mdtext += "- [Additional region annotations](#add-annot-stats)\n"

    if mrna_prof_dic:
        mdtext += "- [mRNA region coverage profile](#mrna-prof-plot)\n"

    if ei_ol_stats_dic:
        mdtext += "- [Exon-intron overlap statistics](#ei-ol-stats)\n"

    if reg2annot_dic:
        mdtext += "- [Intronic regions intron border distances](#intron-border-dist)\n"

    # Upset plot.
    # mdtext += "\n"
    if args.enable_upset_plot:
        mdtext += "- [RBP combinations upset plot](#rbp-comb-upset-plot)\n"

    # If --set-rbp-id given.
    if args.set_rbp_id is not None:
        # if reg2annot_dic is None:
        #     mdtext += "\n"
        mdtext += "- [Set RBP %s motifs distance statistics](#rbp-motif-dist-stats)\n" %(args.set_rbp_id)
        # mdtext += "- [Set RBP %s motif distance plot](#rbp-motif-dist-plot)\n" %(set_rbp_id)
        for idx, motif_id in enumerate(name2ids_dic[args.set_rbp_id]):
            mdtext += "    - [Motif %s distance statistics](#single-motif-%i-dist-stats)\n" %(motif_id, idx)
            # mdtext += "    - [Single motif %s distance plot](#single-motif-%i-dist-plot)\n" %(motif_id, idx)

    if args.run_goa:
        mdtext += "- [GO enrichment analysis results](#goa-results)\n"

    add_head_info = ""
    if args.bed_sc_thr is not None:
        add_head_info = " BED score threshold (--bed-sc-thr) = %s" %(str(args.bed_sc_thr))
        if args.bed_sc_thr_rev_filter:
            add_head_info += " (reverse filtering applied, i.e., the lower the better)."
        else:
            add_head_info += "."

    mdtext += "\nFIMO p-value threshold (--fimo-pval) = %s.%s # of considered input %s = %i. Region extension (upstream, downstream) = (%i, %i).\n" %(str(args.fimo_pval), add_head_info, reg_seq_str, c_in_regions, args.ext_up, args.ext_down)
    mdtext += "\n&nbsp;\n"


    """
    RBP region score motif enrichment statistics

    """

    if not disable_motif_enrich_table:

        # Inform about set alterntive hypothesis for Wilcoxon rank sum test.
        wrs_mode_info1 = "Wilcoxon rank-sum test alternative hypothesis is set to 'greater', i.e., low p-values mean motif-containing regions have significantly higher scores."
        wrs_mode_info2 = "higher"
        if args.wrs_mode == 2:
            wrs_mode_info1 = "Wilcoxon rank-sum test alternative hypothesis is set to 'less', i.e., low p-values mean motif-containing regions have significantly lower scores."
            wrs_mode_info2 = "lower"

        mdtext += """
## RBP region score motif enrichment statistics ### {#rbp-enrich-stats}

**Table:** RBP region score motif enrichment statistics. Given a score for each genomic region (# input regions = %i), 
RBPbench checks whether motif-containing regions have significantly different scores compared to regions without motif hits.
%s
In other words, a low test p-value for a given RBP indicates 
that %s-scoring regions are more likely to contain motif hits of the respective RBP.
**NOTE** that if scores associated with input genomic regions are all the same, p-values become meaningless 
(i.e., they result in p-values of 1.0).
Likewise, the p-value becomes non-informative if most or all input regions have RBP motif hits (i.e., very high hit region percentages).
By default, BED genomic regions input file column 5 is used as the score column (change with --bed-score-col).
%s

""" %(c_in_regions, wrs_mode_info1, wrs_mode_info2, regex_motif_info)


        pval_dic = {}
        for rbp_id in search_rbps_dic:
            wc_pval = search_rbps_dic[rbp_id].wc_pval
            pval_dic[rbp_id] = wc_pval

        # for rbp_id, wc_pval in sorted(pval_dic.items(), key=lambda item: item[1], reverse=False):
        #     wc_pval_str = convert_sci_not_to_decimal(wc_pval)  # Convert scientific notation to decimal string for sorting to work.
        #     c_hit_reg = search_rbps_dic[rbp_id].c_hit_reg
        #     perc_hit_reg = search_rbps_dic[rbp_id].perc_hit_reg
        #     c_uniq_motif_hits = search_rbps_dic[rbp_id].c_uniq_motif_hits
        #     mdtext += '<tr>' + "\n"
        #     mdtext += '<td style="text-align: center;"' + ">%s</td>\n" %(rbp_id)
        #     mdtext += '<td style="text-align: center;"' + ">%i</td>\n" %(c_hit_reg) 
        #     mdtext += '<td style="text-align: center;"' + ">%.2f</td>\n" %(perc_hit_reg) 
        #     mdtext += '<td style="text-align: center;"' + ">%i</td>\n" %(c_uniq_motif_hits) 
        #     mdtext += '<td style="text-align: center;"' + ">" + str(wc_pval) + "</td>\n"
        #     mdtext += '</tr>' + "\n"

        # mdtext += '| RBP ID | # hit regions | % hit regions | # motif hits | p-value |' + " \n"
        # mdtext += "| :-: | :-: | :-: | :-: | :-: |\n"

        mdtext += '<table style="max-width: 800px; width: 100%; border-collapse: collapse; line-height: 0.8;">' + "\n"
        mdtext += "<thead>\n"
        mdtext += "<tr>\n"
        mdtext += "<th>RBP ID</th>\n"
        mdtext += "<th># hit regions</th>\n"
        mdtext += "<th>% hit regions</th>\n"
        mdtext += "<th># motif hits</th>\n"
        mdtext += "<th>p-value</th>\n"
        mdtext += "</tr>\n"
        mdtext += "</thead>\n"
        mdtext += "<tbody>\n"

        for rbp_id, wc_pval in sorted(pval_dic.items(), key=lambda item: item[1], reverse=False):
            wc_pval_str = convert_sci_not_to_decimal(wc_pval)  # Convert scientific notation to decimal string for sorting to work.
            c_hit_reg = search_rbps_dic[rbp_id].c_hit_reg
            perc_hit_reg = search_rbps_dic[rbp_id].perc_hit_reg
            c_uniq_motif_hits = search_rbps_dic[rbp_id].c_uniq_motif_hits
            # mdtext += "| %s | %i | %.2f | %i | %s |\n" %(rbp_id, c_hit_reg, perc_hit_reg, c_uniq_motif_hits, wc_pval)
            mdtext += '<tr>' + "\n"
            mdtext += "<td>" + rbp_id + "</td>\n"
            mdtext += "<td>" + str(c_hit_reg) + "</td>\n"
            mdtext += "<td>%.2f" %(perc_hit_reg) + "</td>\n"
            mdtext += "<td>" + str(c_uniq_motif_hits) + "</td>\n"
            mdtext += "<td>" + str(wc_pval) + "</td>\n"
            mdtext += '</tr>' + "\n"

        mdtext += '</tbody>' + "\n"
        mdtext += '</table>' + "\n"

        mdtext += "\n&nbsp;\n&nbsp;\n"
        mdtext += "\nColumn IDs have the following meanings: "
        mdtext += "**RBP ID** -> RBP ID from database or user-defined (typically RBP name), "
        mdtext += '**# hit regions** -> number of input genomic regions with motif hits (after filtering and optional extension), '
        mdtext += '**% hit regions** -> percentage of hit regions over all regions (i.e., how many input regions contain >= 1 RBP binding motif), '
        mdtext += '**# motif hits** -> number of unique motif hits in input regions (removed double counts), '
        mdtext += '**p-value** -> Wilcoxon rank-sum test p-value.' + "\n"
        mdtext += "\n&nbsp;\n"



    """
    Co-occurrence heat map.

    """

    cooc_plot_plotly =  "co-occurrence_plot.plotly.html"
    cooc_plot_plotly_out = plots_out_folder + "/" + cooc_plot_plotly

    create_cooc_plot_plotly(df_pval, pval_cont_lll, cooc_plot_plotly_out,
                            max_motif_dist=args.max_motif_dist,
                            min_motif_dist=args.min_motif_dist,
                            include_plotlyjs=include_plotlyjs,
                            full_html=plotly_full_html)

    plot_path = plots_folder + "/" + cooc_plot_plotly

    mdtext += """
## RBP motif co-occurrences heat map ### {#cooc-heat-map}

"""
    if args.plotly_js_mode in [5, 6, 7]:
        # Read in plotly code.
        # mdtext += '<div style="width: 1200px; height: 1200px; align-items: center;">' + "\n"
        js_code = read_file_content_into_str_var(cooc_plot_plotly_out)
        js_code = js_code.replace("height:100%; width:100%;", "height:1200px; width:1200px;")
        mdtext += js_code + "\n"
        # mdtext += '</div>'
    else:
        if plotly_embed_style == 1:
            # mdtext += '<div class="container-fluid" style="margin-top:40px">' + "\n"
            mdtext += "<div>\n"
            mdtext += '<iframe src="' + plot_path + '" width="1200" height="1200"></iframe>' + "\n"
            mdtext += '</div>'
        elif plotly_embed_style == 2:
            mdtext += '<object data="' + plot_path + '" width="1200" height="1200"> </object>' + "\n"

    p_val_info = "P-values below %s are considered significant." %(str(args.cooc_pval_thr))

    min_motif_dist_info = ""
    if args.min_motif_dist > 0:
        min_motif_dist_info = " + a mean minimum motif distance >= %i nt " %(args.min_motif_dist)

    if args.cooc_pval_mode == 1:
        p_val_info = "RBP co-occurrences with Benjamini-Hochberg multiple testing corrected p-values below %s %sare considered significant." %(str(args.cooc_pval_thr), min_motif_dist_info)
    elif args.cooc_pval_mode == 2:
        p_val_info = "RBP co-occurrences with p-values below %s (p-value threshold Bonferroni multiple testing corrected) %sare considered significant." %(str(args.cooc_pval_thr), min_motif_dist_info)
    elif args.cooc_pval_mode == 3:
        p_val_info = "RBP co-occurrences with p-values below %s %sare considered significant." %(str(args.cooc_pval_thr), min_motif_dist_info)
    else:
        assert False, "Invalid co-occurrence p-value mode (--cooc-pval-mode) set: %i" %(args.cooc_pval_mode)
    
    p_val_info += " # of RBP co-occurrence comparisons: %i. # of significant co-occurrences: %i (%.2f%%)." %(args.c_all_fisher_pval, args.c_sig_fisher_pval, args.perc_sig_fisher_pval)

    # Inform about set alterntive hypothesis for Fisher exact test on RBP motif co-occurrences.
    fisher_mode_info = "Fisher exact test alternative hypothesis is set to 'greater', i.e., significantly overrepresented co-occurrences are reported."
    if args.fisher_mode == 2:
        fisher_mode_info = "Fisher exact test alternative hypothesis is set to 'two-sided', i.e., significantly over- and underrepresented co-occurrences are reported."
    elif args.fisher_mode == 3:
        fisher_mode_info = "Fisher exact test alternative hypothesis is set to 'less', i.e., significantly underrepresented co-occurrences are reported."

    mdtext += """

**Figure:** Heat map of motif hit co-occurrences (Fisher's exact test p-values) between RBPs. 
RBP motif hit co-occurrences that are not significant are colored white, while
significant co-occurrences are colored according to their -log10 p-value (used as legend color, i.e., the higher the more significant).
%s
%s
Hover box: 
**1)** RBP1.
**2)** RBP2.
**3)** p-value: Fisher's exact test p-value (calculated based on contingency table (6) between RBP1 and RBP2). 
**4)** p-value after filtering: p-value after filtering, i.e., p-value is kept if significant (< %s), otherwise it is set to 1.0.
**5)** RBPs compaired. 
**6)** Counts[]: contingency table of co-occurrence counts (i.e., number of %s regions with/without shared motif hits) between compaired RBPs, 
with format [[A, B], [C, D]], where 
A: RBP1 AND RBP2, 
B: NOT RBP1 AND RBP2
C: RBP1 AND NOT RBP2
D: NOT RBP1 AND NOT RBP2.
**7)** Mean minimum distance of RBP1 and RBP2 motifs (mean over all regions containing RBP1 + RBP2 motifs). Distances measured from motif center positions.
**8)** Over all regions containing RBP1 and RBP2 motif pairs, percentage of regions where RBP1 + RBP2 motifs are within %i nt distance (set via --max-motif-dist).
**9)** Correlation: Pearson correlation coefficient between RBP1 and RBP2.
%s regions are labelled 1 or 0 (RBP motif present or not), resulting in a vector of 1s and 0s for each RBP.
Correlations are then calculated by comparing vectors for every pair of RBPs.
**10)** -log10 of p-value after filtering, used for legend coloring. Using p-value after filtering, all non-significant p-values become 0 
for easier distinction between significant and non-significant co-occurrences.
%s

&nbsp;

""" %(p_val_info, fisher_mode_info, str(args.cooc_pval_thr), site_type, args.max_motif_dist, site_type_uc, regex_motif_info)


    """
    Input sequence lengths distribution violin plot.
    
    """


    # Additional plot if GTF annotations given.
    if seq_len_df is not None:

        mdtext += """
## Input sequence length distribution ### {#seq-len-plot}

"""

        seq_len_plot_plotly =  "seq_len_violin_plot.plotly.html"
        seq_len_plot_plotly_out = plots_out_folder + "/" + seq_len_plot_plotly

        create_len_distr_violin_plot_plotly(seq_len_df, seq_len_plot_plotly_out,
                                include_plotlyjs=include_plotlyjs,
                                full_html=plotly_full_html)

        plot_path = plots_folder + "/" + seq_len_plot_plotly

        if args.plotly_js_mode in [5, 6, 7]:
            # Read in plotly code.
            # mdtext += '<div style="width: 1200px; height: 1200px; align-items: center;">' + "\n"
            js_code = read_file_content_into_str_var(seq_len_plot_plotly_out)
            js_code = js_code.replace("height:100%; width:100%;", "height:1000px; width:1000px;")
            mdtext += js_code + "\n"
            # mdtext += '</div>'
        else:
            if plotly_embed_style == 1:
                # mdtext += '<div class="container-fluid" style="margin-top:40px">' + "\n"
                mdtext += "<div>\n"
                mdtext += '<iframe src="' + plot_path + '" width="1000" height="1000"></iframe>' + "\n"
                mdtext += '</div>'
            elif plotly_embed_style == 2:
                mdtext += '<object data="' + plot_path + '" width="1000" height="1000"> </object>' + "\n"

        mdtext += """

**Figure:** Input sequence lengths (after filtering and optional extension) violin plot.
Hover box over data points shows sequence ID, sequence (50 nt per line), sequence length, and motif hits. 
Motif hit format: motif ID, motif start - motif end.
Violin plot shows density distribution of sequence lengths, including min, max, median, 
and length quartiles (q1: 25th percentile, q3: 75th percentile).

&nbsp;

"""

    if not disable_top_kmers_plot:

        assert reg2seq_dic, "No region sequences supplied to report function"

        mdtext += """
## Input sequences top k-mers plot ### {#seqs-top-kmer-plot}

"""

        seqs_kmer_dic = seqs_dic_count_kmer_freqs(reg2seq_dic, args.kmer_plot_k, 
                                                  rna=False,
                                                  return_ratios=True,
                                                  perc=True,
                                                  report_key_error=False,
                                                  skip_non_dic_keys=True,
                                                  convert_to_uc=True)

        n_top_kmers = 20

        top_kmers_list = []
        top_kmer_perc_list = []

        for kmer, perc in sorted(seqs_kmer_dic.items(), key=lambda x: x[1], reverse=True)[:n_top_kmers]:
            top_kmers_list.append(kmer)
            top_kmer_perc_list.append(perc)

        # Average percentage in seqs_kmer_dic.
        avg_perc = sum(seqs_kmer_dic.values()) / len(seqs_kmer_dic)

        fig, ax = plt.subplots(figsize=(11, 3.5))

        ax.bar(top_kmers_list, top_kmer_perc_list, color='#e5ecf6', zorder=2)  # #e5ecf6 lightgray

        ax.set_ylabel(str(args.kmer_plot_k) + "-mer %")

        # Expected k-mer percentage.
        exp_kmer_perc = calc_exp_kmer_perc(args.kmer_plot_k)

        # Add horizontal line for expected uniform percentage.
        ax.axhline(y=exp_kmer_perc, color='red', linestyle='--', linewidth=1, label=f'Expected uniform % ({exp_kmer_perc:.2f}%)')

        # ax.set_yticks(range(0, 101, 10))

        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['bottom'].set_visible(False)

        ax.grid(axis='y', color='#e5ecf6', linestyle='-', alpha=0.7, linewidth=0.7, zorder=1)

        plt.xticks(rotation=-45, ha='center')

        plt.subplots_adjust(left=0.06, right=0.99, top=0.95, bottom=0.20)
        # plt.subplots_adjust(left=0.05, right=0.95, top=0.9, bottom=0.25)

        top_kmers_plot = "seqs_top_kmers.bar_plot.png"
        top_kmers_plot_out = plots_out_folder + "/" + top_kmers_plot
        plot_path = plots_folder + "/" + top_kmers_plot

        plt.savefig(top_kmers_plot_out, dpi=140)

        if args.plot_pdf and top_kmers_plot_out.endswith('.png'):
            pdf_out = top_kmers_plot_out[:-4] + '.pdf'
            plt.savefig(pdf_out, dpi=140)

        plt.close()
        mdtext += '<image src = "' + plot_path + '" width="1100px"></image>'  + "\n"

        mdtext += """
**Figure:** Top %i %i-mers encountered in input sequences. 
In case of a uniform distribution with all %i-mers present, each %i-mer would have a percentage of %s%% (marked by red line).
Change k via --kmer-plot-k.

&nbsp;

""" %(n_top_kmers, args.kmer_plot_k, args.kmer_plot_k, args.kmer_plot_k, exp_kmer_perc)


    """
    Top vs bottom scoring regions k-mer distribution.

    reg2sc_dic:
        region ID -> score dictionary.

    """
    if not args.disable_kmer_tb_plot:

        assert reg2sc_dic, "No region scores supplied (reg2sc_dic missing) for plotting top vs bottom scoring regions k-mer plot"

        top_seqs_dic = {}
        bottom_seqs_dic = {}

        rev_sort = True  # True if scores (i.e. the higher the better site quality).
        if args.bed_sc_thr_rev_filter:  # If scores are e.g. p-values, reverse filtering.
            rev_sort = False

        top_ids, bottom_ids = split_regions_by_sc(reg2sc_dic, 
                                                  top_n=args.kmer_tb_plot_top_n, 
                                                  bottom_n=args.kmer_tb_plot_bottom_n,
                                                  rev_sort=rev_sort)

        c_top_sites = len(top_ids)
        c_bottom_sites = len(bottom_ids)

        for reg_id in top_ids:
            top_seqs_dic[reg_id] = reg2seq_dic[reg_id]
        for reg_id in bottom_ids:
            bottom_seqs_dic[reg_id] = reg2seq_dic[reg_id]

        top_kmer_dic = seqs_dic_count_kmer_freqs(top_seqs_dic, args.kmer_plot_k, 
                                                 rna=False,
                                                 return_ratios=True,
                                                 perc=True,
                                                 report_key_error=False,
                                                 skip_non_dic_keys=True,
                                                 convert_to_uc=True)
        bottom_kmer_dic = seqs_dic_count_kmer_freqs(bottom_seqs_dic, args.kmer_plot_k, 
                                                    rna=False,
                                                    return_ratios=True,
                                                    perc=True,
                                                    report_key_error=False,
                                                    skip_non_dic_keys=True,
                                                    convert_to_uc=True)


        mdtext += """
## Top vs bottom scoring regions k-mer distribution ### {#kmer-dist}

"""

        plotly_kmer_plot = "plotly_scatter_kmer.html"
        plotly_kmer_plot_out = plots_out_folder + "/" + plotly_kmer_plot

        # Create k-mer plotly scatter plot.
        create_kmer_sc_plotly_scatter_plot(top_kmer_dic, bottom_kmer_dic, args.kmer_plot_k,
                                           plotly_kmer_plot_out, plotly_js_path,
                                           pos_label=f"{args.kmer_plot_k}-mer % top scoring sites",
                                           neg_label=f"{args.kmer_plot_k}-mer % bottom scoring sites",
                                           kmer_label=f"{args.kmer_plot_k}-mer")

        # Plot paths inside html report.
        plotly_kmer_plot_path = plots_folder + "/" + plotly_kmer_plot

        # R2 score.
        r2_kmer = calc_r2_corr_measure(top_kmer_dic, bottom_kmer_dic,
                                       is_dic=True)

        # Expected k-mer percentage.
        exp_kmer_perc = calc_exp_kmer_perc(args.kmer_plot_k)

        if args.plotly_js_mode in [5, 6, 7]:
            js_code = read_file_content_into_str_var(plotly_kmer_plot_out)
            js_code = js_code.replace("height:100%; width:100%;", "height:800px; width:800px;")
            mdtext += js_code + "\n"
        else:
            mdtext += "<div>\n"
            mdtext += '<iframe src="' + plotly_kmer_plot_path + '" width="800" height="800"></iframe>' + "\n"
            mdtext += '</div>'

        mdtext += """

**Figure:** Sequence %i-mer percentages in the top %i scoring and bottom %i scoring input sites. In case of
a uniform distribution with all %i-mers present, each %i-mer would have a percentage of %s%%. R2 = %.6f.
By default, BED genomic regions input file column 5 is used as the score column (change with --bed-score-col).
**NOTE** that this plot is only meaningful if the provided scores are related to binding site quality / binding affinity.
In this case, the plot can hint at RBP binding preferences (corresponding to enriched k-mers in the top input sites). 
Plot can be further modified by specifying top and bottom n scoring regions (--kmer-tb-plot-top-n, --kmer-tb-plot-bottom-n, 
by default top and bottom are split in half), or change k (--kmer-plot-k).

&nbsp;

""" %(args.kmer_plot_k, c_top_sites, c_bottom_sites, args.kmer_plot_k, args.kmer_plot_k, str(exp_kmer_perc), r2_kmer)


    """
    Input sequences k-mer variation plot.

    """

    if reg2seq_dic and not args.disable_kmer_var_plot:

        mdtext += """
## Input sequences k-mer variation plot ### {#seq-var-plot}


"""
        remove_zero_val = True

        kmer2stats_dic = {}
        single_cv_dic, avg_cv = calculate_k_nucleotide_cv(reg2seq_dic,
                                                          k=args.kmer_plot_k,
                                                          nucleotides=['A', 'C', 'G', 'T'],
                                                          kmer2stats_dic=kmer2stats_dic,
                                                          reg2sc_dic=reg2sc_dic,
                                                          only_observed=True)

        seq_var_plot_plotly =  "seq_kmer_var_plot.plotly.html"
        seq_var_plot_plotly_out = plots_out_folder + "/" + seq_var_plot_plotly

        create_seq_var_violin_plot_plotly(kmer2stats_dic, single_cv_dic, avg_cv, 
                                          seq_var_plot_plotly_out,
                                          kmer_size=args.kmer_plot_k,
                                          remove_zero_val=remove_zero_val,
                                          color_mode=args.kmer_var_color_mode,
                                          include_plotlyjs=include_plotlyjs,
                                          full_html=plotly_full_html)

        plot_path = plots_folder + "/" + seq_var_plot_plotly

        if args.plotly_js_mode in [5, 6, 7]:
            js_code = read_file_content_into_str_var(seq_var_plot_plotly_out)
            js_code = js_code.replace("height:100%; width:100%;", "height:1000px; width:1000px;")
            mdtext += js_code + "\n"
        else:
            if plotly_embed_style == 1:
                mdtext += "<div>\n"
                mdtext += '<iframe src="' + plot_path + '" width="1000" height="1000"></iframe>' + "\n"
                mdtext += '</div>'
            elif plotly_embed_style == 2:
                mdtext += '<object data="' + plot_path + '" width="1000" height="1000"> </object>' + "\n"

        mdtext += """

**Figure:** k-mer variations in input sequences. Set k-mer size = %i (change via --kmer-plot-k).
x-axis: k-mer variation, describing the variation of the k-mer over all input sequences (see coefficient of variation (CV) in hover box description).
y-axis: site %%, i.e., percentage of sequences where k-mer is present.
By default, the correlation (Spearman correlation coefficient) between k-mer ratios and input region scores is used for coloring 
(alternatively use --seq-var-mode to change to k-mer %%).
Hover box:
**k-mer** -> observed k-mer.
**Variation** -> coefficient of variation (CV) of the k-mer ratios over all input site sequences. 
CV is defined as the standard deviation of ratios divided by the mean ratio. 
Thus, the lower the CV of a k-mer, the more even its ratios over the input site sequences.
A CV of 0.0 for a given k-mer would mean that all sequences have exactly the same ratio of the k-mer.
In general, we can assume that k-mers with high CVs are less important for the RBP binding, since
they do not occur evenly over the dataset. In contrast, k-mers with low CVs are more 
evenly distributed over the dataset sequences, and thus should be more important for RBP binding 
(assuming a sufficient dataset quality + an affinity of the RBP towards specific k-mer sequences). 
**k-mer #** -> number of times k-mer appears in all input sequences.
**k-mer %%** -> percentage of k-mer over all k-mers in the input sequences.
**Site %%** -> percentage of sequences where k-mer is present.
**Correlation** -> correlation coefficient (Spearman) between k-mer ratios and input region scores.

&nbsp;

""" %(args.kmer_plot_k)


    # Region ID to motifs string for plotting.
    reg2motifs_dic = {}  # If false motif info not added to hover box.

    # Snatch motif hit info from seq_len_df.
    if not args.kmer_pca_plot_no_motifs:
        try:
            assert not seq_len_df.empty, "seq_len_df dataframe is defined but empty"
        except NameError:
            raise NameError("seq_len_df dataframe is not defined")

        reg2motifs_dic = seq_len_df.set_index('Sequence ID')['Motif hits'].to_dict()

    if not args.disable_kmer_pca_plot:

        mdtext += """
## Input sequences by k-mer content plot ### {#input-seqs-kmer}

"""
        seqs_kmer_k = args.kmer_pca_plot_k
        add_seq_comp = True
        if args.kmer_pca_plot_no_comp:
            add_seq_comp = False
        min_max_norm_seq_comp = False

        # Get kmer -> count dictionary.
        # kmer2c_dic = get_kmer_dic(seqs_kmer_k, rna=False)

        reg2ntps_dic = {}
        # Average mono-nucleotide percentages over whole dataset. 
        ntps_dic = {"A" : 0.0,
                    "C" : 0.0,
                    "G" : 0.0,
                    "T" : 0.0}

        reg2kmer_rat_dic = seqs_dic_get_kmer_ratios(reg2seq_dic, seqs_kmer_k,
                                                    rna=False,
                                                    report_key_error=False,
                                                    skip_non_dic_keys=True,
                                                    add_seq_comp=add_seq_comp,
                                                    add_gc_skew=False,
                                                    add_at_skew=False,
                                                    add_at_content=False,
                                                    add_1nt_ratios=False,
                                                    seq_comp_k=args.kmer_pca_plot_comp_k,
                                                    ntps_dic=ntps_dic,
                                                    reg2ntps_dic=reg2ntps_dic,
                                                    convert_to_uc=True)

        # Normalize sequence complexity value (min max normalization).
        if min_max_norm_seq_comp:
            last_values = []
            for reg_id in sorted(reg2kmer_rat_dic):
                last_values.append(reg2kmer_rat_dic[reg_id][-1])
            scaled_last_values = min_max_scale(last_values)
            i = 0
            for reg_id in sorted(reg2kmer_rat_dic):
                reg2kmer_rat_dic[reg_id][-1] = scaled_last_values[i]
                i += 1
        
        # Ratios ready?
        if reg2kmer_rat_dic:

            seqs_kmer_plot_plotly =  "seqs_kmer_plot.plotly.html"
            seqs_kmer_plot_plotly_out = plots_out_folder + "/" + seqs_kmer_plot_plotly

            create_seqs_kmer_plot_plotly(reg2seq_dic, reg2kmer_rat_dic, seqs_kmer_plot_plotly_out,
                                         k=seqs_kmer_k,
                                         annot2color_dic=annot2color_dic,
                                         reg2annot_dic=reg2annot_dic,
                                         reg2motifs_dic=reg2motifs_dic,
                                         reg2ntps_dic=reg2ntps_dic,
                                         include_plotlyjs=include_plotlyjs,
                                         full_html=plotly_full_html)

            plot_path = plots_folder + "/" + seqs_kmer_plot_plotly

            plot_width = 1000
            if reg2annot_dic:
                plot_width = 1200

            if args.plotly_js_mode in [5, 6, 7]:
                # Read in plotly code.
                # mdtext += '<div style="width: 1200px; height: 1200px; align-items: center;">' + "\n"
                js_code = read_file_content_into_str_var(seq_len_plot_plotly_out)
                js_code = js_code.replace("height:100%; width:100%;", "height:1000px; width:%ipx;" %(plot_width))
                mdtext += js_code + "\n"
                # mdtext += '</div>'
            else:
                if plotly_embed_style == 1:
                    # mdtext += '<div class="container-fluid" style="margin-top:40px">' + "\n"
                    mdtext += "<div>\n"
                    # mdtext += '<iframe src="' + plot_path + '" width="1000" height="1000"></iframe>' + "\n"
                    mdtext += '<iframe src="' + plot_path + '" width="%i" height="1000"></iframe>' %(plot_width) + "\n"
                    mdtext += '</div>'
                elif plotly_embed_style == 2:
                    # mdtext += '<object data="' + plot_path + '" width="1000" height="1000"> </object>' + "\n"
                    mdtext += '<iframe src="' + plot_path + '" width="%i" height="1000"></iframe>' %(plot_width) + "\n"

            add_info = ""
            if add_seq_comp:
                add_info = "Apart from the k-mer features, sequence complexity (i.e., Shannon entropy) is used as an additional feature for each sequence."

            add_avg_mono_nt_perc = ""
            if ntps_dic:
                add_avg_mono_nt_perc = "Average mono-nucleotide percentages over all input sequences: A = %.2f%%, C = %.2f%%, G = %.2f%%, T = %.2f%%." %(ntps_dic["A"], ntps_dic["C"], ntps_dic["G"], ntps_dic["T"])

            mdtext += """

**Figure:** Input sequences by k-mer content plot (k = %i, change via --kmer-pca-plot-k). 
Input sequences are represented by their k-mer content and visualized as 2D PCA plot.
The closer two sequences are, the more similar their k-mer content.
%s
%s

&nbsp;

""" %(args.kmer_pca_plot_k, add_info, add_avg_mono_nt_perc)

        else:
            mdtext += """

No sequence k-mers content plot generated since no k-mer contents extracted from input sequences (not enough sequences or sequences shorter than set --kmer-pca-plot-k?).

&nbsp;

"""



    """
    Region annotations per RBP plot.

    """

    if reg2annot_dic:

        annot_stacked_bars_plot =  "annotation_stacked_bars_plot.png"
        annot_stacked_bars_plot_out = plots_out_folder + "/" + annot_stacked_bars_plot

        mdtext += """
## Region annotations per RBP ### {#annot-rbp-plot}

"""
        if no_region_hits:
            mdtext += """

No plot generated since no motif hits found in input regions.

&nbsp;

"""
        else:

            add_all_reg_bar = True
            if args.disable_all_reg_bar:
                add_all_reg_bar = False

            create_search_annotation_stacked_bars_plot(rbp2regidx_dic, reg_ids_list, reg2annot_dic,
                                                       plot_out=annot_stacked_bars_plot_out,
                                                       annot2color_dic=annot2color_dic,
                                                       plot_pdf=args.plot_pdf,
                                                       add_all_reg_bar=add_all_reg_bar)

            plot_path = plots_folder + "/" + annot_stacked_bars_plot

            more_infos = ""
            if add_all_reg_bar:
                more_infos = "**All**: annotations for all input regions (with and without motif hits)."

            mdtext += '<img src="' + plot_path + '" alt="Annotation stacked bars plot"' + "\n"
            # mdtext += 'title="Annotation stacked bars plot" width="800" />' + "\n"
            mdtext += 'title="Annotation stacked bars plot" />' + "\n"
            mdtext += """
**Figure:** For each RBP, a stacked bar shows the corresponding region annotations 
(from --gtf GTF file, see legend for region types) for the genomic regions 
with motif hits for the respective RBP. 
Total bar height equals to the number of genomic regions with >= 1 motif hit for the RBP.
%s
Note that regions termed "intergenic" include all input regions not overlapping with selected transcript region features. 
By default, RBPBench for each gene in the GTF file selects the most prominent transcript (see manual for details).

&nbsp;

""" %(more_infos)


    """
    Additional region annotations.

    add_annot_stats_dic["c_genes"]
    add_annot_stats_dic["c_promoters"]
    add_annot_stats_dic["c_filt_min_tr_len"]
    add_annot_stats_dic["c_filt_mrna_only"]
    add_annot_stats_dic["c_outside_genes"]
    add_annot_stats_dic["c_inside_prom"]
    add_annot_stats_dic["c_add_annot"]
    add_annot_stats_dic["prom_ext_up"]
    add_annot_stats_dic["prom_ext_down"]
    c_in_regions == all considered input regions.

    """

    if add_annot_stats_dic:
        
        add_annot_bed_info = ""
        if args.add_annot_bed:
            add_annot_not_ol = " "
            if args.add_annot_comp:
                add_annot_not_ol = " NOT "
            add_annot_bed_info = 'Also included is the percentage of input regions' + add_annot_not_ol + 'overlapping with additionally provided regions (via --add-annot-bed, ID: ' + '"' + args.add_annot_id + '").'

        only_mrna_info = ""
        if args.prom_mrna_only:
            only_mrna_info = "Only mRNA transcripts were used for promoter region extraction (# transcripts removed = %i)." %(add_annot_stats_dic["c_filt_mrna_only"])
        min_tr_len_info = ""
        if args.prom_min_tr_len:
            min_tr_len_info = "Minimum transcript length for promoter region extraction = %i nt (# transcripts removed = %i)." %(args.prom_min_tr_len, add_annot_stats_dic["c_filt_min_tr_len"])

        mdtext += """
## Additional region annotations ### {#add-annot-stats}

**Table:** Percentages of input regions located outside of gene regions (annotated in provided GTF, # of considered
gene regions = %i) and input regions overlapping with putative promoter regions (taking the regions %i nt upstream 
to %i nt downstream of the transcript start sites). %s %s # of considered promoter regions = %i. %s High percentages 
of input regions located outside gene regions or 
inside promoter regions can point at dataset issues (assuming RBPs bind primarily to gene/transcript regions) 
or distinct protein functions (e.g., RBPs moonlighting as transcription factors). Note that depending 
on the methods used for dataset generation, input regions outside of gene regions might also have
been removed already. Number of considered input regions = %i.

""" %(add_annot_stats_dic["c_genes"], add_annot_stats_dic["prom_ext_up"], add_annot_stats_dic["prom_ext_down"], only_mrna_info, min_tr_len_info, add_annot_stats_dic["c_promoters"], add_annot_bed_info, c_in_regions)

        perc_outside_genes = 0.0
        if add_annot_stats_dic["c_outside_genes"] > 0:
            perc_outside_genes = round( (add_annot_stats_dic["c_outside_genes"] / c_in_regions) * 100, 1)
        perc_inside_prom = 0.0
        if add_annot_stats_dic["c_inside_prom"] > 0:
            perc_inside_prom = round( (add_annot_stats_dic["c_inside_prom"] / c_in_regions) * 100, 1)

        mdtext += '<table style="max-width: 600px; width: 100%; border-collapse: collapse; line-height: 0.8;">' + "\n"
        mdtext += "<thead>\n"
        mdtext += "<tr>\n"
        mdtext += "<th>Region annotation</th>\n"
        mdtext += "<th># input regions</th>\n"
        mdtext += "<th>% input regions</th>\n"
        mdtext += "</tr>\n"
        mdtext += "</thead>\n"
        mdtext += "<tbody>\n"

        mdtext += '<tr>' + "\n"
        mdtext += "<td>Outside genes</td>\n"
        mdtext += "<td>" + str(add_annot_stats_dic["c_outside_genes"]) + "</td>\n"
        mdtext += "<td>" + str(perc_outside_genes) + "</td>\n"
        mdtext += '</tr>' + "\n"

        mdtext += '<tr>' + "\n"
        mdtext += "<td>Inside promoters</td>\n"
        mdtext += "<td>" + str(add_annot_stats_dic["c_inside_prom"]) + "</td>\n"
        mdtext += "<td>" + str(perc_inside_prom) + "</td>\n"
        mdtext += '</tr>' + "\n"

        if args.add_annot_bed:
            perc_add_annot = 0.0
            if add_annot_stats_dic["c_add_annot"] > 0:
                perc_add_annot = round( (add_annot_stats_dic["c_add_annot"] / c_in_regions) * 100, 1)
            mdtext += '<tr>' + "\n"
            mdtext += "<td>" + args.add_annot_id + "</td>\n"
            mdtext += "<td>" + str(add_annot_stats_dic["c_add_annot"]) + "</td>\n"
            mdtext += "<td>" + str(perc_add_annot) + "</td>\n"
            mdtext += '</tr>' + "\n"

        mdtext += '</tbody>' + "\n"
        mdtext += '</table>' + "\n"

        mdtext += "\n&nbsp;\n&nbsp;\n"
        mdtext += "\nColumn IDs have the following meanings: "
        mdtext += "**Region annotation** -> Type of genomic regions overlapping or not overlapping (depending on specified setting) with input regions, "
        mdtext += "**# input regions** -> number of input regions overlapping or not overlapping with region annotations, "
        mdtext += '**% input regions** -> percentage of input regions overlapping or not overlapping with region annotations.'  + "\n"
        mdtext += "\n&nbsp;\n"

    """
    mRNA region coverage plot.

    """

    if mrna_prof_dic:

        mdtext += """
## mRNA region coverage profile ### {#mrna-prof-plot}

"""
        mrna_prof_plot =  "mRNA_region_plot.png"
        mrna_prof_plot_out = plots_out_folder + "/" + mrna_prof_plot

        dataset_id = "rbpbench_search"
        plot_id = "Input region coverage"
        mrna_reg_occ_dic = {}
        mrna_reg_occ_dic[plot_id] = {}
        mrna_reg_occ_dic[plot_id]["5'UTR"] = mrna_prof_dic[dataset_id].utr5_pc_list
        mrna_reg_occ_dic[plot_id]["CDS"] = mrna_prof_dic[dataset_id].cds_pc_list
        mrna_reg_occ_dic[plot_id]["3'UTR"] = mrna_prof_dic[dataset_id].utr3_pc_list
        c_ol_sites = mrna_prof_dic[dataset_id].c_ol_sites
        c_all_sites = mrna_prof_dic[dataset_id].c_all_sites
        utr5_len_norm = mrna_prof_dic[dataset_id].utr5_len_norm
        cds_len_norm = mrna_prof_dic[dataset_id].cds_len_norm
        utr3_len_norm = mrna_prof_dic[dataset_id].utr3_len_norm
        norm_mode = mrna_prof_dic[dataset_id].norm_mode
        c_ol_mrnas = mrna_prof_dic[dataset_id].c_ol_mrnas
        perc_ol_sites = 0.0
        if c_ol_sites and c_all_sites:
            perc_ol_sites = round(c_ol_sites / c_all_sites * 100, 1)
        perc_min_overlap = round(args.gtf_min_mrna_overlap * 100, 1)

        create_mrna_region_occ_plot([plot_id], mrna_reg_occ_dic,
                                    annot2color_dic, mrna_prof_plot_out,
                                    plot_pdf=args.plot_pdf,
                                    rbp_id=False)

        plots_path = plots_folder + "/" + mrna_prof_plot

        mdtext += '<img src="' + plots_path + '" alt="mRNA region coverage plot"' + "\n"
        mdtext += 'title="mRNA region coverage plot" />' + "\n"
        mdtext += """
**Figure:** mRNA region coverage profile for provided input regions (# input regions = %i, # of input regions overlapping with mRNAs = %i). 
Percentage of input regions overlapping with mRNA exons = %s%%. 
**NOTE** that if this percentage is low, it likely means that the input regions (if derived from CLIP-seq or similar protocols) 
originate from an intron binding RBP, making the plot less informative (since the RBP is typically not binding to spliced mRNAs).
Minimum overlap amount with mRNA exons required for input region to be counted as overlapping = %s%% (set via --gtf-min-mrna-overlap).
All overlapping input region positions are used for the coverage calculation.
Only mRNA regions overlapping with input regions are used for plot generation (# mRNAs with input regions = %i).
mRNA region lengths used for plotting are derived from the occupied mRNA regions, using their %s region lengths (5'UTR = %s, CDS = %s, 3'UTR = %s).

&nbsp;

""" %(c_all_sites, c_ol_sites, str(perc_ol_sites), str(perc_min_overlap), c_ol_mrnas, norm_mode, str(utr5_len_norm), str(cds_len_norm), str(utr3_len_norm))



    """
    Exon-intron overlap statistics

    """

    if ei_ol_stats_dic:

        mdtext += """
## Exon-intron overlap statistics ### {#ei-ol-stats}

"""

        ei_ol_stats = ei_ol_stats_dic["rbpbench_search"]
        exon_sites_perc = 0.0
        if ei_ol_stats.c_exon_sites and ei_ol_stats.c_input_sites:
            exon_sites_perc = round(ei_ol_stats.c_exon_sites / ei_ol_stats.c_input_sites * 100, 1)
        intron_sites_perc = 0.0
        if ei_ol_stats.c_intron_sites and ei_ol_stats.c_input_sites:
            intron_sites_perc = round(ei_ol_stats.c_intron_sites / ei_ol_stats.c_input_sites * 100, 1)
        us_ib_sites_perc = 0.0
        if ei_ol_stats.c_us_ib_sites and ei_ol_stats.c_input_sites:
            us_ib_sites_perc = round(ei_ol_stats.c_us_ib_sites / ei_ol_stats.c_input_sites * 100, 1)
        ds_ib_sites_perc = 0.0
        if ei_ol_stats.c_ds_ib_sites and ei_ol_stats.c_input_sites:
            ds_ib_sites_perc = round(ei_ol_stats.c_ds_ib_sites / ei_ol_stats.c_input_sites * 100, 1)
        us_ib_dist_sites_perc = 0.0
        if ei_ol_stats.c_us_ib_dist_sites and ei_ol_stats.c_input_sites:
            us_ib_dist_sites_perc = round(ei_ol_stats.c_us_ib_dist_sites / ei_ol_stats.c_input_sites * 100, 1)
        ds_ib_dist_sites_perc = 0.0
        if ei_ol_stats.c_ds_ib_dist_sites and ei_ol_stats.c_input_sites:
            ds_ib_dist_sites_perc = round(ei_ol_stats.c_ds_ib_dist_sites / ei_ol_stats.c_input_sites * 100, 1)
        first_exon_sites_perc = 0.0
        if ei_ol_stats.c_first_exon_sites and ei_ol_stats.c_input_sites:
            first_exon_sites_perc = round(ei_ol_stats.c_first_exon_sites / ei_ol_stats.c_input_sites * 100, 1)
        last_exon_sites_perc = 0.0
        if ei_ol_stats.c_last_exon_sites and ei_ol_stats.c_input_sites:
            last_exon_sites_perc = round(ei_ol_stats.c_last_exon_sites / ei_ol_stats.c_input_sites * 100, 1)
        single_exon_sites_perc = 0.0
        if ei_ol_stats.c_single_exon_sites and ei_ol_stats.c_input_sites:
            single_exon_sites_perc = round(ei_ol_stats.c_single_exon_sites / ei_ol_stats.c_input_sites * 100, 1)
        eib_sites_perc = 0.0
        if ei_ol_stats.c_eib_sites and ei_ol_stats.c_input_sites:
            eib_sites_perc = round(ei_ol_stats.c_eib_sites / ei_ol_stats.c_input_sites * 100, 1)
        perc_min_overlap = round(args.gtf_feat_min_overlap * 100, 1)

        intron_bl =  args.gtf_intron_border_len

        # Make bar plot.
        categories = ['Exon\nregions', 'Intron\nregions', '%i nt us\nintron\nregions' %(intron_bl), '%i nt ds\nintron\nregions' %(intron_bl), '> %i nt\nus intron\nregions' %(intron_bl), '> %i nt\nds intron\nregions' %(intron_bl), '+/- 50 nt\nexon-intron\nborders', 'First exon', 'Last exon', 'Single exon']
        percentages = [exon_sites_perc, intron_sites_perc, us_ib_sites_perc, ds_ib_sites_perc, us_ib_dist_sites_perc, ds_ib_dist_sites_perc, eib_sites_perc, first_exon_sites_perc, last_exon_sites_perc, single_exon_sites_perc]

        fig, ax = plt.subplots(figsize=(10.5, 4))

        ax.bar(categories, percentages, color='#e5ecf6', zorder=2)

        ax.set_ylabel('Input regions overlapping with %')

        ax.set_yticks(range(0, 101, 10))

        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['bottom'].set_visible(False)

        ax.grid(axis='y', color='#e5ecf6', linestyle='-', alpha=0.7, linewidth=0.7, zorder=1)  # #e5ecf6 lightgray

        plt.subplots_adjust(left=0.06, right=0.99, top=0.95, bottom=0.15)

        eib_stats_plot = "eib_ol_stats.bar_plot.png"
        eib_stats_plot_out = plots_out_folder + "/" + eib_stats_plot
        plot_path = plots_folder + "/" + eib_stats_plot

        plt.savefig(eib_stats_plot_out, dpi=150)

        if args.plot_pdf and eib_stats_plot_out.endswith('.png'):
            pdf_out = eib_stats_plot_out[:-4] + '.pdf'
            plt.savefig(pdf_out, dpi=150)

        plt.close()
        mdtext += '<image src = "' + plot_path + '" width="1200px"></image>'  + "\n"
        # mdtext += '<img src="' + plots_path + '" alt="Exon intron overlap plot"' + "\n"
        # mdtext += 'title="mRNA region occupancy plot" />' + "\n"


        mdtext += r"""
**Figure:** Exon, intron + border region overlap statistics. \# input regions = %i. 
\# input regions overlapping with exon regions = %i.
\# input regions overlapping with intron regions = %i.
Categories:
**Exon regions** -> %% of input regions overlapping with exon regions (== exonic regions).
**Intron regions** -> %% of input regions overlapping with intron regions (== intronic regions).
**%i nt us intron regions** -> %% of intronic regions closer to intron upstream borders + within first %i nt of intron (\# %i).
**%i nt ds intron regions** -> %% of intronic regions closer to intron downstream borders + within first %i nt of intron (\# %i).
**> %i nt us intron regions** -> %% of intronic regions closer to intron upstream borders and > %i nt away from intron borders (\# %i).
**> %i nt ds intron regions** -> %% of intronic regions closer to intron downstream borders and > %i nt away from intron borders (\# %i).
**+/- 50 nt exon intron borders** -> %% of input regions overlapping with exon-intron borders 
(50 nt upstream and downstream of exon-intron borders) (\# %i).
**Exon regions** -> %% of input regions overlapping with exon regions (== exonic regions).
**First exon** -> %% of input regions overlapping with transcript first exons (\# %i).
**Last exon** -> %% of input regions overlapping with transcript last exons (\# %i).
**Single exon** -> %% of input regions overlapping with single exon transcripts (\# %i).
**NOTE** that for upstream/downstream intron end overlaps, the input region is always assigned 
to the closest intron border (based on distance between input region center position and intron ends), 
independent of intron length (i.e., intron length can also be < %i nt).

&nbsp;

""" %(ei_ol_stats.c_input_sites, ei_ol_stats.c_exon_sites, ei_ol_stats.c_intron_sites, intron_bl, intron_bl, ei_ol_stats.c_us_ib_sites, intron_bl, intron_bl, ei_ol_stats.c_ds_ib_sites, intron_bl, intron_bl, ei_ol_stats.c_us_ib_dist_sites, intron_bl, intron_bl, ei_ol_stats.c_ds_ib_dist_sites, ei_ol_stats.c_eib_sites, ei_ol_stats.c_first_exon_sites, ei_ol_stats.c_last_exon_sites, ei_ol_stats.c_single_exon_sites, intron_bl)


    """
    Intronic regions intron border distance plots.

    """

    if reg2annot_dic:

        assert reg2seq_dic, "reg2seq_dic for intronic regions intron border distance plots"

        mdtext += """
## Intronic regions intron border distances ### {#intron-border-dist}

"""
        # reg2annot_dic format: reg_id -> [annot_id, tr_id, border_dist, us_ds_label, annot_reg_len, exon_intron_nr, gene_id, gene_name, tr_biotype]
        intron_regions = False
        for reg_id in reg2annot_dic:
            annot_id = reg2annot_dic[reg_id][0]
            if annot_id == "intron":
                intron_regions = True
                break
                
        c_us_regions = 0
        c_ds_regions = 0  
        for reg_id in reg2annot_dic:
            annot_id = reg2annot_dic[reg_id][0]
            if annot_id == "intron":
                us_ds_label = reg2annot_dic[reg_id][3]
                if us_ds_label == "up":
                    c_us_regions += 1
                elif us_ds_label == "down":
                    c_ds_regions += 1

        c_all_regions = c_us_regions + c_ds_regions
        perc_us_regions = 0.0
        perc_ds_regions = 0.0
        if c_all_regions:
            perc_us_regions = round(c_us_regions / c_all_regions * 100, 1)
            perc_ds_regions = round(c_ds_regions / c_all_regions * 100, 1)

        if intron_regions:

            us_intron_border_dist_plot_plotly =  "us_intron_border_dist_plot.plotly.html"
            us_intron_border_dist_plot_plotly_out = plots_out_folder + "/" + us_intron_border_dist_plot_plotly
            ds_intron_border_dist_plot_plotly =  "ds_intron_border_dist_plot.plotly.html"
            ds_intron_border_dist_plot_plotly_out = plots_out_folder + "/" + ds_intron_border_dist_plot_plotly

            create_intron_border_dist_plot_plotly(reg2annot_dic, reg2seq_dic,
                                                  us_intron_border_dist_plot_plotly_out,
                                                  ds_intron_border_dist_plot_plotly_out,
                                                  reg2motifs_dic=reg2motifs_dic,
                                                  include_plotlyjs=include_plotlyjs,
                                                  full_html=plotly_full_html)


            us_plot_path = plots_folder + "/" + us_intron_border_dist_plot_plotly
            ds_plot_path = plots_folder + "/" + ds_intron_border_dist_plot_plotly

            # Upstream plot.
            if args.plotly_js_mode in [5, 6, 7]:
                js_code = read_file_content_into_str_var(us_intron_border_dist_plot_plotly_out)
                js_code = js_code.replace("height:100%; width:100%;", "height:1000px; width:1000px;")
                mdtext += js_code + "\n"
            else:
                if plotly_embed_style == 1:
                    mdtext += "<div>\n"
                    mdtext += '<iframe src="' + us_plot_path + '" width="1200" height="800"></iframe>' + "\n"
                    mdtext += '</div>'
                elif plotly_embed_style == 2:
                    mdtext += '<iframe src="' + us_plot_path + '" width="1200" height="800"></iframe>' + "\n"


            mdtext += """

**Figure:** Distances of intronic input regions to their next intron border. Only input regions closer to the upstream 
border are shown (# regions = %i, percentage =  %s%%). The input region center position is taken to calculate the distance.

&nbsp;

""" %(c_us_regions, str(perc_us_regions))

            # Downstream plot.
            if args.plotly_js_mode in [5, 6, 7]:
                js_code = read_file_content_into_str_var(ds_intron_border_dist_plot_plotly_out)
                js_code = js_code.replace("height:100%; width:100%;", "height:800px; width:1200px;")
                mdtext += js_code + "\n"
            else:
                if plotly_embed_style == 1:
                    mdtext += "<div>\n"
                    mdtext += '<iframe src="' + ds_plot_path + '" width="1200" height="800"></iframe>' + "\n"
                    mdtext += '</div>'
                elif plotly_embed_style == 2:
                    mdtext += '<iframe src="' + ds_plot_path + '" width="1200" height="800"></iframe>' + "\n"

            mdtext += """

**Figure:** Distances of intronic input regions to their next intron border. Only input regions closer to the downstream 
border are shown (# regions = %i, percentage =  %s%%). The input region center position is taken to calculate the distance.

&nbsp;

""" %(c_ds_regions, str(perc_ds_regions))


        else:
            mdtext += """

No intronic regions intron border distance plot generated since there are no input regions overlapping with introns.

&nbsp;

"""



    """
    RBP region occupancies upset plot.

    """

    if args.enable_upset_plot:

        mdtext += """
## RBP combinations upset plot ### {#rbp-comb-upset-plot}

"""

        rbp_reg_occ_upset_plot =  "rbp_region_occupancies.upset_plot.png"
        rbp_reg_occ_upset_plot_out = plots_out_folder + "/" + rbp_reg_occ_upset_plot

        upset_plot_nr_included_rbps = len(rbp2regidx_dic)
        if args.upset_plot_max_rbp_rank is not None:
            if args.upset_plot_max_rbp_rank < upset_plot_nr_included_rbps:
                upset_plot_nr_included_rbps = args.upset_plot_max_rbp_rank

        plotted, reason, count = create_rbp_reg_occ_upset_plot(rbp2regidx_dic, reg_ids_list, 
                                    reg2annot_dic=reg2annot_dic,
                                    min_degree=args.upset_plot_min_degree,
                                    max_degree=args.upset_plot_max_degree,
                                    min_subset_size=args.upset_plot_min_subset_size,
                                    max_subset_rank=args.upset_plot_max_subset_rank,
                                    min_rbp_count=args.upset_plot_min_rbp_count,
                                    max_rbp_rank=args.upset_plot_max_rbp_rank,
                                    plot_pdf=args.plot_pdf,
                                    plot_out=rbp_reg_occ_upset_plot_out)


        plot_path = plots_folder + "/" + rbp_reg_occ_upset_plot

        if plotted:
            mdtext += '<img src="' + plot_path + '" alt="RBP region occupancies upset plot"' + "\n"
            mdtext += 'title="RBP region occupancies upset plot" />' + "\n"
            mdtext += """

**Figure:** Upset plot of RBP combinations found in the given set of %s regions (# of regions = %i). 
Intersection size == how often a specific RBP combination is found in the regions dataset.
For example, if two regions in the input set contain motif hits for RBP1, RBP3, and RBP5, then the RBP combination "RBP1,RBP3,RBP5" will get a count (i.e., Intersection size) of 2.
Minimum occurrence number for a combination to be reported = %i (command line parameter: --upset-plot-min-subset-size). 
How many RBPs a combination must at least contain to be reported = %i (command line parameter: --upset-plot-min-degree).
Maximum rank of a combination (w.r.t. to its intersection size) to be reported = %i (command line parameter: --upset-plot-max-subset-rank).
Number of top RBPs included in statistic + plot (ranked by # input regions with hits) = %i (command line parameter: --upset-plot-max-rbp-rank).
To be included in the statistic + plot, RBPs need to have motif hits in >= %i input regions (command line parameter: --upset-plot-min-rbp-count).
The numbers on the left side for each RBP tell how many %s regions have motif hits of the respective RBP. 
If a GTF file was given, bar charts become stacked bar charts, showing what GTF annotations the regions overlap with (see legend for region types).
NOTE that upsetplot currently (v0.9) only supports distinct mode (no intersect or union modes available). 
In distinct mode, the reported intersection/subset size for a combination is the number of times 
this distinct RBP combination shows up in the data. For example, if the combination "RBP1,RBP3,RBP5" has a count of 2,
it means that there are two regions in the input set which contain motif hits exclusively by RBP1, RBP2, and RBP5 (and no other RBPs!).
This is why the more RBPs are selected, the smaller intersection sizes typically become for specific combinations (down to counts of 1).

&nbsp;

""" %(site_type, c_regions, args.upset_plot_min_subset_size, args.upset_plot_min_degree, args.upset_plot_max_subset_rank, upset_plot_nr_included_rbps, args.upset_plot_min_rbp_count, site_type)

        else:

            if reason == "min_degree":

                mdtext += """

No upset plot generated since set --upset-plot-min-degree > maximum degree found in the RBP combination set. Please use lower number for --upset-plot-min-degree parameter.
Also NOTE that upsetplot currently (v0.9) only supports distinct mode (no intersect or union modes available). 
In distinct mode, the reported intersection/subset size for a combination is the number of times 
this distinct RBP combination shows up in the data. For example, if the combination "RBP1,RBP3,RBP5" has a count of 2,
it means that there are two regions in the input set which contain motif hits exclusively by RBP1, RBP2, and RBP5 (and no other RBPs!).
This is why the more RBPs are selected, the smaller intersection sizes typically become for specific combinations (down to counts of 1).

&nbsp;

"""

            elif reason == "min_subset_size":

                mdtext += """

No upset plot generated since set --upset-plot-min-subset-size (%i) > maximum subset size (%i) found in the RBP combination set. Please use lower number for --upset-plot-min-subset-size parameter.
Also NOTE that upsetplot currently (v0.9) only supports distinct mode (no intersect or union modes available). 
In distinct mode, the reported intersection/subset size for a combination is the number of times 
this distinct RBP combination shows up in the data. For example, if the combination "RBP1,RBP3,RBP5" has a count of 2,
it means that there are two regions in the input set which contain motif hits exclusively by RBP1, RBP2, and RBP5 (and no other RBPs!).
This is why the more RBPs are selected, the smaller intersection sizes typically become for specific combinations (down to counts of 1).

&nbsp;

""" %(args.upset_plot_min_subset_size, count)

            elif reason == "len(rbp_id_list) == 1":

                mdtext += """

No plot generated since number of selected RBPs == 1.

&nbsp;

"""

            elif reason == "no_region_hits":

                mdtext += """

No plot generated since no motif hits found in input regions.

&nbsp;

"""

            elif reason == "min_rbp_count":

                mdtext += """

No plot generated since set --upset-plot-min-rbp-count results in no RBPs remaining for upset plot.

&nbsp;

"""

            else:
                assert False, "invalid reason given for not plotting upset plot"



    """
    Set RBP motif distance stats + plot.

    """

    # If --set-rbp-id and no motif hits, disable set_rbp_id again.
    set_rbp_id = args.set_rbp_id
    if args.set_rbp_id is not None and no_region_hits:
        mdtext += """
## Set RBP %s motifs distance statistics ### {#rbp-motif-dist-stats}

No motif distance statistics and plots generated since no motif hits found in input regions.


""" %(args.set_rbp_id)
        
        set_rbp_id = None

    if set_rbp_id is not None:

        rbp_motif_dist_plot_plotly =  "%s.rbp_level_dist.plotly.html" %(set_rbp_id)
        rbp_motif_dist_plot_plotly_out = plots_out_folder + "/" + rbp_motif_dist_plot_plotly

        plotted, pc_dic, in_dic, out_dic = plot_motif_dist_rbp_level(set_rbp_id,
            region_rbp_motif_pos_dic,
            id2name_dic,
            name2ids_dic,
            reg2pol_dic,
            html_out=rbp_motif_dist_plot_plotly_out,
            include_plotlyjs=include_plotlyjs,
            full_html=plotly_full_html,
            line_plot_range=args.motif_distance_plot_range,
            min_pair_count=args.rbp_min_pair_count)

        mdtext += """
## Set RBP %s motifs distance statistics ### {#rbp-motif-dist-stats}

**Table:** Motif distance statistics between set RBP ID (%s) motifs and other RBP ID motifs (including itself). 
Note that the statistics are generated by focussing (i.e., centering) on the highest-scoring motif of the set RBP 
(lowest p-value for sequence motifs, highest bit score for structure motifs) in each input region. 
If a regex is included in the analysis (--regex), the first regex hit in each region is used for the statistics 
(as all regex hits have same score).
In case of an empty table, try to lower --rbp-min-pair-count (current value: %i).

""" %(set_rbp_id, set_rbp_id, args.rbp_min_pair_count)

        # mdtext += '| Set RBP ID | Other RBP ID | Pair count | # near motifs | # distant motifs |' + " \n"
        # mdtext += "| :-: | :-: | :-: | :-: | :-: |\n"

        mdtext += '<table style="max-width: 800px; width: 100%; border-collapse: collapse; line-height: 0.8;">' + "\n"
        mdtext += "<thead>\n"
        mdtext += "<tr>\n"
        mdtext += "<th>Set RBP ID</th>\n"
        mdtext += "<th>Other RBP ID</th>\n"
        mdtext += "<th>Pair count</th>\n"
        mdtext += "<th># near motifs</th>\n"
        mdtext += "<th># distant motifs</th>\n"
        mdtext += "</tr>\n"
        mdtext += "</thead>\n"
        mdtext += "<tbody>\n"

        top_c = 0
        for rbp_id, pair_c in sorted(pc_dic.items(), key=lambda item: item[1], reverse=True):
            if pair_c >= args.rbp_min_pair_count:
                top_c += 1
                # mdtext += "| %s | %s | %i | %i | %i |\n" %(set_rbp_id, rbp_id, pair_c, in_dic[rbp_id], out_dic[rbp_id])
                mdtext += '<tr>' + "\n"
                mdtext += "<td>" + set_rbp_id + "</td>\n"
                mdtext += "<td>" + rbp_id + "</td>\n"
                mdtext += "<td>" + str(pair_c) + "</td>\n"
                mdtext += "<td>" + str(in_dic[rbp_id]) + "</td>\n"
                mdtext += "<td>" + str(out_dic[rbp_id]) + "</td>\n"
                mdtext += '</tr>' + "\n"

        mdtext += '</tbody>' + "\n"
        mdtext += '</table>' + "\n"

        mdtext += "\n&nbsp;\n&nbsp;\n"
        mdtext += "\nColumn IDs have the following meanings: "
        mdtext += "**Set RBP ID** -> set RBP ID (specified via --set-rbp-id), "
        mdtext += "**Other RBP ID** -> other RBP ID that set RBP ID is compared to, "
        mdtext += "**Pair count** -> number of input regions with motif hits from both RBP IDs (minimum pair count to be reported: %i (set by --rbp-min-pair-count)), " %(args.rbp_min_pair_count)
        mdtext += '**# near motifs** -> ' + "number of other RBP ID motifs within specified distance (--motif-distance-plot-range %i) of the centered set RBP ID motifs, " %(args.motif_distance_plot_range)
        mdtext += '**# distant motifs** -> ' + "number of other RBP ID motifs outside specified distance (--motif-distance-plot-range %i) of the centered set RBP ID motifs." %(args.motif_distance_plot_range)
        # mdtext += "\n"
        mdtext += "\n&nbsp;\n"

#         mdtext += """
# ## %s motif distance plot ### {#rbp-motif-dist-plot}

# """ %(set_rbp_id)

        if plotted:

            if args.plotly_js_mode in [5, 6, 7]:
                js_code = read_file_content_into_str_var(rbp_motif_dist_plot_plotly_out)
                js_code = js_code.replace("height:100%; width:100%;", "height:700px; width:1200px;")
                mdtext += js_code + "\n"

            else:

                plot_path = plots_folder + "/" + rbp_motif_dist_plot_plotly

                if plotly_embed_style == 1:
                    mdtext += "<div>\n"
                    mdtext += '<iframe src="' + plot_path + '" width="1200" height="700"></iframe>' + "\n"
                    mdtext += '</div>'
                elif plotly_embed_style == 2:
                    mdtext += '<object data="' + plot_path + '" width="1200" height="700"> </object>' + "\n"

            mdtext += """

**Figure:** Line plot showing motif distances between set RBP ID (%s) motifs (using highest-scoring %s motif for each input region) and other RBP ID motifs (including itself).
Coverage corresponds to positions occupied by RBP motifs. 
The plot is centered on the highest-scoring %s motif (lowest p-value for sequence or highest bit score) for each region containing >= 1 %s motif.
Each RBP with a pair count (definition see table above) of >= %i is shown, and the coverage of the RBP is the accumulation of its individual motif coverages.

&nbsp;

""" %(set_rbp_id, set_rbp_id, set_rbp_id, set_rbp_id, args.rbp_min_pair_count)

        else:

            mdtext += """

<em>No motif distance plot generated for set RBP %s. Try to lower --rbp-min-pair-count (current value: %i). Other reasons: only one motif hit found in total or no regions with > 1 motif.</em>

&nbsp;

""" %(set_rbp_id, args.rbp_min_pair_count)

        """
        Set RBP single motif distance stats + plots.
        
        """

        for idx, motif_id in enumerate(name2ids_dic[set_rbp_id]):

            motif_id_plot_str = motif_id
            if set_rbp_id == args.regex_id:
                motif_id_plot_str = args.regex_id

            single_motif_dist_plot_plotly =  "%s.motif_level_dist.plotly.html" %(motif_id_plot_str)
            single_motif_dist_plot_plotly_out = plots_out_folder + "/" + single_motif_dist_plot_plotly

            plotted, pc_dic, in_dic, out_dic = plot_motif_dist_motif_level(motif_id,
                region_rbp_motif_pos_dic,
                name2ids_dic,
                reg2pol_dic,
                html_out=single_motif_dist_plot_plotly_out,
                include_plotlyjs=include_plotlyjs,
                full_html=plotly_full_html,
                line_plot_range=args.motif_distance_plot_range,
                min_pair_count=args.motif_min_pair_count)


            mdtext += """
### Motif %s distance statistics ### {#single-motif-%i-dist-stats}

""" %(motif_id, idx)

            # Plot motif (if sequence motif).
            if motif_id in seq_motif_blocks_dic:

                motif_plot = "%s.%s.png" %(set_rbp_id, motif_id)
                motif_plot_out = plots_out_folder + "/" + motif_plot
                plot_path = plots_folder + "/" + motif_plot

                # Check if motif in motif database folder.
                if args.motif_db_str:
                    db_motif_path = benchlib_path + "/content/motif_plots/%s" %(motif_plot)
                    if os.path.exists(db_motif_path):
                        shutil.copy(db_motif_path, motif_plot_out)
                        if args.plot_pdf:
                            create_motif_plot(motif_id, seq_motif_blocks_dic,
                                            motif_plot_out,
                                            plot_pdf=True,
                                            plot_png=False)

                if not os.path.exists(motif_plot_out):
                    create_motif_plot(motif_id, seq_motif_blocks_dic,
                                      motif_plot_out,
                                      plot_pdf=args.plot_pdf,
                                      plot_png=True)

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
If a regex is included in the analysis (--regex), the first regex hit in each region is used for the statistics 
(as all regex hits have same score).
In case of an empty table, try to lower --motif-min-pair-count (current value: %i).

""" %(motif_id, set_rbp_id, motif_id, args.motif_min_pair_count)

            # mdtext += '| Motif ID | Other motif ID | &nbsp; Other motif ID plot &nbsp; | Pair count | # near motifs | # distant motifs |' + " \n"
            # mdtext += "| :-: | :-: | :-: | :-: | :-: | :-: |\n"

            mdtext += '<table style="max-width: 1000px; width: 100%; border-collapse: collapse; line-height: 0.8;">' + "\n"
            mdtext += "<thead>\n"
            mdtext += "<tr>\n"
            mdtext += "<th>Motif ID</th>\n"
            mdtext += "<th>Other motif ID</th>\n"
            mdtext += "<th>Other motif ID plot</th>\n"
            mdtext += "<th>Pair count</th>\n"
            mdtext += "<th># near motifs</th>\n"
            mdtext += "<th># distant motifs</th>\n"
            mdtext += "</tr>\n"
            mdtext += "</thead>\n"
            mdtext += "<tbody>\n"

            for other_motif_id, pair_c in sorted(pc_dic.items(), key=lambda item: item[1], reverse=True):

                if pair_c >= args.motif_min_pair_count:

                    # Plot motif (if sequence motif).
                    plot_str = "-"
                    if other_motif_id in seq_motif_blocks_dic:

                        other_rbp_id = id2name_dic[other_motif_id]
                        motif_plot = "%s.%s.png" %(other_rbp_id, other_motif_id)
                        motif_plot_out = plots_out_folder + "/" + motif_plot
                        plot_path = plots_folder + "/" + motif_plot

                        # Check if motif in motif database folder.
                        if args.motif_db_str:
                            db_motif_path = benchlib_path + "/content/motif_plots/%s" %(motif_plot)
                            if os.path.exists(db_motif_path):
                                shutil.copy(db_motif_path, motif_plot_out)
                                if args.plot_pdf:
                                    create_motif_plot(motif_id, seq_motif_blocks_dic,
                                                    motif_plot_out,
                                                    plot_pdf=True,
                                                    plot_png=False)

                        if not os.path.exists(motif_plot_out):
                            create_motif_plot(other_motif_id, seq_motif_blocks_dic,
                                              motif_plot_out,
                                              plot_pdf=args.plot_pdf,
                                              plot_png=True)

                        plot_str = '<image src = "' + plot_path + '" width="300px"></image>'

                    # else:
                    #     print("Motif ID %s not in seq_motif_blocks_dic ... " %(other_motif_id))

                    # mdtext += "| %s | %s | %s | %i | %i | %i |\n" %(motif_id, other_motif_id, plot_str, pair_c, in_dic[other_motif_id], out_dic[other_motif_id])

                    mdtext += '<tr>' + "\n"
                    mdtext += "<td>" + motif_id + "</td>\n"
                    mdtext += "<td>" + other_motif_id + "</td>\n"
                    mdtext += "<td>" + plot_str + "</td>\n"
                    mdtext += "<td>" + str(pair_c) + "</td>\n"
                    mdtext += "<td>" + str(in_dic[other_motif_id]) + "</td>\n"
                    mdtext += "<td>" + str(out_dic[other_motif_id]) + "</td>\n"
                    mdtext += '</tr>' + "\n"
            
            mdtext += '</tbody>' + "\n"
            mdtext += '</table>' + "\n"

            mdtext += "\n&nbsp;\n&nbsp;\n"
            mdtext += "\nColumn IDs have the following meanings: "
            mdtext += "**Motif ID** -> motif ID (motif belonging to set RBP), "
            mdtext += "**Other motif ID** -> other motif ID that set RBP motif ID is compared to, "
            mdtext += "**Other motif ID plot** -> other motif ID sequence motif plot (if motif is sequence motif), "
            mdtext += "**Pair count** -> number of input regions containing hits for both motifs (minimum pair count to be reported: %i (set by --motif-min-pair-count)), " %(args.motif_min_pair_count)
            mdtext += '**# near motifs** -> ' + "number of other motifs within specified distance (--motif-distance-plot-range %i) of the centered motif belonging to set RBP, " %(args.motif_distance_plot_range)
            mdtext += '**# distant motifs** -> ' + "number of other motifs outside specified distance (--motif-distance-plot-range %i) of the centered motif belonging to set RBP." %(args.motif_distance_plot_range)
            # mdtext += "\n"
            mdtext += "\n&nbsp;\n"

            if plotted:

                if args.plotly_js_mode in [5, 6, 7]:
                    js_code = read_file_content_into_str_var(single_motif_dist_plot_plotly_out)
                    js_code = js_code.replace("height:100%; width:100%;", "height:700px; width:1200px;")
                    mdtext += js_code + "\n"

                else:
                    plot_path = plots_folder + "/" + single_motif_dist_plot_plotly

                    if plotly_embed_style == 1:
                        mdtext += '<div class=class="container-fluid" style="margin-top:40px">' + "\n"
                        mdtext += '<iframe src="' + plot_path + '" width="1200" height="700"></iframe>' + "\n"
                        mdtext += '</div>'
                    elif plotly_embed_style == 2:
                        mdtext += '<object data="' + plot_path + '" width="1200" height="700"> </object>' + "\n"

                mdtext += """

**Figure:** Line plot showing motif distances between motif %s (using highest-scoring hit of motif %s for each input region) and other motifs (from same and different RBPs).
Coverage corresponds to positions occupied by motifs.
The plot is centered on the highest-scoring hit of motif %s (lowest p-value for sequence or highest bit score) for each region containing >= 1 hit of motif %s.
Only motifs with a pair count of >= %i appear in the plot.

&nbsp;

""" %(motif_id, motif_id, motif_id, motif_id, args.motif_min_pair_count)


            else:

                mdtext += """

<em>No motif distance plot generated for motif %s. Try to lower --motif-min-pair-count (current value: %i). Other reasons: only one motif hit found in total or no regions with > 1 motif.</em>

&nbsp;

""" %(motif_id, args.motif_min_pair_count)



    """
    GOA results.

    """

    if args.run_goa:

        mdtext += """
## GO enrichment analysis results ### {#goa-results}

"""
        c_goa_results = 0
        if isinstance(goa_results_df, pd.DataFrame) and not goa_results_df.empty:
            c_goa_results = len(goa_results_df)

        filter_purified_info = " GO terms with significantly higher and lower concentration ([e,p]) in study group are shown."
        filter_purified_info2 = "significant"
        if args.goa_filter_purified:
            filter_purified_info = " Only GO terms with significantly higher concentration in study group are shown."
            filter_purified_info2 = "significantly enriched"
            c_goa_results = len(goa_results_df[goa_results_df["enrichment"] == "e"])
        filter_further_info = ""
        if args.goa_max_child is not None: 
            filter_further_info += " Only GO terms with <= %i children are shown." %(args.goa_max_child)
        if args.goa_min_depth is not None:
            filter_further_info += " Only GO terms with >= %i depth are shown." %(args.goa_min_depth)
        if filter_further_info:
            filter_further_info += " Note that additional filters (children + depth) can result in an empty table. For all significant GO terms (i.e., unfiltered results) check *goa_results.tsv* output table."
        filter_only_cooc_info = ""
        if args.goa_cooc_mode == 2:
            filter_only_cooc_info = " Only target genes are considered which contain regions with motif hits from any specified RBP (including regex)."
        elif args.goa_cooc_mode == 3:
            filter_only_cooc_info = " Only target genes are considered which contain regions with motif hits from all specified RBPs (including regex)."

        if c_goa_results > 0:

            mdtext += """
**Table:** GO enrichment analysis results. # of %s GO terms found: %i. Filter p-value threshold (on corrected p-value) = %s. # of target genes used for GOA: %i. # of background genes used for GOA: %i.
%s %s %s

""" %(filter_purified_info2, c_goa_results, str(goa_stats_dic["pval_thr"]), goa_stats_dic["c_target_genes_goa"], goa_stats_dic["c_background_genes_goa"], filter_only_cooc_info, filter_purified_info, filter_further_info)

            mdtext += '<table style="max-width: 1200px; width: 100%; border-collapse: collapse; line-height: 0.9;">' + "\n"
            mdtext += "<thead>\n"
            mdtext += "<tr>\n"
            mdtext += "<th>GO</th>\n"
            mdtext += "<th>Term</th>\n"
            mdtext += "<th>Class</th>\n"
            mdtext += "<th>p-value</th>\n"
            mdtext += "<th>[e,p]</th>\n"
            mdtext += "<th>Depth</th>\n"
            mdtext += "<th># child</th>\n"
            mdtext += "<th># genes</th>\n"
            mdtext += "<th># study</th>\n"
            mdtext += "<th>% genes</th>\n"
            mdtext += "</tr>\n"
            mdtext += "</thead>\n"
            mdtext += "<tbody>\n"

            for index, row in goa_results_df.iterrows():

                go_id = row['GO']
                go_term = row['term']
                go_class = row['class']
                # go_p = row['p']
                go_p_corr = row['p_corr']
                go_enrichment = row['enrichment']
                go_depth = row['depth']
                go_n_genes = row['n_genes']
                go_n_study = row['n_study']
                go_perc_genes = row['perc_genes']
                go_n_children = row['n_children']

                if args.goa_filter_purified:
                    if go_enrichment == "p":
                        continue

                if args.goa_max_child is not None:
                    if go_n_children > args.goa_max_child:
                        continue
                if args.goa_min_depth is not None:
                    if go_depth < args.goa_min_depth:
                        continue

                mdtext += '<tr>' + "\n"
                mdtext += "<td>" + go_id + "</td>\n"
                mdtext += "<td>" + go_term + "</td>\n"
                mdtext += "<td>" + go_class + "</td>\n"
                mdtext += "<td>" + str(go_p_corr) + "</td>\n"
                mdtext += "<td>" + go_enrichment + "</td>\n"
                mdtext += "<td>" + str(go_depth) + "</td>\n"
                mdtext += "<td>" + str(go_n_children) + "</td>\n"
                mdtext += "<td>" + str(go_n_genes) + "</td>\n"
                mdtext += "<td>" + str(go_n_study) + "</td>\n"
                mdtext += "<td>" + str(go_perc_genes) + "</td>\n"
                mdtext += '</tr>' + "\n"

            mdtext += '</tbody>' + "\n"
            mdtext += '</table>' + "\n"
            
            mdtext += "\n&nbsp;\n&nbsp;\n"
            mdtext += "\nColumn IDs have the following meanings: "
            mdtext += "**GO** -> gene ontology (GO) ID, "
            mdtext += "**Term** -> GO term / name, "
            mdtext += "**Class** -> GO term class (biological_process, molecular_function, or cellular_component), "
            mdtext += "**p-value** -> multiple testing corrected (BH) p-value, "
            mdtext += "**[e,p]** -> e: enriched, i.e., GO term with significantly higher concentration, p: purified, GO term with significantly lower concentration), "
            mdtext += "**Depth** -> depth / level of GO term in GO hierarchy (the higher number, the more specific), "
            mdtext += "**# child** -> number of GO term children, "
            mdtext += "**# genes** -> number of genes associated with GO term, "
            mdtext += "**# study** -> number of genes in study (i.e., target genes), "
            mdtext += "**% genes** -> percentage of study genes associated with GO term." + "\n"
            mdtext += "\n&nbsp;\n"

        else:

            if "c_target_genes_goa" in goa_stats_dic:

                mdtext += """

No %s GO terms found given p-value threshold of %s. # of target genes used for GOA: %i. # of background genes used for GOA: %i.

&nbsp;

""" %(filter_purified_info2, str(goa_stats_dic["pval_thr"]), goa_stats_dic["c_target_genes_goa"], goa_stats_dic["c_background_genes_goa"])

            else:

                mdtext += """

No significant GO terms found due to no GO IDs associated with target genes. # of initial target genes (i.e., genes overlapping with --in regions): %i.

&nbsp;

""" %(goa_stats_dic["c_target_genes_pre_filter"])


    # Convert mdtext to html.
    md2html = markdown(mdtext, extensions=['attr_list', 'tables'])

    # OUTMD = open(md_out,"w")
    # OUTMD.write("%s\n" %(mdtext))
    # OUTMD.close()

    html_content = html_head + md2html + html_tail

    OUTHTML = open(html_out,"w")
    OUTHTML.write("%s\n" %(html_content))
    OUTHTML.close()


################################################################################

def plot_motif_dist_rbp_level(set_rbp_id,
                              region_rbp_motif_pos_dic,
                              id2name_dic,
                              rbp2mids_dic,
                              reg2pol_dic,
                              html_out="rbp_motif_distances.plotly.html",
                              include_plotlyjs="cdn",
                              full_html=False,
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
            motif_id, s, e, p = motif_str.split(":")
            if id2name_dic[motif_id] != set_rbp_id:
                continue
            pval = float(p)
            if pval < best_motif_pval:  # for regex with all same p-value (0.0), this means first gets selected.
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
            motif_id, s, e, p = motif_str.split(":")
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
        full_html=full_html,
        include_plotlyjs=include_plotlyjs)

    return True, set_rbp_other_rbp_pair_count_dic, rbp_in_range_c_dic, rbp_out_range_c_dic


################################################################################

def plot_motif_dist_motif_level(set_motif_id,
                              region_rbp_motif_pos_dic,
                              rbp2mids_dic,
                              reg2pol_dic,
                              html_out="motif_distances.plotly.html",
                              include_plotlyjs="cdn",
                              full_html=False,
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
            motif_id, s, e, p = motif_str.split(":")
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
            motif_id, s, e, p = motif_str.split(":")
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
        full_html=full_html,
        include_plotlyjs=include_plotlyjs)

    return True, set_motif_other_motif_pair_count_dic, motif_in_range_c_dic, motif_out_range_c_dic


################################################################################

def create_rbp_reg_occ_upset_plot(rbp2regidx_dic, reg_ids_list,
                                  reg2annot_dic=False,
                                  min_degree=2,
                                  max_degree=None,
                                  min_subset_size=10,
                                  max_subset_rank=25,
                                  min_rbp_count=0,
                                  max_rbp_rank=None,
                                  annot2color_dic=False,
                                  plot_pdf=False,
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

    # Filter rbp2regidx_dic if min_rbp_count or max_rbp_rank set.
    if min_rbp_count > 0 or max_rbp_rank is not None:
        rbp2regidx_dic = filter_rbp2regidx_dic(rbp2regidx_dic,
                            min_rbp_count=min_rbp_count,
                            max_rbp_rank=max_rbp_rank)
        # If all RBPs filtered out due to min_rbp_count > max_count.
        if not rbp2regidx_dic:
            return False, "min_rbp_count", 0

    # Check if there are regions that have motif hits.
    no_region_hits = True
    for rbp_id in rbp2regidx_dic:
        if rbp2regidx_dic[rbp_id]:
            no_region_hits = False 
            break
    if no_region_hits:
        return False, "no_region_hits", 0

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

    if reg2annot_dic:
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

        if not annot2color_dic:

            # Get annotation ID -> hex color dictionary.
            annot2color_dic = {}

            # Get all annotation IDs in dataset.
            annot_dic = {}
            for reg_id in reg2annot_dic:
                annot = reg2annot_dic[reg_id][0]
                if annot not in annot_dic:
                    annot_dic[annot] = 1
                else:
                    annot_dic[annot] += 1


            # hex_colors = get_hex_colors_list(min_len=len(annot_with_hits_dic))
            hex_colors = get_hex_colors_list(min_len=len(annot_dic))

            idx = 0
            for annot in sorted(annot_dic, reverse=False):
                # hc = hex_colors[idx]
                # print("Assigning hex color %s to annotation %s ... " %(hc, annot))
                annot2color_dic[annot] = hex_colors[idx]
                idx += 1

            # idx = 0
            # for annot in sorted(annot_with_hits_dic, reverse=False):
            #     # hc = hex_colors[idx]
            #     # print("Assigning hex color %s to annotation %s ... " %(hc, annot))
            #     annot2color_dic[annot] = hex_colors[idx]
            #     idx += 1
        
        # Set indices (== RBP ID columns).
        df = df.set_index(rbp_id_list)
        # Move on to plotting.
        print("Plotting upset plot with GTF annotations ... ")
        upset = UpSet(df, orientation='horizontal', 
                        min_degree=min_degree,  # number of RBPs in set (e.g. 2 -> at least 2 RBP pairs, not single RBPs)
                        max_degree=max_degree,
                        min_subset_size=min_subset_size,  # min size of a set to be reported.
                        max_subset_rank=max_subset_rank,
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
                      max_subset_rank=max_subset_rank,
                      show_counts=True,
                      sort_by="cardinality")
        upset.plot()

    plt.savefig(plot_out, dpi=125, bbox_inches='tight')

    if plot_pdf and plot_out.endswith('.png'):
        pdf_out = plot_out[:-4] + '.pdf'
        plt.savefig(pdf_out, dpi=125, bbox_inches='tight')

    plt.close()
    return True, "yowza", 0


################################################################################

def create_enmo_annotation_bar_plot(reg2annot_dic, 
                                    annot2color_dic=False,
                                    data_id="Input",
                                    plot_pdf=False,
                                    plot_out="enmo_annotation_bar_plot.png"):
    """
    Plot gene region annotations bar plot for input or backgroud set.

    """
    fheight = 0.8
    fwidth = 10

    # Get all annotation IDs in dataset.
    annot_dic = {}
    for reg_id in reg2annot_dic:
        annot = reg2annot_dic[reg_id][0]
        if annot not in annot_dic:
            annot_dic[annot] = 1
        else:
            annot_dic[annot] += 1

    data_dic = {}
    for annot in sorted(annot_dic, reverse=True):
        data_dic[annot] = [0]

    data_dic["data_id"] = [data_id]

    for reg_id in reg2annot_dic:
        annot = reg2annot_dic[reg_id][0]
        data_dic[annot][0] += 1

    df = pd.DataFrame(data_dic)
    # # Remove annotation columns with no counts.
    # df = df.loc[:, (df != 0).any(axis=0)]

    if not annot2color_dic:
        annot2color_dic = {}
        hex_colors = get_hex_colors_list(min_len=len(annot_dic))
        idx = 0
        for annot in sorted(annot_dic, reverse=False):
            annot2color_dic[annot] = hex_colors[idx]
            idx += 1

    ax = df.set_index('data_id').plot(kind='barh', stacked=True, legend=False, color=annot2color_dic, edgecolor="none", figsize=(fwidth, fheight))

    plt.xlabel('Annotation overlap')
    ax.set_ylabel('')
    ax.yaxis.grid(False)
    ax.xaxis.grid(True)
    ax.set_axisbelow(True)

    # Remove border lines.
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)

    plt.legend(bbox_to_anchor=(1.01, 1), loc='upper left')

    if plot_pdf and plot_out.endswith('.png'):
        pdf_out = plot_out[:-4] + '.pdf'
        plt.savefig(pdf_out, dpi=110, bbox_inches='tight')

    plt.savefig(plot_out, dpi=110, bbox_inches='tight')
    plt.close()


################################################################################

def create_search_annotation_stacked_bars_plot(rbp2regidx_dic, reg_ids_list, reg2annot_dic,
                                               plot_out="annotation_stacked_bars_plot.png",
                                               annot2color_dic=False,
                                               add_all_reg_bar=True,
                                               plot_pdf=False,
                                               all_regions_id="All"):
    """
    Create a stacked bars plot, with each bar showing the annotations for one RBP,
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
    if add_all_reg_bar:
        c_ids += 1
    fheight = 0.6 * c_ids
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

    # Add all regions ID to rbp_id_list.
    if add_all_reg_bar:
        rbp_id_list.append(all_regions_id)

    data_dic = {}
    for annot in sorted(annot_dic, reverse=True):  # make reverse sort to have same color coding as upset plot.
        data_dic[annot] = []
        for rbp_id in rbp_id_list:
            data_dic[annot].append(0)

    # Add region annotation counts for all regions.
    if add_all_reg_bar:
        for annot in annot_dic:
            annot_c = annot_dic[annot]
            data_dic[annot][len(rbp_id_list)-1] = annot_c

    data_dic["rbp_id"] = []

    for rbp_id in rbp_id_list:
        data_dic["rbp_id"].append(rbp_id)

    # annot_with_hits_dic = {}

    for rbp_id in rbp2regidx_dic:
        for hit_idx in rbp2regidx_dic[rbp_id]:
            reg_id = reg_ids_list[hit_idx]
            annot = reg2annot_dic[reg_id][0]
            # annot_with_hits_dic[annot] = 1
            rbp_idx = rbp2idx_dic[rbp_id]
            data_dic[annot][rbp_idx] += 1

    # data_dic:
    # {'intron': [0, 0, 1], 'CDS': [5, 4, 7], "5'UTR": [0, 0, 1], "3'UTR": [0, 2, 2], 'rbp_id': ['RBFOX2', 'PUM2', 'PUM1']}

    df = pd.DataFrame(data_dic)
    # Remove annotation columns with no counts.
    df = df.loc[:, (df != 0).any(axis=0)]

    # print("data_dic:")
    # print(data_dic)
    # print("df:")
    # print(df)

    if not annot2color_dic:
        # Get annotation ID -> hex color mapping for plotting.
        annot2color_dic = {}

        # Colors for all annotations.
        hex_colors = get_hex_colors_list(min_len=len(annot_dic))
        # hex_colors = get_hex_colors_list(min_len=len(annot_with_hits_dic))

        idx = 0
        for annot in sorted(annot_dic, reverse=False):
            annot2color_dic[annot] = hex_colors[idx]
            idx += 1

    # idx = 0
    # for annot in sorted(annot_with_hits_dic, reverse=False):
    #     # hc = hex_colors[idx]
    #     # print("Assigning hex color %s to annotation %s ... " %(hc, annot))
    #     annot2color_dic[annot] = hex_colors[idx]
    #     idx += 1

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

    plt.savefig(plot_out, dpi=110, bbox_inches='tight')

    if plot_pdf and plot_out.endswith('.png'):
        pdf_out = plot_out[:-4] + '.pdf'
        plt.savefig(pdf_out, dpi=110, bbox_inches='tight')

    plt.close()


################################################################################

def create_intron_border_dist_plot_plotly(reg2annot_dic, reg2seq_dic,
                                          us_intron_border_dist_plot_out,
                                          ds_intron_border_dist_plot_out,
                                          reg2motifs_dic=False,
                                          include_plotlyjs="cdn",
                                          full_html=False):
    """
    Create two plotly plots showing distances of intronic regions to upstream
    and downstream intron borders.

    reg2annot_dic format: reg_id -> 
    [annot_id, tr_id, border_dist, us_ds_label, annot_reg_len, exon_intron_nr, gene_id, gene_name, tr_biotype]

    """
    assert reg2annot_dic, "given reg2annot_dic empty"

    import plotly.graph_objects as go

    reg_ids_list = []
    for reg_id in sorted(reg2annot_dic):
        annot_id = reg2annot_dic[reg_id][0]
        if annot_id == "intron":
            reg_ids_list.append(reg_id)

    df = pd.DataFrame(reg_ids_list, columns=["Region ID"])

    df['Transcript ID'] = [reg2annot_dic[reg_id][1] for reg_id in reg_ids_list]
    df['Distance'] = [reg2annot_dic[reg_id][2] for reg_id in reg_ids_list]
    df['Category'] = [reg2annot_dic[reg_id][3] for reg_id in reg_ids_list]  # down or up.
    df['annot_reg_len'] = [reg2annot_dic[reg_id][4] for reg_id in reg_ids_list]
    df['exon_intron_nr'] = [reg2annot_dic[reg_id][5] for reg_id in reg_ids_list]
    df['Gene ID'] = [reg2annot_dic[reg_id][6] for reg_id in reg_ids_list]
    df['Gene name'] = [reg2annot_dic[reg_id][7] for reg_id in reg_ids_list]
    df['Transcript biotype'] = [reg2annot_dic[reg_id][8] for reg_id in reg_ids_list]
    df['Sequence'] = [insert_line_breaks(reg2seq_dic[reg_id], line_len=50) for reg_id in reg_ids_list]
    df['Sequence length'] = [len(reg2seq_dic[reg_id]) for reg_id in reg_ids_list]

    if reg2motifs_dic:
        df['Motif hits'] = [reg2motifs_dic[reg_id] for reg_id in reg_ids_list]

    df_up = df[df['Category'] == 'up']
    df_down = df[df['Category'] == 'down']

    hover_data = ['Region ID', 'Sequence', 'Sequence length', 'Transcript ID', 'Gene ID', 'Gene name', 'Transcript biotype', 'Distance', 'annot_reg_len', 'exon_intron_nr']

    motif_hover_info = ''
    if reg2motifs_dic:
        hover_data += ['Motif hits']
        motif_hover_info = '<br>Motif hits:<br>%{customdata[10]}'

    """
    Create upstream plot.

    """
    fig = go.Figure()

    fig.add_trace(go.Violin(
        x=df_up['Distance'],
        y=df_up['Category'],
        box_visible=False,
        points="all",
        side='positive',  # Show upwards.
        bandwidth=100.0,  # Higher values == smoothing.
        spanmode='hard',  # Ensure density starts at 0.
        orientation='h',  # Make the plot horizontal.
        # name='down',
        jitter=0.1,  # Reduce jitter to move points closer.
        pointpos=-0.3,  # Move points further down.
        line_color='#2ca02c',  # Set the color of the points and the line
        fillcolor='#2ca02c'  # Set the color of the filled area
    ))

    fig.update_layout(
        xaxis_title="Distance to upstream intron border",
        yaxis=dict(
            showticklabels=False  # Remove y-axis labels.
        ),
        # yaxis_title="Category",
        violingap=0.1,  # Reduce gap between violins
        violingroupgap=0.1  # Reduce group gap
    )

    # hover_data:
    # ['Region ID', 'Sequence', 'Sequence length', 'Transcript ID', 'Gene ID', 'Gene name', 
    # 'Transcript biotype', 'Distance', 'annot_reg_len', 'exon_intron_nr', 'Motif hits']

    fig.update_traces(
        width=0.3,
        customdata=df_up[hover_data].values,
        hovertemplate=(
            '<span style="font-family: \'Courier New\', monospace;">>%{customdata[0]}</span><br>'
            '<span style="font-family: \'Courier New\', monospace;">%{customdata[1]}</span><br>'
            'Sequence Length: %{customdata[2]}<br>'
            'Transcript ID: %{customdata[3]}<br>'
            'Gene ID: %{customdata[4]}<br>'
            'Gene name: %{customdata[5]}<br>'
            'Transcript biotype: %{customdata[6]}<br>'
            'Border distance: %{customdata[7]}<br>'
            'Intron length: %{customdata[8]}<br>'
            'Intron number: %{customdata[9]}' + motif_hover_info + '<extra></extra>'
        )
    )

    fig.write_html(us_intron_border_dist_plot_out,
                   full_html=full_html,
                   include_plotlyjs=include_plotlyjs)

    """
    Create downstream plot.

    """

    fig = go.Figure()

    fig.add_trace(go.Violin(
        x=df_down['Distance'],
        y=df_down['Category'],
        box_visible=False,
        points="all",
        side='positive',  # Show upwards.
        bandwidth=100.0,  # Higher values == smoothing.
        spanmode='hard',  # Ensure density starts at 0.
        orientation='h',  # Make the plot horizontal.
        # name='down',
        jitter=0.1,  # Reduce jitter to move points closer.
        pointpos=-0.3,  # Move points further down.
        line_color='#2ca02c',  # Set the color of the points and the line
        fillcolor='#2ca02c'  # Set the color of the filled area
    ))

    fig.update_layout(
        xaxis_title="Distance to downstream intron border",
        yaxis=dict(
            showticklabels=False  # Remove y-axis labels.
        ),
        # yaxis_title="Category",
        xaxis=dict(
            autorange='reversed'  # Reverse the x-axis, for downstream data part.
        ),
        violingap=0.1,  # Reduce gap between violins
        violingroupgap=0.1  # Reduce group gap
    )

    fig.update_traces(
        width=0.3,
        customdata=df_down[hover_data].values,
        hovertemplate=(
            '<span style="font-family: \'Courier New\', monospace;">>%{customdata[0]}</span><br>'
            '<span style="font-family: \'Courier New\', monospace;">%{customdata[1]}</span><br>'
            'Sequence Length: %{customdata[2]}<br>'
            'Transcript ID: %{customdata[3]}<br>'
            'Gene ID: %{customdata[4]}<br>'
            'Gene name: %{customdata[5]}<br>'
            'Transcript biotype: %{customdata[6]}<br>'
            'Border distance: %{customdata[7]}<br>'
            'Intron length: %{customdata[8]}<br>'
            'Intron number: %{customdata[9]}' + motif_hover_info + '<extra></extra>'
        )
    )

    fig.write_html(ds_intron_border_dist_plot_out,
                   full_html=full_html,
                   include_plotlyjs=include_plotlyjs)


################################################################################

def create_seqs_kmer_plot_plotly(reg2seq_dic, reg2kmer_rat_dic, plot_out,
                                 k=3,
                                 annot2color_dic=False,
                                 reg2annot_dic=False,
                                 reg2motifs_dic=False,
                                 reg2ntps_dic=False,
                                 include_plotlyjs="cdn",
                                 full_html=False):
    """
    Create a plotly scatter plot with PCA of sequence k-mer ratios.

    """

    if reg2annot_dic:
        assert annot2color_dic, "reg2annot_dic set but given annot2color_dic empty"
    if annot2color_dic:
        assert reg2annot_dic, "annot2color_dic set but given reg2annot_dic empty"

    reg_ids_list = []
    reg_kmer_rat_ll = []
    for reg_id in sorted(reg2kmer_rat_dic):
        reg_ids_list.append(reg_id)
        reg_kmer_rat_ll.append(reg2kmer_rat_dic[reg_id])

    pca = PCA(n_components=2)
    data_2d_pca = pca.fit_transform(reg_kmer_rat_ll)
    df = pd.DataFrame(data_2d_pca, columns=['PC1', 'PC2'])
    explained_variance = pca.explained_variance_ratio_ * 100

    df['Region ID'] = reg_ids_list
    df['Sequence'] = [insert_line_breaks(reg2seq_dic[reg_id], line_len=50) for reg_id in reg_ids_list]
    df['Sequence length'] = [len(reg2seq_dic[reg_id]) for reg_id in reg_ids_list]

    if reg2motifs_dic:
        df['Motif hits'] = [reg2motifs_dic[reg_id] for reg_id in reg_ids_list]
    if reg2ntps_dic:
        df['Mono-nucleotide percentages'] = [reg2ntps_dic[reg_id] for reg_id in reg_ids_list]
    motif_hover_info = ''
    hover_data = ['Region ID', 'Sequence', 'Sequence length']

    if reg2annot_dic:

        # reg2annot_dic format: 
        # reg_id -> [annot_id, tr_id, border_dist, us_ds_label, annot_reg_len, exon_intron_nr, gene_id, gene_name, tr_biotype]

        df['Annotation'] = [str(reg2annot_dic[reg_id][0]) for reg_id in reg_ids_list]
        df['Transcript ID'] = [reg2annot_dic[reg_id][1] for reg_id in reg_ids_list]
        df['Gene ID'] = [reg2annot_dic[reg_id][6] for reg_id in reg_ids_list]
        df['Gene name'] = [reg2annot_dic[reg_id][7] for reg_id in reg_ids_list]
        df['Transcript biotype'] = [reg2annot_dic[reg_id][8] for reg_id in reg_ids_list]

        hover_data += ['Annotation', 'Transcript ID', 'Gene ID', 'Gene name', 'Transcript biotype']

        if reg2motifs_dic:
            hover_data += ['Motif hits']
            motif_hover_info = '<br>Motif hits:<br>%{customdata[8]}'
        if reg2ntps_dic:
            hover_data += ['Mono-nucleotide percentages']
            if reg2motifs_dic:
                motif_hover_info += '<br>Mono-nucleotide percentages:<br>%{customdata[9]}'
            else:
                motif_hover_info = '<br>Mono-nucleotide percentages:<br>%{customdata[8]}'

        fig = px.scatter(
            df,
            x='PC1',
            y='PC2',
            color='Annotation',
            color_discrete_map=annot2color_dic,
            labels={
                'PC1': f'PC1 ({explained_variance[0]:.2f}% variance)',
                'PC2': f'PC2 ({explained_variance[1]:.2f}% variance)'
            },
            # hover_name='Region ID',
            hover_data=hover_data
        )

        # Adding '<extra></extra>' removes side bars on left/right of hover boxes with annotation labels.
        fig.update_traces(
            # customdata=df[['Region ID', 'Sequence', 'Sequence length', 'Annotation', 'Transcript ID', 'Gene ID', 'Gene name', 'Transcript biotype']].values, 
            hovertemplate=(
                '<span style="font-family: \'Courier New\', monospace;">>%{customdata[0]}</span><br>'
                '<span style="font-family: \'Courier New\', monospace;">%{customdata[1]}</span><br>'
                # '%{customdata[1]}<br>'
                'Sequence Length: %{customdata[2]}<br>'
                'Annotation: %{customdata[3]}<br>'
                'Transcript ID: %{customdata[4]}<br>'
                'Gene ID: %{customdata[5]}<br>'
                'Gene name: %{customdata[6]}<br>'
                'Transcript biotype: %{customdata[7]}' + motif_hover_info + '<extra></extra>'
            )
            # hoverinfo='skip'  # Disable default hover info
        )

    else:

        if reg2motifs_dic:
            hover_data += ['Motif hits']
            motif_hover_info = '<br>Motif hits:<br>%{customdata[3]}'
        if reg2ntps_dic:
            hover_data += ['Mono-nucleotide percentages']
            if reg2motifs_dic:
                motif_hover_info += '<br>Mono-nucleotide percentages:<br>%{customdata[4]}'
            else:
                motif_hover_info = '<br>Mono-nucleotide percentages:<br>%{customdata[3]}'

        dot_col = "#2b7bba"

        fig = px.scatter(
            df,
            x='PC1',
            y='PC2',
            labels={
                'PC1': f'PC1 ({explained_variance[0]:.2f}% variance)',
                'PC2': f'PC2 ({explained_variance[1]:.2f}% variance)'
            },
            hover_data=hover_data,
            color_discrete_sequence=[dot_col]
        )

        fig.update_traces(
            hovertemplate=(
                '>%{customdata[0]}<br>'
                '%{customdata[1]}<br>'
                'Sequence Length: %{customdata[2]}' + motif_hover_info + '<extra></extra>'
            )
        )

    fig.update_traces(marker=dict(size=8, line=dict(width=0.5, color='white')))
    fig.write_html(plot_out, full_html=full_html, include_plotlyjs=include_plotlyjs)


    """

    fig = px.violin(in_df, y='Sequence Length', box=True, points='all')
    # Set customdata and hovertemplate
    fig.update_traces(customdata=in_df[['Sequence ID', 'Sequence', 'Motif hits']].values, hovertemplate='>%{customdata[0]}<br>%{customdata[1]}<br>Sequence Length: %{y}<br>Motif hits:<br>%{customdata[2]}')
    fig.update_layout(xaxis_title='Density', yaxis_title='Sequence Length', violinmode='overlay')
    fig.write_html(plot_out,
                   full_html=full_html,
                   include_plotlyjs=include_plotlyjs)




    fig.update_traces(customdata=df[['Region ID', 'Sequence', 'Sequence length', 'Annotation', 'Transcript ID', 'Gene ID', 'Gene name', 'Transcript biotype']].values, 
                      hovertemplate='>%{customdata[0]}<br>%{customdata[1]}<br>Sequence Length: %{y}<br>Motif hits:<br>%{customdata[2]}')


    fig.update_traces(
        hovertemplate='<b>%{hovertext}</b><br>Annotation percentages: %{customdata[1]}<extra></extra>'
    )

    # Get highest percentage for every dataset.
    highest_perc_annot_list = []
    perc_annot_str_list = []
    for annots_l in annots_ll:
        perc_annot_string = "<br>"
        highest_perc_annot = "-"
        highest_freq = 0.0
        for idx, annot_freq in enumerate(annots_l):
            annot = annot_ids_list[idx]
            if annot_freq > highest_freq:
                highest_freq = annot_freq
                highest_perc_annot = annot
            # Make percentage out of frequency and round to two decimal places.
            annot_perc = round(annot_freq * 100, 2)
            perc_annot_string += annot + ": " + str(annot_perc) + '%<br>'
        highest_perc_annot_list.append(highest_perc_annot)
        perc_annot_str_list.append(perc_annot_string)
    

    pca = PCA(n_components=2)  # Reduce data to 2 dimensions.
    data_2d_pca = pca.fit_transform(annots_ll)
    df = pd.DataFrame(data_2d_pca, columns=['PC1', 'PC2'])
    df['Dataset ID'] = dataset_ids_list
    df['Highest % annotation'] = highest_perc_annot_list
    df['Annotation percentages'] = perc_annot_str_list

    explained_variance = pca.explained_variance_ratio_ * 100

    fig = px.scatter(
        df,  # Use the DataFrame directly
        x='PC1',
        y='PC2',
        color='Highest % annotation',
        color_discrete_map=annot2color_dic,
        # title='2D Visualization with Dataset IDs',
        labels={
            'PC1': f'PC1 ({explained_variance[0]:.2f}% variance)',
            'PC2': f'PC2 ({explained_variance[1]:.2f}% variance)'
        },
        hover_name='Dataset ID',
        hover_data=['Highest % annotation', 'Annotation percentages']
    )

    fig.update_traces(
        hovertemplate='<b>%{hovertext}</b><br>Annotation percentages: %{customdata[1]}<extra></extra>'
    )

    fig.update_traces(marker=dict(size=8, line=dict(width=0.5, color='white')))
    # fig.update_layout(margin=dict(l=0, r=0, b=0, t=0))
    fig.write_html(plot_out, full_html=full_html, include_plotlyjs=include_plotlyjs)

    """


################################################################################

def seqs_dic_get_kmer_ratios(seqs_dic, k,
                             rna=False,
                             report_key_error=False,
                             skip_non_dic_keys=True,
                             add_seq_comp=False,
                             add_gc_skew=False,
                             add_at_skew=False,
                             add_at_content=False,
                             add_1nt_ratios=False,
                             seq_comp_k=1,
                             reg2ntps_dic=None,
                             ntps_dic=False,
                             convert_to_uc=True):
    """
    Given a dictionary with mappings: sequence ID -> sequence,
    return dictionary with mapping:
    sequence ID -> list of k-mer ratios for the sequence.

    ntps_dic:
        If set, calculate average mono-nucleotide percentages of seqs_dic.

    >>> seqs_dic = {'seq1': 'AACGN', 'seq2': 'ACGT', 'seq3': 'tttt'}
    >>> seqs_dic_get_kmer_ratios(seqs_dic, 1)
    {'seq1': [0.5, 0.25, 0.25, 0.0], 'seq2': [0.25, 0.25, 0.25, 0.25], 'seq3': [0.0, 0.0, 0.0, 1.0]}
    >>> seqs_dic = {'seq1': 'AACTT'}
    >>> seqs_dic_get_kmer_ratios(seqs_dic, 2)
    {'seq1': [0.25, 0.25, 0.0, 0.0, 0.0, 0.0, 0.0, 0.25, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.25]}

    """
    # Checks.
    assert seqs_dic, "given dictinary seqs_dic empty"
    assert k, "invalid k given"
    assert k > 0, "invalid k given"

    reg2kmer_rat_dic = {}

    for seq_id in seqs_dic:

        seq = seqs_dic[seq_id]

        if convert_to_uc:
            seq = seq.upper()

        count_dic = get_kmer_dic(k, rna=rna)

        total_c = 0

        # Count k-mers in sequence.
        for i in range(len(seq)-k+1):
            kmer = seq[i:i+k]
            if skip_non_dic_keys:
                if kmer not in count_dic:
                    continue
            else:
                if report_key_error:
                    assert kmer in count_dic, "k-mer \"%s\" not in count_dic" %(kmer)
            if kmer in count_dic:
                count_dic[kmer] += 1
                total_c += 1
        
        # This can happen if sequence length < k.
        if total_c == 0:
            continue

        for kmer in count_dic:
            ratio = count_dic[kmer] / total_c
            count_dic[kmer] = ratio

        # Convert count_dic to list.
        kmer_rat_list = []
        for kmer in count_dic:
            kmer_rat_list.append(count_dic[kmer])

        # Single-nt counts.
        ntc_dic = get_ntc_dic(seq, rna=rna)
        eff_seq_l = 0
        for nt in ntc_dic:
            eff_seq_l += ntc_dic[nt]

        assert eff_seq_l, "effective sequence length is 0 for sequence %s" %(seq_id)

        # region ID -> nucleotide percentages string.
        if reg2ntps_dic is not None:
            ntps_list = []
            for nt in ntc_dic:
                nt_perc = (ntc_dic[nt] / eff_seq_l) * 100
                if ntps_dic:
                    ntps_dic[nt] += nt_perc
                ntps = "%s: %.2f%%" %(nt, nt_perc)
                ntps_list.append(ntps)
            reg2ntps_dic[seq_id] = ", ".join(ntps_list)

        if add_seq_comp:
            seq_comp = calc_seq_entropy(eff_seq_l, ntc_dic,
                                        k=seq_comp_k)
            kmer_rat_list.append(seq_comp)
        
        if add_gc_skew:
            gc_skew = calc_seq_gc_skew_ntc(ntc_dic)
            kmer_rat_list.append(gc_skew)
        
        if add_at_skew:
            at_skew = calc_seq_at_skew_ntc(ntc_dic, rna=rna)
            kmer_rat_list.append(at_skew)

        if add_at_content:
            at_content = calc_seq_at_content_ntc(eff_seq_l, ntc_dic, rna=rna)
            kmer_rat_list.append(at_content)

        if add_1nt_ratios:
            for nt in ntc_dic:
                nt_ratio = ntc_dic[nt] / eff_seq_l
                kmer_rat_list.append(nt_ratio)

        reg2kmer_rat_dic[seq_id] = kmer_rat_list

    if ntps_dic:
        for nt in ntps_dic:
            ntps_dic[nt] = ntps_dic[nt] / len(seqs_dic)

    return reg2kmer_rat_dic


################################################################################

def seqs_dic_calc_entropies(seqs_dic,
                            rna=True,
                            k=1,
                            return_dic=False):
    """
    Given a dictionary of sequences, calculate entropies for each sequence
    and return list of entropy values.

    seqs_dic:
    Dictionary with sequences.
    k:
    k-mer size for entropy calculation.

    rna:
    Use RNA alphabet for counting (uppercase chars only)

    >>> seqs_dic = {'seq1': 'AAAAAAAA', 'seq2': 'AAAACCCC', 'seq3': 'AACCGGUU'}
    >>> seqs_dic_calc_entropies(seqs_dic)
    [0, 0.5, 1.0]

    """
    assert seqs_dic, "given dictionary seqs_dic empty"
    entr_list = []
    if return_dic:
        entr_dic = {}
    for seq_id in seqs_dic:
        seq = seqs_dic[seq_id]
        seq_l = len(seq)
        # Make uppercase (otherwise seq_l not correct).
        seq = seq.upper()
        # Get nt count dic.
        count_dic = {}
        if k == 1:
            count_dic = seq_count_nt_freqs(seq, rna=rna)
        else:
            count_dic = get_kmer_counts_dic(seq, k, rna=rna)

        # Calculate sequence entropy.
        seq_entr = calc_seq_entropy(seq_l, count_dic,
                                    k=k)
        #if seq_entr > 0.5:
        #    print("Entropy: %.2f" %(seq_entr))
        #    print("%s: %s" %(seq_id, seq))
        if return_dic:
            entr_dic[seq_id] = seq_entr
        else:
            entr_list.append(seq_entr)
    if return_dic:
        return entr_dic
    else:
        return entr_list


################################################################################

def get_bint_perc_from_ntc_dic(ntc_dic):
    """
    Given a DNA nucleotide counts ntc_dic with format:
    {'A': 4, 'C': 3, 'G': 2, 'T': 1}
    Get AC, AG, AT, CG, CT, GT percentages.
    4 elements, select 2, no order, no repeated elements:
    Binomial coefficient -> 6

    >>> ntc_dic = {'A': 4, 'C': 3, 'G': 2, 'T': 1}
    >>> get_bint_perc_from_ntc_dic(ntc_dic)
    {'AC': 70.0, 'AG': 60.0, 'AT': 50.0, 'CG': 50.0, 'CT': 40.0, 'GT': 30.0}
    >>> ntc_dic = {'A': 4, 'C': 6, 'G': 0, 'T': 0}
    >>> get_bint_perc_from_ntc_dic(ntc_dic)
    {'AC': 100.0, 'AG': 40.0, 'AT': 40.0, 'CG': 60.0, 'CT': 60.0, 'GT': 0.0}

    """

    assert ntc_dic, "given dictionary ntc_dic empty"
    # Get total number.
    total_n = 0
    for nt in ntc_dic:
        total_n += ntc_dic[nt]
    bint_perc_dic = {}

    perc_ac = ((ntc_dic["A"] + ntc_dic["C"]) / total_n ) * 100
    perc_ag = ((ntc_dic["A"] + ntc_dic["G"]) / total_n ) * 100
    perc_at = ((ntc_dic["A"] + ntc_dic["T"]) / total_n ) * 100
    perc_cg = ((ntc_dic["C"] + ntc_dic["G"]) / total_n ) * 100
    perc_ct = ((ntc_dic["C"] + ntc_dic["T"]) / total_n ) * 100
    perc_gt = ((ntc_dic["G"] + ntc_dic["T"]) / total_n ) * 100

    bint_perc_dic = {}

    bint_perc_dic["AC"] = perc_ac
    bint_perc_dic["AG"] = perc_ag
    bint_perc_dic["AT"] = perc_at
    bint_perc_dic["CG"] = perc_cg
    bint_perc_dic["CT"] = perc_ct
    bint_perc_dic["GT"] = perc_gt

    return bint_perc_dic


################################################################################

def seqs_dic_count_nt_freqs(seqs_dic,
                            rna=False,
                            convert_to_uc=False,
                            count_dic=False):
    """
    Given a dictionary with sequences seqs_dic, count how many times each
    nucleotide is found in all sequences (== get nt frequencies).
    Return nucleotide frequencies count dictionary.

    By default, a DNA dictionary (A,C,G,T) is used, counting only these
    characters (note they are uppercase!).

    rna:
    Instead of DNA dictionary, use RNA dictionary (A,C,G,U) for counting.

    convert_to_uc:
    Convert sequences to uppercase before counting.

    count_dic:
    Supply a custom dictionary for counting only characters in
    this dictionary + adding counts to this dictionary.

    >>> seqs_dic = {'s1': 'AAAA', 's2': 'CCCGGT'}
    >>> seqs_dic_count_nt_freqs(seqs_dic)
    {'A': 4, 'C': 3, 'G': 2, 'T': 1}
    >>> seqs_dic_count_nt_freqs(seqs_dic, rna=True)
    {'A': 4, 'C': 3, 'G': 2, 'U': 0}

    """
    assert seqs_dic, "given dictionary seqs_dic empty"
    if not count_dic:
        count_dic = {'A': 0, 'C': 0, 'G': 0, 'T': 0}
        if rna:
            count_dic = {'A': 0, 'C': 0, 'G': 0, 'U': 0}
    for seq_id in seqs_dic:
        seq = seqs_dic[seq_id]
        if convert_to_uc:
            seq = seq.upper()
        seq_count_nt_freqs(seq, rna=rna, count_dic=count_dic)
    return count_dic


################################################################################

def seq_count_nt_freqs(seq,
                       rna=False,
                       count_dic=False):
    """
    Count nucleotide (character) frequencies in given sequence seq.
    Return count_dic with frequencies.
    If count_dic is given, add count to count_dic.

    rna:
    Instead of DNA dictionary, use RNA dictionary (A,C,G,U) for counting.

    count_dic:
    Supply a custom dictionary for counting only characters in
    this dictionary + adding counts to this dictionary.

    >>> seq = 'AAAACCCGGT'
    >>> seq_count_nt_freqs(seq)
    {'A': 4, 'C': 3, 'G': 2, 'T': 1}
    >>> seq = 'acgtacgt'
    >>> seq_count_nt_freqs(seq)
    {'A': 0, 'C': 0, 'G': 0, 'T': 0}

    """

    assert seq, "given sequence string seq empty"
    if not count_dic:
        count_dic = {'A': 0, 'C': 0, 'G': 0, 'T': 0}
        if rna:
            count_dic = {'A': 0, 'C': 0, 'G': 0, 'U': 0}
    # Conver to list.
    seq_list = list(seq)
    for nt in seq_list:
        if nt in count_dic:
            count_dic[nt] += 1
    return count_dic


################################################################################

def calc_seq_gc_skew(seq):
    """
    Calculate GC skew of a given sequence seq.

    >>> seq = 'AATTGGCC'
    >>> calc_seq_gc_skew(seq)
    0.0
    >>> seq = 'AATTCC'
    >>> calc_seq_gc_skew(seq)
    -1.0

    """
    assert seq, "given sequence string seq empty"
    seq = seq.upper()
    g_c = seq.count('G')
    c_c = seq.count('C')
    if g_c + c_c == 0:
        return 0
    else:
        return (g_c - c_c) / (g_c + c_c)


################################################################################

def calc_seq_at_skew(seq):
    """
    Calculate AT skew of a given sequence seq.

    >>> seq = 'AATTGGCC'
    >>> calc_seq_at_skew(seq)
    0.0
    >>> seq = 'AAGGCC'
    >>> calc_seq_at_skew(seq)
    1.0

    """
    assert seq, "given sequence string seq empty"
    seq = seq.upper()
    a_c = seq.count('A')
    t_c = seq.count('T')
    if a_c + t_c == 0:
        return 0
    else:
        return (a_c - t_c) / (a_c + t_c)



################################################################################

def calc_seq_gc_skew_ntc(ntc_dic):
    """
    Calculate GC skew of a given a nucleotides count dictionary.

    >>> ntc_dic = {'A': 2, 'C': 2, 'G': 2, 'T': 2}
    >>> calc_seq_gc_skew_ntc(ntc_dic)
    0.0
    >>> ntc_dic = {'A': 2, 'C': 2, 'G': 0, 'T': 2}
    >>> calc_seq_gc_skew_ntc(ntc_dic)
    -1.0

    """
    assert ntc_dic, "given ntc_dic empty"

    g_c = ntc_dic["G"]
    c_c = ntc_dic["C"]
    if g_c + c_c == 0:
        return 0
    else:
        return (g_c - c_c) / (g_c + c_c)


################################################################################

def calc_seq_at_skew_ntc(ntc_dic,
                         rna=False):
    """
    Calculate AT skew of a given a nucleotides count dictionary.

    >>> ntc_dic = {'A': 2, 'C': 2, 'G': 2, 'T': 2}
    >>> calc_seq_at_skew_ntc(ntc_dic)
    0.0
    >>> ntc_dic = {'A': 2, 'C': 2, 'G': 2, 'T': 0}
    >>> calc_seq_at_skew_ntc(ntc_dic)
    1.0

    """
    assert ntc_dic, "given ntc_dic empty"

    a_c = ntc_dic["A"]
    if rna:
        t_c = ntc_dic["U"]
    else:
        t_c = ntc_dic["T"]
    if a_c + t_c == 0:
        return 0
    else:
        return (a_c - t_c) / (a_c + t_c)


################################################################################

def calc_seq_at_content_ntc(seq_l, ntc_dic,
                            rna=False):
    """
    Calculate AT content of a given a nucleotides count dictionary.

    >>> ntc_dic = {'A': 2, 'C': 2, 'G': 2, 'T': 2}
    >>> calc_seq_at_content_ntc(8, ntc_dic)
    0.5
    >>> ntc_dic = {'A': 1, 'C': 2, 'G': 2, 'T': 0}
    >>> calc_seq_at_content_ntc(5, ntc_dic)
    0.2

    """
    assert ntc_dic, "given ntc_dic empty"
    assert seq_l, "given sequence length seq_l 0 or empty"

    a_c = ntc_dic["A"]
    if rna:
        t_c = ntc_dic["U"]
    else:
        t_c = ntc_dic["T"]

    return (a_c + t_c) / seq_l


################################################################################

def calc_seq_entropy(seq_l, ntc_dic,
                     k=1):
    """
    Given a dictionary of nucleotide counts for a sequence ntc_dic and
    the length of the sequence seq_l, compute the Shannon entropy of
    the sequence.

    Formula (see CE formula) taken from:
    https://www.ncbi.nlm.nih.gov/pubmed/15215465

    >>> seq_l = 8
    >>> ntc_dic = {'A': 8, 'C': 0, 'G': 0, 'U': 0}
    >>> calc_seq_entropy(seq_l, ntc_dic)
    0
    >>> ntc_dic = {'A': 4, 'C': 4, 'G': 0, 'U': 0}
    >>> calc_seq_entropy(seq_l, ntc_dic)
    0.5
    >>> ntc_dic = {'A': 2, 'C': 2, 'G': 2, 'U': 2}
    >>> calc_seq_entropy(seq_l, ntc_dic)
    1.0
    >>> seq_l = 4
    >>> ntc_dic = {'AA': 3, 'AC': 0, 'AG': 0, 'AT': 0, 'CA': 0, 'CC': 0, 'CG': 0, 'CT': 0, 'GA': 0, 'GC': 0, 'GG': 0, 'GT': 0, 'TA': 0, 'TC': 0, 'TG': 0, 'TT': 0}
    >>> calc_seq_entropy(seq_l, ntc_dic, k=2)
    0
    >>> seq_l = 17
    >>> ntc_dic = {'AA': 1, 'AC': 1, 'AG': 1, 'AT': 1, 'CA': 1, 'CC': 1, 'CG': 1, 'CT': 1, 'GA': 1, 'GC': 1, 'GG': 1, 'GT': 1, 'TA': 1, 'TC': 1, 'TG': 1, 'TT': 1}
    >>> calc_seq_entropy(seq_l, ntc_dic, k=2)
    1.0

    """
    assert seq_l, "invalid sequence length seq_l given"

    # Shannon entropy.
    ce = 0
    n_kmers = len(ntc_dic)  # number of k-mers, .e.g. 4 for k=1, 16 for k=2.
    for nt in ntc_dic:
        c = ntc_dic[nt]
        total_c = seq_l - k + 1 # total number of k-mers in sequence.

        if c != 0:
            ce += (c/total_c) * log((c/total_c), n_kmers)
    if ce == 0:
        return 0
    else:
        return -1*ce


################################################################################

def ntc_dic_to_ratio_dic(ntc_dic,
                         perc=False):
    """
    Given a dictionary of nucleotide counts, return dictionary of nucleotide
    ratios (count / total nucleotide number).

    perc:
    If True, make percentages out of ratios (*100).

    >>> ntc_dic = {'A': 5, 'C': 2, 'G': 2, 'T': 1}
    >>> ntc_dic_to_ratio_dic(ntc_dic)
    {'A': 0.5, 'C': 0.2, 'G': 0.2, 'T': 0.1}

    """
    assert ntc_dic, "given dictionary ntc_dic empty"
    # Get total number.
    total_n = 0
    for nt in ntc_dic:
        total_n += ntc_dic[nt]
    ntr_dic = {}
    for nt in ntc_dic:
        ntc = ntc_dic[nt]
        ntr = ntc / total_n
        if perc:
            ntr = ntr*100
        ntr_dic[nt] = ntr
    return ntr_dic


################################################################################

def seqs_dic_count_kmer_freqs(seqs_dic, k,
                              rna=False,
                              perc=False,
                              return_ratios=False,
                              report_key_error=True,
                              skip_non_dic_keys=False,
                              convert_to_uc=False):
    """
    Given a dictionary with sequences seqs_dic, count how many times each
    k-mer is found over all sequences (== get k-mer frequencies).
    Return k-mer frequencies count dictionary.
    By default, a DNA dictionary is used, and key errors will be reported.

    rna:
    Instead of DNA dictionary, use RNA dictionary (ACGU) for counting
    k-mers.
    perc:
    If True, make percentages out of ratios (*100).
    return_ratios:
    Return k-mer ratios instead of frequencies (== counts).
    report_key_error:
    If True, report key error (di-nucleotide not in count_dic).
    convert_to_uc:
    Convert sequences to uppercase before counting.
    skip_non_dic_keys:
    Skip k-mers not in count_dic. These usually are N-containing k-mers.
    By default, key errors are reported.


    >>> seqs_dic = {'seq1': 'AACGTC', 'seq2': 'GGACT'}
    >>> seqs_dic_count_kmer_freqs(seqs_dic, 2)
    {'AA': 1, 'AC': 2, 'AG': 0, 'AT': 0, 'CA': 0, 'CC': 0, 'CG': 1, 'CT': 1, 'GA': 1, 'GC': 0, 'GG': 1, 'GT': 1, 'TA': 0, 'TC': 1, 'TG': 0, 'TT': 0}
    >>> seqs_dic = {'seq1': 'AAACGT'}
    >>> seqs_dic_count_kmer_freqs(seqs_dic, 2, return_ratios=True, perc=True)
    {'AA': 40.0, 'AC': 20.0, 'AG': 0.0, 'AT': 0.0, 'CA': 0.0, 'CC': 0.0, 'CG': 20.0, 'CT': 0.0, 'GA': 0.0, 'GC': 0.0, 'GG': 0.0, 'GT': 20.0, 'TA': 0.0, 'TC': 0.0, 'TG': 0.0, 'TT': 0.0}

    """
    # Checks.
    assert seqs_dic, "given dictinary seqs_dic empty"
    assert k, "invalid k given"
    assert k > 0, "invalid k given"
    # Get k-mer dictionary.
    count_dic = get_kmer_dic(k, rna=rna)
    # Count k-mers for all sequences in seqs_dic.
    total_c = 0
    for seq_id in seqs_dic:
        seq = seqs_dic[seq_id]
        if convert_to_uc:
            seq = seq.upper()
        for i in range(len(seq)-k+1):
            kmer = seq[i:i+k]
            if skip_non_dic_keys:
                if kmer not in count_dic:
                    continue
            else:
                if report_key_error:
                    assert kmer in count_dic, "k-mer \"%s\" not in count_dic" %(kmer)
            if kmer in count_dic:
                count_dic[kmer] += 1
                total_c += 1
    assert total_c, "no k-mers counted for given seqs_dic (sequence lengths < set k ?)"

    # Calculate ratios.
    if return_ratios:
        for kmer in count_dic:
            ratio = count_dic[kmer] / total_c
            if perc:
                count_dic[kmer] = ratio*100
            else:
                count_dic[kmer] = ratio
    # Return k-mer counts or ratios.
    return count_dic


################################################################################

def get_ntc_dic(seq, rna=False):
    """
    Get single nucleotide count dictionary.

    >>> seq = 'ACGTACGT'
    >>> get_ntc_dic(seq)
    {'A': 2, 'C': 2, 'G': 2, 'T': 2}

    """

    ntc_dic = get_kmer_dic(1, rna=rna)

    for nt in seq:
        if nt in ntc_dic:
            ntc_dic[nt] += 1

    return ntc_dic


################################################################################

def get_kmer_counts_dic(seq, k, rna=False):
    """
    Get k-mer counts dictionary for sequence.

    >>> seq = 'ACGT'
    >>> get_kmer_counts_dic(seq, 2)
    {'AA': 0, 'AC': 1, 'AG': 0, 'AT': 0, 'CA': 0, 'CC': 0, 'CG': 1, 'CT': 0, 'GA': 0, 'GC': 0, 'GG': 0, 'GT': 1, 'TA': 0, 'TC': 0, 'TG': 0, 'TT': 0}
    >>> seq = 'ACGTA'
    >>> get_kmer_counts_dic(seq, 1)
    {'A': 2, 'C': 1, 'G': 1, 'T': 1}

    """

    count_dic = get_kmer_dic(k, rna=rna)
    total_c = 0

    for i in range(len(seq)-k+1):
        kmer = seq[i:i+k]
        if kmer in count_dic:
            count_dic[kmer] += 1
            total_c += 1
    
    assert total_c, "no k-mers counted for given sequence \"%s\" (sequence lengths < set k ?)" %(seq)

    return count_dic


################################################################################

def get_kmer_dic(k,
                 fill_idx=False,
                 rna=False):
    """
    Return a dictionary of k-mers. By default, DNA alphabet is used (ACGT).
    Value for each k-mer key is set to 0.

    rna:
    Use RNA alphabet (ACGU).

    >>> get_kmer_dic(1)
    {'A': 0, 'C': 0, 'G': 0, 'T': 0}
    >>> get_kmer_dic(2, rna=True)
    {'AA': 0, 'AC': 0, 'AG': 0, 'AU': 0, 'CA': 0, 'CC': 0, 'CG': 0, 'CU': 0, 'GA': 0, 'GC': 0, 'GG': 0, 'GU': 0, 'UA': 0, 'UC': 0, 'UG': 0, 'UU': 0}
    >>> get_kmer_dic(1, fill_idx=True)
    {'A': 1, 'C': 2, 'G': 3, 'T': 4}
    >>> get_kmer_dic(2, rna=True, fill_idx=True)
    {'AA': 1, 'AC': 2, 'AG': 3, 'AU': 4, 'CA': 5, 'CC': 6, 'CG': 7, 'CU': 8, 'GA': 9, 'GC': 10, 'GG': 11, 'GU': 12, 'UA': 13, 'UC': 14, 'UG': 15, 'UU': 16}

    """
    # Check.
    assert k, "invalid k given"
    assert k > 0, "invalid k given"
    # Dictionary.
    mer2c_dic = {}
    # Alphabet.
    nts = ["A", "C", "G", "T"]
    if rna:
        nts = ["A", "C", "G", "U"]
    # Recursive k-mer dictionary creation.
    def fill(i, seq, mer2c_dic):
        if i:
            for nt in nts:
                fill(i-1, seq+nt, mer2c_dic)
        else:
            mer2c_dic[seq] = 0
    fill(k, "", mer2c_dic)
    if fill_idx:
        idx = 0
        for kmer,c in sorted(mer2c_dic.items()):
            idx += 1
            mer2c_dic[kmer] = idx
    return mer2c_dic


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
        math domain error. 2.2e-308 pval -> 307.6575773191778 -log10pval.

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
        # assert len(df) == len(rbp_list), "len(df) != len(rbp_list) (%i != %i)" %(len(df), len(rbp_list)) 
        # for i,rbp_i in enumerate(rbp_list):
        #     for j,rbp_j in enumerate(rbp_list):
        #         if df.loc[rbp_i][rbp_j] is not None:
        #             pv = df.loc[rbp_i][rbp_j]
        #             if convert_zero_pv:
        #                 if pv == 0:
        #                     pv = min_pv
        #             ltf_pval = log_tf_pval(pv)
        #             df.loc[rbp_i, rbp_j] = ltf_pval

        for i in range(len(df)):
            for j in range(len(df)):
                if df.iloc[i, j] is not None:
                    pv = df.iloc[i, j]
                    if convert_zero_pv:
                        if pv == 0:
                            pv = min_pv
                    ltf_pval = log_tf_pval(pv)
                    df.iloc[i, j] = ltf_pval


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

def search_generate_html_motif_plots(args, search_rbps_dic, 
                                     seq_motif_blocks_dic, str_motif_blocks_dic,
                                     benchlib_path, motif2db_dic,
                                     rbp2motif2annot2c_dic=False,
                                     rbp2motif2annot2normc_dic=False,
                                     annot2color_dic=False,
                                     mrna_reg_occ_dic=False,
                                     norm_mrna_reg_dic=False,
                                     html_report_out="motif_plots.rbpbench_search.html",
                                     rbpbench_mode="search --plot-motifs",
                                     reg_seq_str="regions",
                                     goa_results_df=False,
                                     goa_stats_dic=False,
                                     goa_results_tsv="goa_results.tsv",
                                     id2pids_dic=False,
                                     id2exp_dic=False,
                                     match_c_dic=False,
                                     match_c_total_dic=False,
                                     plots_subfolder="html_motif_plots"):
    """
    Create motif plots for selected RBPs.

    """

    out_folder = args.out_folder
    # Use absolute paths?
    if args.plot_abs_paths:
        out_folder = os.path.abspath(out_folder)

    # Version string.
    version_str = "v" + args.version

    plots_folder = plots_subfolder
    plots_out_folder = out_folder + "/" + plots_folder
    if args.plot_abs_paths:
        plots_folder = plots_out_folder

    # Delete plots if already present.
    if os.path.exists(plots_out_folder):
        shutil.rmtree(plots_out_folder)
    os.makedirs(plots_out_folder)

    html_out = out_folder + "/" + "motif_plots.rbpbench_search.html"
    md_out = out_folder + "/" + "motif_plots.rbpbench_search.md"
    if html_report_out:
        html_out = html_report_out

    site_type_uc = "Genomic"
    site_type = "genomic"
    tr_modes = ["searchlongrna", "searchrna --plot-motifs"]
    seq_modes = ["searchseq --plot-motifs"]
    if rbpbench_mode in tr_modes:
        site_type_uc = "Transcript"
        site_type = "transcript"
    if rbpbench_mode in seq_modes:
        site_type_uc = "Sequence"
        site_type = "sequence"


    """
    Setup sorttable.js to make tables in HTML sortable.

    """
    sorttable_js_path = benchlib_path + "/content/sorttable.js"
    assert os.path.exists(sorttable_js_path), "sorttable.js not at %s" %(sorttable_js_path)
    sorttable_js_html = '<script src="' + sorttable_js_path + '" type="text/javascript"></script>'
    if args.sort_js_mode == 2:
        shutil.copy(sorttable_js_path, plots_out_folder)
        sorttable_js_path = plots_folder + "/sorttable.js"
        sorttable_js_html = '<script src="' + sorttable_js_path + '" type="text/javascript"></script>'
    elif args.sort_js_mode == 3:
        js_code = read_file_content_into_str_var(sorttable_js_path)
        sorttable_js_html = "<script>\n" + js_code + "\n</script>\n"

    # Logo path.
    logo_path_html = plots_folder + "/logo.png"
    if not os.path.exists(logo_path_html):
        logo_path = benchlib_path + "/content/logo.png"
        assert os.path.exists(logo_path), "logo.png not found in %s" %(logo_path)
        shutil.copy(logo_path, plots_out_folder)

    # HTML head section.
    html_head = """<!DOCTYPE html>
<html>
<head>
<title>RBPBench - Motif Plots and Hit Statistics</title>

<style>
    th, td {
        border: 1px solid black;
        padding: 8px;
        text-align: center;
    }
    th {
        background-color: #f2f2f2;
    }
    .page {
        page-break-after: always;
    }
    .title-container {
        display: flex;
        align-items: center;
    }
    .title-container img {
        margin-right: 10px; /* Adjust the spacing as needed */
    }
</style>

</head>

<div class="title-container">
    <img src="%s" alt="Logo" width="175">
    <h1>Motif Plots and Hit Statistics</h1>
</div>


<body>
""" %(logo_path_html)


    # HTML tail section.
    html_tail = """
%s
</body>
</html>
""" %(sorttable_js_html)

    # Markdown part.
    mdtext = """

List of available motif hit statistics and motif plots generated
by RBPBench (%s, rbpbench %s):

- [Motif hit statistics](#motif-hit-stats)
""" %(version_str, rbpbench_mode)

    if args.run_goa_tr:
        mdtext += "- [Motif hit GO enrichment analysis results](#goa-results)\n"

    motif_plot_ids_dic = {}
    idx = 0
    for rbp_id, rbp in sorted(search_rbps_dic.items()):
        idx += 1
        tab_id = "plot-%i" %(idx)
        motif_plot_ids_dic[rbp_id] = tab_id
        mdtext += "- [%s motifs](#%s)\n" %(rbp_id, tab_id)

    add_head_info = ""
    if args.bed_sc_thr is not None:
        add_head_info = " BED score threshold (--bed-sc-thr) = %s" %(str(args.bed_sc_thr))
        if args.bed_sc_thr_rev_filter:
            add_head_info += " (reverse filtering applied, i.e., the lower the better)."
        else:
            add_head_info += "."

    mdtext += "\nFIMO p-value threshold (--fimo-pval) = %s.%s # considered input %s = %i." %(str(args.fimo_pval), add_head_info, reg_seq_str, args.c_regions)
    if args.ext_up is not None:
        mdtext += " Region extension (upstream, downstream) = (%i, %i).\n" %(args.ext_up, args.ext_down)
    mdtext += "\n"

    mdtext += "\n&nbsp;\n"


    """
    Motif hit statistics table.

    """

    gh_info = ""
    if args.greatest_hits:
        gh_info = "--greatest-hits enabled, i.e., only highest scoring hits are reported for each input region."

    mdtext += """
## Motif hit statistics ### {#motif-hit-stats}

**Table:** RBP motif hit statistics with RBP ID, motif ID, motif database ID (set to "user" if user-supplied motif, otherwise internal database ID), 
and respective number of motif hits found in supplied %s regions.
%s

""" %(site_type, gh_info)

    # mdtext += '| &nbsp; RBP ID &nbsp; | &nbsp; Motif ID &nbsp; | Motif database | # motif hits |' + " \n"
    # mdtext += "| :-: | :-: | :-: | :-: |\n"

    mdtext += '<table style="max-width: 750px; width: 100%; border-collapse: collapse; line-height: 0.8;">' + "\n"
    mdtext += "<thead>\n"
    mdtext += "<tr>\n"
    mdtext += "<th>RBP ID</th>\n"
    mdtext += "<th>Motif ID</th>\n"
    mdtext += "<th>Motif database</th>\n"
    mdtext += "<th># motif hits</th>\n"
    mdtext += "</tr>\n"
    mdtext += "</thead>\n"
    mdtext += "<tbody>\n"

    for rbp_id, rbp in sorted(search_rbps_dic.items()):
        for idx, motif_id in enumerate(rbp.seq_motif_ids):
            c_motif_hits = rbp.seq_motif_hits[idx]
            motif_db = motif2db_dic[motif_id]
            # mdtext += "| %s | %s | %s | %i |\n" %(rbp_id, motif_id, motif_db, c_motif_hits)
            mdtext += "<tr>\n"
            mdtext += "<td>%s</td>\n" %(rbp_id)
            mdtext += "<td>%s</td>\n" %(motif_id)
            mdtext += "<td>%s</td>\n" %(motif_db)
            mdtext += "<td>%i</td>\n" %(c_motif_hits)
            mdtext += "</tr>\n"

        for idx, motif_id in enumerate(rbp.str_motif_ids):
            c_motif_hits = rbp.str_motif_hits[idx]
            motif_db = motif2db_dic[motif_id]
            # mdtext += "| %s | %s | %s | %i |\n" %(rbp_id, motif_id, motif_db, c_motif_hits)
            mdtext += "<tr>\n"
            mdtext += "<td>%s</td>\n" %(rbp_id)
            mdtext += "<td>%s</td>\n" %(motif_id)
            mdtext += "<td>%s</td>\n" %(motif_db)
            mdtext += "<td>%i</td>\n" %(c_motif_hits)
            mdtext += "</tr>\n"
            
    mdtext += "</tbody>\n"
    mdtext += "</table>\n"

    mdtext += "\n&nbsp;\n&nbsp;\n"
    mdtext += "\nColumn IDs have the following meanings: "
    mdtext += "**RBP ID** -> RBP ID from database or user-defined (typically RBP name), "
    mdtext += "**Motif ID** -> Motif ID from database or user-defined, "
    mdtext += "**Motif database** -> Motif database used for search run, "
    mdtext += '**# motif hits** -> number of unique individual motif hits (i.e., unique hits for motif with motif ID).' + "\n"
    mdtext += "\n&nbsp;\n"

    """
    GOA results on transcripts (underlying genes) with motif hits.

    """

    if args.run_goa_tr:

        mdtext += """
## Motif hit GO enrichment analysis results ### {#goa-results}

"""
        c_goa_results = 0
        if isinstance(goa_results_df, pd.DataFrame) and not goa_results_df.empty:
            c_goa_results = len(goa_results_df)

        filter_purified_info = "GO terms with significantly higher and lower concentration ([e,p]) in study group are shown."
        filter_purified_info2 = "significant"
        if args.goa_filter_purified:
            filter_purified_info = "Only GO terms with significantly higher concentration in study group are shown."
            filter_purified_info2 = "significantly enriched"
            c_goa_results = len(goa_results_df[goa_results_df["enrichment"] == "e"])

        filter_further_info = ""
        if args.goa_max_child is not None: 
            filter_further_info += " Only GO terms with <= %i children are shown." %(args.goa_max_child)
        if args.goa_min_depth is not None:
            filter_further_info += " Only GO terms with >= %i depth are shown." %(args.goa_min_depth)
        if filter_further_info:
            filter_further_info += " Note that additional filters (children + depth) can result in an empty table. For all significant GO terms (i.e., unfiltered results) check *%s* output table." %(goa_results_tsv)

        goa_rna_region_info = "genes"
        if args.goa_rna_region == 1:
            goa_rna_region_info = "transcripts"
        elif args.goa_rna_region == 2:
            goa_rna_region_info = "3'UTR regions"
        elif args.goa_rna_region == 3:
            goa_rna_region_info = "CDS regions"
        elif args.goa_rna_region == 4:
            goa_rna_region_info = "5'UTR regions"

        goa_only_cooc_info = ""
        if args.goa_cooc_mode == 2:
            goa_only_cooc_info = " Only target genes with motif hits from any specified RBP (including regex) are considered."
        elif args.goa_cooc_mode == 3:
            goa_only_cooc_info = " Only target genes with motif hits from all specified RBPs (including regex) are considered."


        if c_goa_results > 0:

            mdtext += """
**Table:** GO enrichment analysis results for %s with motif hits, taking the corresponding genes for analysis. # of %s GO terms found: %i. Filter p-value threshold (on corrected p-value) = %s. # of target genes used for GOA: %i. # of background genes used for GOA: %i. 
%s %s %s

""" %(goa_rna_region_info, filter_purified_info2, c_goa_results, str(goa_stats_dic["pval_thr"]), goa_stats_dic["c_target_genes_goa"], goa_stats_dic["c_background_genes_goa"], goa_only_cooc_info, filter_purified_info, filter_further_info)

            mdtext += '<table style="max-width: 1200px; width: 100%; border-collapse: collapse; line-height: 0.9;">' + "\n"
            mdtext += "<thead>\n"
            mdtext += "<tr>\n"
            mdtext += "<th>GO</th>\n"
            mdtext += "<th>Term</th>\n"
            mdtext += "<th>Class</th>\n"
            mdtext += "<th>p-value</th>\n"
            mdtext += "<th>[e,p]</th>\n"
            mdtext += "<th>Depth</th>\n"
            mdtext += "<th># child</th>\n"
            mdtext += "<th># genes</th>\n"
            mdtext += "<th># study</th>\n"
            mdtext += "<th>% genes</th>\n"
            mdtext += "</tr>\n"
            mdtext += "</thead>\n"
            mdtext += "<tbody>\n"

            for index, row in goa_results_df.iterrows():

                go_id = row['GO']
                go_term = row['term']
                go_class = row['class']
                go_p = row['p']
                go_p_corr = row['p_corr']
                go_enrichment = row['enrichment']
                go_depth = row['depth']
                go_n_genes = row['n_genes']
                go_n_study = row['n_study']
                go_perc_genes = row['perc_genes']
                go_n_children = row['n_children']

                if args.goa_filter_purified:
                    if go_enrichment == "p":
                        continue
                if args.goa_max_child is not None:
                    if go_n_children > args.goa_max_child:
                        continue
                if args.goa_min_depth is not None:
                    if go_depth < args.goa_min_depth:
                        continue

                mdtext += '<tr>' + "\n"
                mdtext += "<td>" + go_id + "</td>\n"
                mdtext += "<td>" + go_term + "</td>\n"
                mdtext += "<td>" + go_class + "</td>\n"
                mdtext += "<td>" + str(go_p_corr) + "</td>\n"
                mdtext += "<td>" + go_enrichment + "</td>\n"
                mdtext += "<td>" + str(go_depth) + "</td>\n"
                mdtext += "<td>" + str(go_n_children) + "</td>\n"
                mdtext += "<td>" + str(go_n_genes) + "</td>\n"
                mdtext += "<td>" + str(go_n_study) + "</td>\n"
                mdtext += "<td>" + str(go_perc_genes) + "</td>\n"
                mdtext += '</tr>' + "\n"

            mdtext += '</tbody>' + "\n"
            mdtext += '</table>' + "\n"
            
            mdtext += "\n&nbsp;\n&nbsp;\n"
            mdtext += "\nColumn IDs have the following meanings: "
            mdtext += "**GO** -> gene ontology (GO) ID, "
            mdtext += "**Term** -> GO term / name, "
            mdtext += "**Class** -> GO term class (biological_process, molecular_function, or cellular_component), "
            mdtext += "**p-value** -> multiple testing corrected (BH) p-value, "
            mdtext += "**[e,p]** -> e: enriched, i.e., GO term with significantly higher concentration, p: purified, GO term with significantly lower concentration), "
            mdtext += "**Depth** -> depth / level of GO term in GO hierarchy (the higher number, the more specific), "
            mdtext += "**# child** -> number of GO term children, "
            mdtext += "**# genes** -> number of genes associated with GO term, "
            mdtext += "**# study** -> number of genes in study (i.e., target genes), "
            mdtext += "**% genes** -> percentage of study genes associated with GO term." + "\n"
            mdtext += "\n&nbsp;\n"

        else:

            if "c_target_genes_goa" in goa_stats_dic:

                mdtext += """

No %s GO terms found given p-value threshold of %s. # of target genes used for GOA: %i. # of background genes used for GOA: %i.

&nbsp;

""" %(filter_purified_info2, str(goa_stats_dic["pval_thr"]), goa_stats_dic["c_target_genes_goa"], goa_stats_dic["c_background_genes_goa"])

            else:

                mdtext += """

No significant GO terms found due to no GO IDs associated with target genes. # of initial target genes (i.e., genes overlapping with --in regions): %i.

&nbsp;

""" %(goa_stats_dic["c_target_genes_pre_filter"])


    """
    Motif plots.

    """

    for rbp_id, rbp in sorted(search_rbps_dic.items()):
        tab_id = motif_plot_ids_dic[rbp_id]

        # Count number of total motif hits for RBP.
        c_total_rbp_hits = 0

        motif_ids_list = []

        for idx, motif_id in enumerate(rbp.seq_motif_ids):
            c_total_rbp_hits += rbp.seq_motif_hits[idx]
            motif_ids_list.append(motif_id)
        for idx, motif_id in enumerate(rbp.str_motif_ids):
            c_total_rbp_hits += rbp.str_motif_hits[idx]
            motif_ids_list.append(motif_id)

        seq_motif_info = "RBP \"%s\" sequence motif plots." %(rbp_id)
        if rbp2motif2annot2c_dic and c_total_rbp_hits:
            seq_motif_info = "RBP \"%s\" sequence motif plots and %s region annotations for motif hits." %(rbp_id, site_type)
        if mrna_reg_occ_dic and c_total_rbp_hits:
            seq_motif_info = "RBP \"%s\" sequence motif plots and mRNA region annotations for motif hits." %(rbp_id)

        # RBP has sequence motifs?
        if rbp.seq_motif_ids and rbp_id != args.regex_id:
            mdtext += """
## %s motifs ### {#%s}

%s

""" %(rbp_id, tab_id, seq_motif_info)

        elif rbp.seq_motif_ids and rbp_id == args.regex_id:

            motif_id = rbp.seq_motif_ids[0]
            c_motif_hits = rbp.seq_motif_hits[0]

            mdtext += """
## %s motifs ### {#%s}

Motif ID / Regex: "%s". Number of unique motif hits in supplied %s regions: %i.

""" %(rbp_id, tab_id, motif_id, site_type, c_motif_hits)

        else:  # structure motifs.
            mdtext += """
## %s motifs ### {#%s}

RBP "%s" only contains structure motifs, which are currently not available for plotting.

""" %(rbp_id, tab_id, rbp_id)

        # If mRNA annotations given via mrna_reg_occ_dic.
        if mrna_reg_occ_dic and c_total_rbp_hits:

            assert annot2color_dic, "given mrna_reg_occ_dic but annot2color_dic is empty"

            mrna_occ_stacked_plot =  "mRNA_region_occ_stacked_plot.%s.png" %(rbp_id)
            mrna_occ_stacked_plot_out = plots_out_folder + "/" + mrna_occ_stacked_plot

            create_mrna_region_occ_plot(motif_ids_list, mrna_reg_occ_dic, 
                                        annot2color_dic, mrna_occ_stacked_plot_out,
                                        plot_pdf=args.plot_pdf,
                                        rbp_id=rbp_id)

            plots_path = plots_folder + "/" + mrna_occ_stacked_plot

            mdtext += '<img src="' + plots_path + '" alt="mRNA region occupancy stacked plot"' + "\n"
            mdtext += 'title="mRNA region occupancy stacked plot" />' + "\n"
            mdtext += """
**Figure:** mRNA region motif hit coverage profiles for RBP "%s" motif hits.
Motif hit coverage profiles are shown for all motifs of RBP "%s" combined, as well as single motifs (unless there is only one motif), over 5'UTR, CDS, and 3'UTR regions of mRNA.
x-axis is the motif hit coverage, i.e., how many motif hits found over the mRNA regions.
Only motif hit center positions are used for the annotation and coverage profile.
mRNA region lengths used for plotting are the %s region lengths obtained from the GTF file (5'UTR = %s, CDS = %s, 3'UTR = %s).
Number of mRNA sequences used for prediction and plot generation: %i.

&nbsp;

""" %(rbp_id, rbp_id, norm_mrna_reg_dic["mode"], str(norm_mrna_reg_dic["5'UTR"]), str(norm_mrna_reg_dic["CDS"]), str(norm_mrna_reg_dic["3'UTR"]), norm_mrna_reg_dic["c_mrna_seqs"])

        # If there are motif hit region annotations and hits for the RBP.
        if rbp2motif2annot2c_dic and c_total_rbp_hits:

            assert annot2color_dic, "given rbp2motif2annot2c_dic but annot2color_dic is empty"

            annot_stacked_bars_plot =  "annotation_stacked_bars_plot.%s.png" %(rbp_id)
            annot_stacked_bars_plot_out = plots_out_folder + "/" + annot_stacked_bars_plot

            create_annotation_stacked_bars_plot(rbp_id, rbp2motif2annot2c_dic, annot2color_dic,
                                                annot_stacked_bars_plot_out,
                                                plot_pdf=args.plot_pdf,
                                                x_label="Annotation overlap")

            plot_path = plots_folder + "/" + annot_stacked_bars_plot

            mdtext += '<img src="' + plot_path + '" alt="Region annotations plot"' + "\n"
            # mdtext += 'title="Annotation stacked bars plot" width="800" />' + "\n"
            mdtext += 'title="Region annotations plot" width="1050" />' + "\n"
            mdtext += """
**Figure:** Genomic region annotations for RBP "%s" motif hits.
%s motif hit regions are overlapped with genomic regions from GTF file and genomic region feature with highest overlap
is assigned to each motif hit region. "intergenic" feature means none of the used GTF region features overlap with the motif hit region.
Genomic annotations are shown for all motifs of RBP "%s" combined, as well as for the single motifs (unless there is only one motif).

&nbsp;

""" %(rbp_id, site_type_uc, rbp_id)


        # If there are motif hit region annotations and hits for the RBP.
        if rbp2motif2annot2normc_dic and c_total_rbp_hits:

            assert annot2color_dic, "given rbp2motif2annot2normc_dic but annot2color_dic is empty"

            annot_stacked_bars_plot =  "annotation_stacked_bars_plo.normc.%s.png" %(rbp_id)
            annot_stacked_bars_plot_out = plots_out_folder + "/" + annot_stacked_bars_plot

            create_annotation_stacked_bars_plot(rbp_id, rbp2motif2annot2normc_dic, annot2color_dic,
                                                annot_stacked_bars_plot_out,
                                                no_x_labels=True,
                                                plot_pdf=args.plot_pdf,
                                                x_label="Normalized annotation overlap")

            plot_path = plots_folder + "/" + annot_stacked_bars_plot

            mdtext += '<img src="' + plot_path + '" alt="Region annotations normalized counts plot"' + "\n"
            # mdtext += 'title="Annotation stacked bars plot" width="800" />' + "\n"
            mdtext += 'title="Region annotations normalized counts plot"  width="1050" />' + "\n"
            mdtext += """
**Figure:** Genomic region annotations for RBP "%s" motif hits, normalized by annotation region lengths in input regions.
I.e., annotation counts from the above figure are normalized depending on how much the annotation covers the input regions.
This removes annotation region length biases introduced in long genomic input regions and can give a better idea of 
motif prevalences (given a reasonably large input size/number) in certain genomic regions (e.g. the motif tends to occur 
more often in intron, 3'UTR etc.). Unique input regions size (nt): %i (i.e., overlapping regions merged).

&nbsp;

""" %(rbp_id, args.eff_in_reg_size)

        # Regular expression section.
        if rbp_id == args.regex_id:
            if match_c_dic and match_c_total_dic:
                motif_id = rbp.seq_motif_ids[0]
                mdtext += get_match_seqs_html_table(rbp_id, motif_id, 
                                                    match_c_dic, match_c_total_dic,
                                                    add_motif_id_header=False,
                                                    top_n=args.top_n_matched)
            continue

        for idx, motif_id in enumerate(rbp.seq_motif_ids):

            if match_c_dic and match_c_total_dic:
                mdtext += get_match_seqs_html_table(rbp_id, motif_id, 
                                                    match_c_dic, match_c_total_dic,
                                                    add_motif_id_header=False,
                                                    top_n=args.top_n_matched)

            c_motif_hits = rbp.seq_motif_hits[idx]
            motif_db = motif2db_dic[motif_id]
            motif_plot = "%s.%s.png" %(rbp_id, motif_id)
            motif_plot_out = plots_out_folder + "/" + motif_plot
            plot_path = plots_folder + "/" + motif_plot

            # Check if motif in motif database folder.
            if args.motif_db_str:
                db_motif_path = benchlib_path + "/content/motif_plots/%s" %(motif_plot)
                if os.path.exists(db_motif_path):
                    shutil.copy(db_motif_path, motif_plot_out)
                    if args.plot_pdf:
                        create_motif_plot(motif_id, seq_motif_blocks_dic,
                                          motif_plot_out,
                                          plot_pdf=True,
                                          plot_png=False)

            if not os.path.exists(motif_plot_out):
                create_motif_plot(motif_id, seq_motif_blocks_dic,
                                  motif_plot_out,
                                  plot_pdf=args.plot_pdf,
                                  plot_png=True)

            motif_pids_info = "-"
            if id2pids_dic:
                pid_list = []
                if motif_id in id2pids_dic:
                    for pid in id2pids_dic[motif_id]:
                        if pid == "-":
                            continue
                        elif pid.startswith("http"):
                            pid_str = '<a href="%s">DOI</a>' %(pid)
                            pid_list.append(pid_str)
                        else:
                            try:
                                int(pid)
                                pid_str = '<a href="https://pubmed.ncbi.nlm.nih.gov/%s">%s</a>' %(pid, pid)
                                pid_list.append(pid_str)
                            except ValueError:
                                continue

                if pid_list:
                    motif_pids_info = ",".join(pid_list)

            motif_exp_info = "-"
            if id2exp_dic:
                if motif_id in id2exp_dic:
                    motif_exp_list = id2exp_dic[motif_id]
                    motif_exp_info = ",".join(motif_exp_list)

            mdtext += '<img src="' + plot_path + '" alt="' + "sequence motif plot %s" %(motif_id) + "\n"
            mdtext += 'title="' + "sequence motif plot %s" %(motif_id) + '" width="500" />' + "\n"
            mdtext += """

**Figure:** Sequence motif plot for motif ID "%s" (RBP ID: %s, motif database ID: %s). X-axis: motif position. Y-axis: nucleotide probability. Number of %s unique motif hits in supplied %s regions: %i.
Motif references (PubMed, DOI): %s. Motif source database / experiments: %s.

&nbsp;

""" %(motif_id, rbp_id, motif_db, motif_id, site_type, c_motif_hits, motif_pids_info, motif_exp_info)


        for idx, motif_id in enumerate(rbp.str_motif_ids):

            if match_c_dic and match_c_total_dic:
                mdtext += get_match_seqs_html_table(rbp_id, motif_id, 
                                                    match_c_dic, match_c_total_dic,
                                                    add_motif_id_header=False,
                                                    top_n=args.top_n_matched)

            c_motif_hits = rbp.str_motif_hits[idx]
            # NO STRUCTURE MOTIF PLOTTING IMPLEMENTED YET.
            # TO DO ...

    # print("Generate motif plots HTML ... ")

    # Convert mdtext to html.
    md2html = markdown(mdtext, extensions=['attr_list', 'tables'])

    # OUTMD = open(md_out,"w")
    # OUTMD.write("%s\n" %(mdtext))
    # OUTMD.close()

    html_content = html_head + md2html + html_tail
    OUTHTML = open(html_out,"w")
    OUTHTML.write("%s\n" %(html_content))
    OUTHTML.close()


################################################################################

def create_motif_plot(motif_id, 
                      seq_motif_blocks_dic,
                      motif_plot_out,
                      plot_pdf=False,
                      plot_png=True):
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
    if plot_png:
        plt.savefig(motif_plot_out, dpi=100)
    if plot_pdf and motif_plot_out.endswith('.png'):
        pdf_out = motif_plot_out[:-4] + '.pdf'
        plt.savefig(pdf_out, dpi=100)
    plt.close()


################################################################################

def get_match_c_total_dic(match_c_dic):
    """
    Get total counts for motif_id in match_c_dic.
    Formats: 
    match_c_total_dic[rbp_id][motif_id] = count
    match_c_dic[rbp_id][motif_id][matched_seq] = count

    >>> match_c_dic = {'RBP1': {'motif1': {'seq1': 2, 'seq2': 3}, 'motif2': {'seq3': 1}}}
    >>> get_match_c_total_dic(match_c_dic)
    {'RBP1': {'motif1': 5, 'motif2': 1}}
    
    """
    assert match_c_dic, "match_c_dic is empty"

    match_c_total_dic = {}
    for rbp_id in match_c_dic:
        match_c_total_dic[rbp_id] = {}
        for motif_id in match_c_dic[rbp_id]:
            if motif_id not in match_c_total_dic[rbp_id]:
                match_c_total_dic[rbp_id][motif_id] = 0
            for matched_seq in match_c_dic[rbp_id][motif_id]:
                match_c_total_dic[rbp_id][motif_id] += match_c_dic[rbp_id][motif_id][matched_seq]

    return match_c_total_dic
    

################################################################################

def get_match_seqs_html_table(rbp_id, motif_id, match_c_dic, match_c_total_dic,
                              add_motif_id_header=False,
                              top_n=10):
    """
    Create matched sequences count statistics HTML table.
    
    Formats: 
    match_c_total_dic[rbp_id][motif_id] = count
    match_c_dic[rbp_id][motif_id][matched_seq] = count
    """

    assert match_c_dic, "match_c_dic is empty"
    assert match_c_total_dic, "match_c_total_dic is empty"
    assert rbp_id in match_c_dic, "RBP ID %s not found in match_c_dic" %(rbp_id)
    assert motif_id in match_c_dic[rbp_id], "motif ID %s not found in match_c_dic" %(motif_id)
    assert rbp_id in match_c_total_dic, "RBP ID %s not found in match_c_total_dic" %(rbp_id)
    assert motif_id in match_c_total_dic[rbp_id], "motif ID %s not found in match_c_total_dic" %(motif_id)

    # If no hits, do not create table.
    if not match_c_dic[rbp_id][motif_id]:
        return ""

    header_text = ""
    if add_motif_id_header:
        # header_text = "### Motif ID \"%s\"" % (motif_id)
        header_text = "### %s" % (motif_id)
    mdtext = """
%s

**Table:** Matched sequence statistics for RBP ID "%s" + motif ID "%s" combination. 
Match count: number of times matched sequence appears in input regions.
Match percentage: match count divided by total number of motif ID "%s" motif hits.
Only top %i matched sequences are shown.

""" %(header_text, rbp_id, motif_id, motif_id, top_n)


    mdtext += '<table style="max-width: 750px; width: 100%; border-collapse: collapse; line-height: 0.8;">' + "\n"
    mdtext += "<thead>\n"
    mdtext += "<tr>\n"
    mdtext += "<th>RBP ID</th>\n"
    mdtext += "<th>Motif ID</th>\n"
    mdtext += "<th>Matched sequence</th>\n"
    mdtext += "<th>Match count</th>\n"
    mdtext += "<th>Match %</th>\n"
    mdtext += "</tr>\n"
    mdtext += "</thead>\n"
    mdtext += "<tbody>\n"

    total_c = match_c_total_dic[rbp_id][motif_id]

    count = 0
    for matched_seq, match_c in sorted(match_c_dic[rbp_id][motif_id].items(), key=lambda x: x[1], reverse=True):
        count += 1
        if count > top_n:
            break
        match_perc = 0.0
        if match_c > 0:
            match_perc = (float(match_c)/float(total_c)) * 100.0
            mdtext += "<tr>\n"
            mdtext += "<td>%s</td>\n" %(rbp_id)
            mdtext += "<td>%s</td>\n" %(motif_id)
            mdtext += "<td>%s</td>\n" %(matched_seq)
            mdtext += "<td>%i</td>\n" %(match_c)
            mdtext += "<td>%.2f</td>\n" %(match_perc)
            mdtext += "</tr>\n"

    mdtext += "</tbody>\n"
    mdtext += "</table>\n"

    # mdtext += "\n&nbsp;\n&nbsp;\n"
    # mdtext += "\nColumn IDs have the following meanings: "
    # mdtext += "**RBP ID** -> RBP ID from database or user-defined (typically RBP name), "
    # mdtext += "**Motif ID** -> Motif ID from database or user-defined, "
    # mdtext += "**Motif database** -> Motif database used for search run, "
    # mdtext += '**# motif hits** -> number of unique individual motif hits (i.e., unique hits for motif with motif ID).' + "\n"
    mdtext += "\n&nbsp;\n"

    return mdtext


################################################################################

def compare_generate_html_report(args,
                                 compare_methods_dic, compare_datasets_dic,
                                 rbp_stats_dic, motif_stats_dic,
                                 benchlib_path,
                                 html_report_out="report.rbpbench_compare.html",
                                 plots_subfolder="html_plots"):
    """
    Create comparison statistics and HTML report.

    """
    out_folder = args.out_folder
    plot_abs_paths = args.plot_abs_paths
    sort_js_mode = args.sort_js_mode

    # Version string.
    version_str = "v" + args.version

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

    """
    Setup sorttable.js to make tables in HTML sortable.

    """
    sorttable_js_path = benchlib_path + "/content/sorttable.js"
    assert os.path.exists(sorttable_js_path), "sorttable.js not at %s" %(sorttable_js_path)
    sorttable_js_html = '<script src="' + sorttable_js_path + '" type="text/javascript"></script>'
    if sort_js_mode == 2:
        shutil.copy(sorttable_js_path, plots_out_folder)
        sorttable_js_path = plots_folder + "/sorttable.js"
        sorttable_js_html = '<script src="' + sorttable_js_path + '" type="text/javascript"></script>'
    elif sort_js_mode == 3:
        js_code = read_file_content_into_str_var(sorttable_js_path)
        sorttable_js_html = "<script>\n" + js_code + "\n</script>\n"

    # Logo path.
    logo_path_html = plots_folder + "/logo.png"
    if not os.path.exists(logo_path_html):
        logo_path = benchlib_path + "/content/logo.png"
        assert os.path.exists(logo_path), "logo.png not found in %s" %(logo_path)
        shutil.copy(logo_path, plots_out_folder)

    # HTML head section.
    html_head = """<!DOCTYPE html>
<html>
<head>
<title>RBPBench - Motif Search Comparison Report</title>

<style>
    th, td {
        border: 1px solid black;
        padding: 8px;
        text-align: center;
    }
    th {
        background-color: #f2f2f2;
    }
    .page {
        page-break-after: always;
    }
    .title-container {
        display: flex;
        align-items: center;
    }
    .title-container img {
        margin-right: 10px; /* Adjust the spacing as needed */
    }
</style>

</head>

<div class="title-container">
    <img src="%s" alt="Logo" width="175">
    <h1>Search Comparison Report</h1>
</div>

<body>
""" %(logo_path_html)

    # HTML tail section.
    html_tail = """
%s
</body>
</html>
""" %(sorttable_js_html)

    # Markdown part.
    mdtext = """

List of available comparison statistics generated
by RBPBench (%s, rbpbench compare):

""" %(version_str)

    # Comparisons based on method ID.
    method_ids_dic = {}
    idx = 0
    for comp_id, data in sorted(compare_methods_dic.items()):
        # mdtext += "\n"
        if len(data) < 2:
            continue
        idx += 1
        tab_id = "method-%i" %(idx)
        method_ids_dic[comp_id] = tab_id
        mdtext += "- [%s method ID comparison](#%s)\n" %(comp_id, tab_id)

    # Comparisons based on data ID.
    data_ids_dic = {}
    idx = 0
    for comp_id, data in sorted(compare_datasets_dic.items()):
        # mdtext += "\n"
        if len(data) < 2:
            continue
        idx += 1
        tab_id = "data-%i" %(idx)
        data_ids_dic[comp_id] = tab_id
        mdtext += "- [%s data ID comparison](#%s)\n" %(comp_id, tab_id)

    mdtext += "\n&nbsp;\n"


    """
    Method comparisons.

    """

    for comp_id, method_dic in sorted(compare_methods_dic.items()):
        if len(method_dic) < 2:
            continue
        tab_id = method_ids_dic[comp_id]
        mdtext += """
## %s method ID comparison ### {#%s}

**Table:** RBP motif hit statistics for combined ID "%s" (includes data ID, motif database ID, RBP ID) over different methods (method ID column).

""" %(comp_id, tab_id, comp_id)

        # mdtext += '| Method ID | # regions | # motif hits | % regions with motifs | % motif nucleotides | # motif hits per 1000 nt |' + " \n"
        # mdtext += "| :-: | :-: | :-: | :-: | :-: | :-: |\n"

        mdtext += '<table style="max-width: 1200px; width: 100%; border-collapse: collapse; line-height: 0.8;">' + "\n"
        mdtext += "<thead>\n"
        mdtext += "<tr>\n"
        mdtext += "<th>Method ID</th>\n"
        mdtext += "<th># regions</th>\n"
        mdtext += "<th># motif hits</th>\n"
        mdtext += "<th>% regions with motifs</th>\n"
        mdtext += "<th>% motif nucleotides</th>\n"
        mdtext += "<th># motif hits per 1000 nt</th>\n"
        mdtext += "</tr>\n"
        mdtext += "</thead>\n"
        mdtext += "<tbody>\n"

        for method_id in method_dic:
            int_id = method_dic[method_id]
            c_regions = rbp_stats_dic[int_id].c_regions
            c_uniq_motif_hits = rbp_stats_dic[int_id].c_uniq_motif_hits
            perc_reg_with_hits = rbp_stats_dic[int_id].perc_reg_with_hits
            perc_uniq_motif_nts_eff_reg = rbp_stats_dic[int_id].perc_uniq_motif_nts_eff_reg
            uniq_motif_hits_cal_1000nt = rbp_stats_dic[int_id].uniq_motif_hits_cal_1000nt
            # mdtext += "| %s | %i | %i | %.2f | %.2f | %.2f |\n" %(method_id, c_regions, c_uniq_motif_hits, perc_reg_with_hits, perc_uniq_motif_nts_eff_reg, uniq_motif_hits_cal_1000nt)
            mdtext += "<tr>\n"
            mdtext += "<td>%s</td>\n" %(method_id)
            mdtext += "<td>%i</td>\n" %(c_regions)
            mdtext += "<td>%i</td>\n" %(c_uniq_motif_hits)
            mdtext += "<td>%.2f</td>\n" %(perc_reg_with_hits)
            mdtext += "<td>%.2f</td>\n" %(perc_uniq_motif_nts_eff_reg)
            mdtext += "<td>%.2f</td>\n" %(uniq_motif_hits_cal_1000nt)
            mdtext += "</tr>\n"
        
        mdtext += "</tbody>\n"
        mdtext += "</table>\n"

        mdtext += "\n&nbsp;\n&nbsp;\n"
        mdtext += "\nColumn IDs have the following meanings: "
        mdtext += "**Method ID** -> method ID set for dataset (typically peak calling method ID), "
        mdtext += '**# regions** -> number of peak regions used for motif search, '
        mdtext += '**# motif hits** -> number of unique motif hits in peak regions (removed double counts), '
        mdtext += '**% regions with motifs** -> percentage of peak regions with motif hits, '
        mdtext += '**% motif nucleotides** -> percentage of unique motif nucleotides over effective peak region size (overlapping regions merged), '
        mdtext += '**# motif hits per 1000 nt** -> number of motif hits over 1000 nt of called peak region size (overlapping regions NOT merged).' + "\n"
        mdtext += "\n&nbsp;\n"

        """
        Venn diagram for method ID comparison.

        """

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
                            plot_pdf=args.plot_pdf,
                            set1_label=method_ids[0],
                            set2_label=method_ids[1])
        elif len(method_ids) == 3:
            create_venn3_diagram(int_ids[0], int_ids[1], int_ids[2],
                            motif_stats_dic,
                            venn_plot_out,
                            plot_pdf=args.plot_pdf,
                            set1_label=method_ids[0],
                            set2_label=method_ids[1],
                            set3_label=method_ids[2])
        elif len(method_ids) > 3 and len(method_ids) <= 24:
            create_vennx_diagram(int_ids, method_ids,
                                 motif_stats_dic, venn_plot_out,
                                 plot_pdf=args.plot_pdf)
        else:
            assert False, "two many methods to compare (comp_id: %s). Please use less methods for plotting (current limit: 24)" %(comp_id)

        mdtext += '<img src="' + plot_path + '" alt="' + "dataset comparison plot %s" %(comp_id) + "\n"
        mdtext += 'title="' + "dataset comparison plot %s" %(comp_id) + '" width="600" />' + "\n"
        mdtext += """

**Figure:** Venn diagram of motif hit occurrences for the %i different methods (%s) with identical combined ID "%s" + corresponding percentages 
of total motif hits for each region (method exclusive and intersection(s)).
Any given motif hit can either be found only by one method, or be identified by any set (>=2) of methods (intersection areas).

&nbsp;

""" %(c_methods, method_ids_str, comp_id)


    """
    Data comparisons.

    """

    for comp_id, data_dic in sorted(compare_datasets_dic.items()):
        if len(data_dic) < 2:
            continue
        tab_id = data_ids_dic[comp_id]
        mdtext += """
## %s data ID comparison ### {#%s}

**Table:** RBP motif hit statistics for combined ID "%s" (includes method ID, motif database ID, RBP ID) over different datasets (data ID column).

""" %(comp_id, tab_id, comp_id)
        # mdtext += '| Data ID | # regions | # motif hits | % regions with motifs | % motif nucleotides | # motif hits per 1000 nt |' + " \n"
        # mdtext += "| :-: | :-: | :-: | :-: | :-: | :-: |\n"

        mdtext += '<table style="max-width: 1200px; width: 100%; border-collapse: collapse; line-height: 0.8;">' + "\n"
        mdtext += "<thead>\n"
        mdtext += "<tr>\n"
        mdtext += "<th>Data ID</th>\n"
        mdtext += "<th># regions</th>\n"
        mdtext += "<th># motif hits</th>\n"
        mdtext += "<th>% regions with motifs</th>\n"
        mdtext += "<th>% motif nucleotides</th>\n"
        mdtext += "<th># motif hits per 1000 nt</th>\n"
        mdtext += "</tr>\n"
        mdtext += "</thead>\n"
        mdtext += "<tbody>\n"

        for dataset_id in data_dic:
            int_id = data_dic[dataset_id]
            c_regions = rbp_stats_dic[int_id].c_regions
            c_uniq_motif_hits = rbp_stats_dic[int_id].c_uniq_motif_hits
            perc_reg_with_hits = rbp_stats_dic[int_id].perc_reg_with_hits
            perc_uniq_motif_nts_eff_reg = rbp_stats_dic[int_id].perc_uniq_motif_nts_eff_reg
            uniq_motif_hits_cal_1000nt = rbp_stats_dic[int_id].uniq_motif_hits_cal_1000nt
            # mdtext += "| %s | %i | %i | %.2f | %.2f | %.2f |\n" %(dataset_id, c_regions, c_uniq_motif_hits, perc_reg_with_hits, perc_uniq_motif_nts_eff_reg, uniq_motif_hits_cal_1000nt)
            mdtext += "<tr>\n"
            mdtext += "<td>%s</td>\n" %(dataset_id)
            mdtext += "<td>%i</td>\n" %(c_regions)
            mdtext += "<td>%i</td>\n" %(c_uniq_motif_hits)
            mdtext += "<td>%.2f</td>\n" %(perc_reg_with_hits)
            mdtext += "<td>%.2f</td>\n" %(perc_uniq_motif_nts_eff_reg)
            mdtext += "<td>%.2f</td>\n" %(uniq_motif_hits_cal_1000nt)
            mdtext += "</tr>\n"
        
        mdtext += "</tbody>\n"
        mdtext += "</table>\n"

        mdtext += "\n&nbsp;\n&nbsp;\n"
        mdtext += "\nColumn IDs have the following meanings: "
        mdtext += "**Data ID** -> data ID set for dataset (typically describing CLIP data, e.g. CLIP method + cell type combination), "
        mdtext += '**# regions** -> number of peak regions used for motif search, '
        mdtext += '**# motif hits** -> number of unique motif hits in peak regions (removed double counts), '
        mdtext += '**% regions with motifs** -> percentage of peak regions with motif hits, '
        mdtext += '**% motif nucleotides** -> percentage of unique motif nucleotides over effective peak region size (overlapping regions merged), '
        mdtext += '**# motif hits per 1000 nt** -> number of motif hits over 1000 nt of called peak region size (overlapping regions NOT merged).' + "\n"
        mdtext += "\n&nbsp;\n"

        """
        Venn diagram for data ID comparison.

        """

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
                            plot_pdf=args.plot_pdf,
                            set1_label=data_ids[0],
                            set2_label=data_ids[1])
        elif len(data_ids) == 3:
            create_venn3_diagram(int_ids[0], int_ids[1], int_ids[2],
                            motif_stats_dic,
                            venn_plot_out,
                            plot_pdf=args.plot_pdf,
                            set1_label=data_ids[0],
                            set2_label=data_ids[1],
                            set3_label=data_ids[2])
        elif len(data_ids) > 3 and len(data_ids) <= 24:
            create_vennx_diagram(int_ids, data_ids,
                                 motif_stats_dic, venn_plot_out,
                                 plot_pdf=args.plot_pdf)
        else:
            assert False, "two many datasets to compare (comp_id: %s). Please use less datasets for plotting (current limit: 24)" %(comp_id)

        method_id = comp_id.split(",")[0]

        mdtext += '<img src="' + plot_path + '" alt="' + "dataset comparison plot %s" %(comp_id) + "\n"
        mdtext += 'title="' + "dataset comparison plot %s" %(comp_id) + '" width="600" />' + "\n"
        mdtext += """

**Figure:** Venn diagram of motif hit occurrences for the %i different datasets (%s) with identical combined ID "%s" + corresponding percentages 
of total motif hits for each region (dataset exclusive and intersection(s)).
Any given motif hit can either be found only in one dataset, or be common to any >= 2 datasets (intersection areas).

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

    # OUTMD = open(md_out,"w")
    # OUTMD.write("%s\n" %(mdtext))
    # OUTMD.close()

    html_content = html_head + md2html + html_tail
    OUTHTML = open(html_out,"w")
    OUTHTML.write("%s\n" %(html_content))
    OUTHTML.close()


################################################################################

def create_venn2_diagram(int_id1, int_id2,
                         motif_stats_dic,
                         out_plot,
                         plot_pdf=False,
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

    plt.savefig(out_plot, dpi=150, bbox_inches='tight')

    if plot_pdf and out_plot.endswith('.png'):
        pdf_out = out_plot[:-4] + '.pdf'
        plt.savefig(pdf_out, dpi=150, bbox_inches='tight')

    plt.close()


################################################################################

def create_venn3_diagram(int_id1, int_id2, int_id3,
                         motif_stats_dic,
                         out_plot,
                         alpha=0.5,
                         plot_pdf=False,
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
    plt.savefig(out_plot, dpi=150, bbox_inches='tight')

    if plot_pdf and out_plot.endswith('.png'):
        pdf_out = out_plot[:-4] + '.pdf'
        plt.savefig(pdf_out, dpi=150, bbox_inches='tight')

    # plt.clf()
    plt.close()


################################################################################

def create_vennx_diagram(int_ids, set_labels,
                         motif_stats_dic, out_plot,
                         plot_pdf=False,
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
    plt.savefig(out_plot, dpi=150, bbox_inches='tight')

    if plot_pdf and out_plot.endswith('.png'):
        pdf_out = out_plot[:-4] + '.pdf'
        plt.savefig(pdf_out, dpi=150, bbox_inches='tight')

    plt.close()


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


"""


   "We could not understand because we were too far and could not remember
    because we were travelling in the night of first ages, of those ages
    that are gone, leaving hardly a sign - and no memories." J.C.



                 ^OMQQOO6|^OQQM6MOMMMOQQQQMMMMOMIQMQMQOOO6QO6O6QQMQQQMO6QMMOOQMMQQMQQOI66IOOQQMMM|QQO66
                 OOQQM!I|6OMOQMOQOQQQOMQQOMMMMQMQQMQOQOOQI^QOMOQMQMMQMMOOQOQQOQMQMQMMOOOOQQMOMQMQOQQQOO!
             .   .MOQQQM6OOQMMMMMMMOMQQQ66IQMMMMMQQOOMOIOO^6|QQMMQMQMOMQOOQOQQOOMQQMM6IOOQQMMQMMQMQQOQ6I
             .|. I6Q6MQQMQOMMMMMMMMMMMMOO|6MMQQMMQ!O6MO|IQ6O6MMMMQMMI6I6QMQOQMOO6MMQQQMI6QOQQQMMMMQQOQOO
                6|IQQOMMMMMMMMMMMMMQMQMO6O6MM6MQO|IOMOO6IMQQQMMMMQOOOQ6MQQMQMMMMQMMMMQMMQQQQMMQQQMMQOQOO
             ..Q^OIIOMMMQQMQMMMMMMMMMMMQ6OOQOOQQQQ|MQQOQMMQMMMQMQO6|6OOMMQMMMMMMMMMMMMQMMMQOMMMQMMQQQQMQ|
            |! QOMO6MMMQQ6QMQMMMMMQMMMQQMMQQQMMOQO||OQOQOMMMMMMOMQ6OQQMMMMMQMMMMQMMMMMMMQQMOQQQMMMMQMMQQ!^
            OQMQQMMMMMMMQOQMMMQMMMMMQMOQOM6OQQMQ!!I66OOO66QMMMOQ|QOOMMMMMMMMMMMMMMMMMMMMMMMMQMMMMMMQMQOQOO      ..
       .   !6 6MMMMMMMMMMMMQMMQQMQMQMMMMMQQQQOMO^I|IQOO^!|6QMMQOOO6OQMMMMMMMMQMMMMMMMMMMMMMMMMMMQMMQMOMQII^  Q^
       .| . .Q6MMMMMMIQMMMQQMMMOQMQMMMMQMQMOQQ6I6!6QO|..^|6OQOOI6OMMQMQMMMMMMMMMQQ6I6OMQMMMQMMMMMQMMQQQMQQO6
         M.OI666MMMMQOQMMQ6OMQOQQMMQQQMQMMQQQO!.OQQ|O..^!IOIIIOOQMMOQMQMMMMMMMQMMMMMMMMQQMQMMMMMMMMMQQ6OQMQ
           QM6MMMMQMQQMMMMOMOQQMQMMMMMM6OQMO6|OI66O|^^!||II!^^!!6O6O6QMMMQMMMMMQQMQMMMMMMMMQI6QMMMMOQMIM.QQQ
            |MMMOMQQ6MMMMQQMMMQQMMMMMMQOQ6|^!|6OI6O|!!|6I|!!^^^^!|I6OOOMQMQMMOQMMMMMMMMQMMQMQMQMMMMMQMQ|!6QII
           O .M66MOQMMMMMMMMOMQMMMMQMOQO6I6|^I6I||666III|!^^.^^.^!|I6OOQMQMMMQQQMMMMMMOOOMMQMQMMMMMMMQMQQ  QM
           I.6Q6OOOQQMMMMMMMQMMMMMMMMOQ6II6O66QQOOQOQOOI!^^^^.^^!I6OMQQMQQQMQQMMQMMMMMMQOMMMQMMMMMMMMQMMMQMOI
           66MQQQOQQ6MMMMMOMMMMMMQMMQ6OO6OOMMQMMMMQMMOO6|^^.^^^!6QQMMMMMMMQMMMQQQMMMMMMQMQQQMMMMMMMMMMMMMMQQOI
      |O^   QMQOQQOQ66QMQOOMQMQMQMMM6QO6666OO6OQMMO6OMQQI^.^.!|6QMMMMMMMMMOQMMQMMMMQQQQMMQQMMMMMMMMMMMMMMQMMOQ
           |MQMQOMOOOO6QMQMMMMMMMMMQQM|^^.!6^^..I|^I6OOII....!6MMMMMQI!66II6OQMQOOQOQQM6OMQMMMMMMMMMMMMMMMQMQO^
         66QQMMQQOQQMMMOQMMMMMMQMMQQMO!.. ..!!^!|||I!^^!^. .^!6OO66II|I6|66OO6666O6OQMQMMQMMMMMMMMMMQMMMMMMMQQ.
      ^    OMMMQOQ6QOMOQMMMMMQMMMMQQM|^. . ...........^^..  .!66!|!^^^.^^!|!!|!||6OOQMMMMMMMMMMMMMMMMMMMMMMQMQ.
           QQMMMMOQMQQOOMQ6MQOOQOQMMM|^..  ..  .. .. .^.... ^!I|!^!^..^^.^^^^^!|I66OQMQQMMMMMMMMQMMMMMMMMMMQMOO
           OOOMMMQMMMMQMQQQMMMMMMMQMQ!........ .  . ...^..  .!I|I!^^.^^^^^^!^^!!|6OQQMMMMOOMMQMQMMMMMMMMMMMQM QO .
           6|I6MMQMMQMOQQQMMMQMMMMMM6!^^^....... . ..^.^.^ ..^II!|!^^^^^^^!^!^!|6OOOQMQMMMMMMMMMMMMMMMMMQQMIMQ M^
           I .66QQQMMMMQOOO6QQMMMMMMO|!^^...^.... ...^^^.^..^!I6|I|^^^^^^^!!^!|I6OOOQOQQMMMMMMMMQMMMMMOMMMQM OQ.6
           .6  .6OOMMMMMQMMMMQMMMMMQQ|!^^^.^..^......^.^.^. ^!II|I6^^^^^!^!!!|IIOOOQOMMMMMMMMMMMMMMMMMMMMMOM| !M Q
            O|Q. .MMMMQMMMMMMMMMMMMQQ|!^.^..^^..^...^.^^..^.^!||||I^.^!^^!||!II6OQQQQMOMMMMMMMMMMMMMMMMMMMMQO  |I|6
          . .Q| ^^OMMMMMMMMMMMMMMMMMM|!^!^^^..^^...^.6QMI!^^!6OOQOO|^.^^^!!||6O6OQOMMMMMMMQMMMMMMMMMMMQMMMMQO   Q.Q
         !  6.QOQMQQMMMMMMMMMMMMMMQMQ!!!^^^^^..^.....!6||OIOOQMMMMM6^^^^!^!!|6OOOQQMMMMMMMMMMMMMMQQMMMQMOMQOI   Q|!.
              IOM. OMMMMMMMMMMMMQMQMM!!!!!^.^^... ...^^!^^IQMMMMQO6|^^^!!||I66OQOQQQMMMMMMMMMMMMMMMMMOQ6I|OMI QQQ
             I^|6O.QMMMMMMMMMMMMMMMMM!^!||^^^... ...^^^^..!III|6I|^^^^^^!|I|66QQQOQQMMMMMMMMMMMMMMM|QOOI6Q6| |M Q
         !.. I ^QMMMMMMMMMMMMMMMMMMMQ!^^||!!.^.^.^^...^..^^^..^!!^!^.^!!||6OOQOQOOQMMMMMMMMMMQMMQQOQMQQII|.|6Q.O^
        .^^|IO6QOQMMMMQMMMMMMMMMMMMMQM|!!!!|!^^^.^^^..^.!^....|||^!!^!!|IIO6OQMQQOQMMMMMMMMMQQMQMQMMMQQMOI6QOM^
         ^  66OQQQQMMMMMQOQMMMMMMMMMMMM|!!!||^!!!6QMQQQ66QOMMMMMMMQQQQQO6OOQQQQOOQQMMMMMMMMQQOQQMQQOOQMQQOO!O
          !  O6Q6QMMQMMOMOMMMMMMMMMMMMMM|!!|I|!I!!^^.^^^. .....^|!|||III6OQOQQQOQQMMMMMMMMMMMQQMQMMQMMMQQMM|
             OQQOQQQQO66QMMMMMMMQMMMMMMMM6|!|II!!^!^!^^!.^.^^^^|||6OO66OQOOQMQOQQMMMMMMQMMMQMQMMMQMQMMQMQ.M
              QOM..^IQQMMQMMM!QMMMMMMMMMMMQI!||II|||6OOQQOOQOQQQQMQQQQOOOQQQQMQMMMMMMMMMMMQQQMOMQMO^OO QM.
                6MQQOQQMMMMMIOMMMMMMMMMMMMMQO||!!||||!|I6OQQMQQMMQO6O6OOQQQMQQMMMMMMMMMMMMMMMQOMMMMM6QM
                    .   QQQM6QMMMMMMMOQQMMMOOQ6||!^!^^....^^!^I|I!||I|OOOQQQMMMMMMMMMMMMMQMMMMMMQIIQOM6                     .^^!^!^^^!
                        |OQMOQMMMMMQMM6QMMMIOOQQ|!!!!^^^^^^!!^^!!!|6I6OOQOMMMMMMMMMMMMMMMMMMMMMMMMMQM6!          .|I|!!!^^^^^^^^^.^.^^
                         MMMMIMQMMMQMQ6QMMM6666OOI||!!!!!^||!|!|I|6IIOQQMMMMMMMMMMMMMMMMMMMMMMMQMOQMM6.       QO6|!^^^.^^.^..^...^^...
                         QQMQQMMQMQIOOOQMMMOI6II6OO6I|IIIIII6II6II6OOQQMMMMMMMMMMMMMMMMMMMMMMMMMMMMM.      IQO6!^^^.^..^..^...^...^^.^
                          MQOQMQMMMQQOMMMMMQ|I||I6OOQQOOQQOOOOQQQQQQMMMMMMMMMMMMMMQMMMMMMMMQMMMMMMMMI    6O66|!^.^^..^...........^...^
.^^^^                   |I!MMMMQMMMQ6MOMMMMMIII|||I6OOQMMMQMMMQMMMMMMMMMMMMMMMMQQQQQMMMMMMMMMMMMMMMMQMO6QO6I|..^.........^......^.^..^
^^^.^^^^^^^^               !QMMMMQOMMQQQMQMMO||!!!|IOOOQMQMQMMMMMMQMMMMMMMMMMMQQOOOOQQMQQMMMMMMMMMMMMQQQO6I!^.^...............^.....^^
. .^...^^..^^^^^             .^|Q6MOMMMMQQMQI!!|!^!|II66OQQQQMQMQMQMMMMMMMMMQQOOOOOOOOOQ6OQMQQQMMMMQMMQ6||!^.......... ........!..^^.^
  .........^.^.^!^^|I!           ^OMMQMMMMOQ!!!!!^!^|!|66OQQQQQMQMMMMMMMMMMQOQO6O6O66OOOQ6O66QQOOQQ6I|!!^....... ....  .......^..^.^^!
.. . ....  ....^.^^^!|III6|          6QQQOQ6!!!!^^!!^!|I66OQQQMQQMMMMMMMMQQOOOOO6O66666^QI66I66I|!^^^...^^.^.^^^...... . ......^^.^!^!
  ..  . ... ......^^^^!!||II6|.   !QMMMMOQOI|^^^!^!^!!!!I666QOQQOQMQMQMQQOOO66O6666I6III6|IIIII||!^!^^^!!!^.^.......... .....^..^.^!^!
 . . .  ..  .....^.^^^^!^!||||6O6OMQMMQMOQ6I!!!^!^^^^^^!^|I66O6OOOOQOOOO6666I6I6IIII|I 6I|I||I|I||!!!!!!^^^^..^................^.^^^^!
. . .  ..   . ......^^^!!!!^!!!|IOOQQQ6Q6|||!!^!!^^^^^^!^!||II666O6666666I666II|I!|||IO|!|I|!||!!!^!^^^^^^....... this is the end
                                                                                                                  beautiful friend ...


"""
