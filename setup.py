#!/usr/bin/env python3

from setuptools import setup


""" Setup RBPBench """

setup(
    name='rbpbench',
    version='1.0.6',
    description='Evaluate CLIP-seq and other genomic region data using a comprehensive collection of known RBP binding motifs',
    long_description=open('README.md').read(),
    url='https://github.com/michauhl/RBPBench',
    license='MIT',
    scripts=['bin/rbpbench', 
             'bin/batch_table_wrapper_rbpbench.py', 
             'bin/gtf_extract_gene_region_bed.py', 
             'bin/gtf_get_mpt_nt_freqs.py', 
             'bin/gtf_get_mpt_with_introns_nt_freqs.py', 
             'bin/gtf_get_gene_region_nt_freqs.py', 
             'bin/gtf_extract_exon_intron_region_bed.py', 
             'bin/bed_shift_regions.py', 
             'bin/gtf_extract_mpt_region_bed.py', 
             'bin/bed_print_first_n_pos.py', 
             'bin/bed_print_last_n_pos.py',
             'bin/gtf_extract_exon_intron_border_bed.py',
             'bin/create_custom_meme_motif_db.py',
             'bin/gtf_extract_tr_feat_bed.py',
             'bin/bed_merge_ol_reg.py',
             'bin/get_genomic_conservation_scores.py',
             'bin/batch_get_common_dataset_gene_ids.py'
             ],
    packages=['rbpbench'],
    package_data={'rbpbench': ['content/*', 'content/motif_plots/*']},
    include_package_data=True,
    zip_safe=False,
)
