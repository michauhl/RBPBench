#!/usr/bin/env python3

from setuptools import setup


""" Setup RBPBench """

setup(
    name='rbpbench',
    version='0.91',
    description='Evaluate CLIP-seq and other genomic region data using a comprehensive collection of known RBP binding motifs',
    long_description=open('README.md').read(),
    author='Michael Uhl',
    author_email='uhlm@informatik.uni-freiburg.de',
    url='https://github.com/michauhl/RBPBench',
    license='MIT',
    scripts=['bin/rbpbench', 'bin/batch_table_wrapper_rbpbench.py', 'bin/gtf_extract_gene_region_bed.py'],
    packages=['rbpbench'],
    package_data={'rbpbench': ['content/*', 'content/catrapid.omics.v2.1.human.6plus_motif_plots/*']},
    zip_safe=False,
)
