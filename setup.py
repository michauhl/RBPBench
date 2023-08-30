#!/usr/bin/env python3

from setuptools import setup


""" Setup RBPBench """

setup(
    name='rbpbench',
    version='0.2',
    description='Evaluate CLIP-seq and other genomic region data using a comprehensive collection of known RBP binding motifs',
    long_description=open('README.md').read(),
    author='Michael Uhl',
    author_email='uhlm@informatik.uni-freiburg.de',
    url='https://github.com/michauhl/RBPBench',
    license='MIT',
    scripts=['bin/rbpbench'],
    packages=['rbpbench'],
    package_data={'rbpbench': ['content/*']},
    zip_safe=False,
)
