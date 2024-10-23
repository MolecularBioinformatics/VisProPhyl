#!/usr/bin/env python3

import re
import setuptools

long_description = open('README.md').read()

version = re.search(r'__version__\s*=\s*[\'"]([^\'"]*)[\'"]',
                    open('protaxovis/__init__.py').read()).group(1)

setuptools.setup(
    name='ProTaxoVis',
    version=version,
    author='Mathias Bockwoldt',
    author_email='mathias.bockwoldt@gmail.com',
    description='Map Blast results on a common-knowledge taxonomix (phylogenetic) tree',
    long_description=long_description,
    long_description_content_type='text/markdown',
    url='https://github.com/MolecularBioinformatics/ProTaxoVis',
    packages=setuptools.find_packages(),
	package_data = {'protaxovis': ['templates/*']},
    entry_points={'console_scripts': [
                                    'taxovis = protaxovis.cli:main',
                                    'taxotree = protaxovis.taxotree:main',
                                    'blast2fasta = protaxovis.blast2fasta:main'
                                    ]},
    install_requires=[
        'PyQt5 >= 5.11.3, < 6',
        'wheel >= 0.33.0, < 1',
        'numpy >= 1.15.1, < 2',
        'scipy >= 1.11.0, < 2',
        'matplotlib >= 3.1.1, < 4',
        'pandas == 1.*',
        'Pillow == 10.*',
        'biopython >= 1.74, <= 1.84',
        'ete3 == 3.1.*',
        'taxfinder >= 0.0.1'
    ],
    classifiers=[
        'Programming Language :: Python :: 3 :: Only',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'Natural Language :: English',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
    ],
    python_requires='>=3.6',
)
