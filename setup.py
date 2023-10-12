#!/usr/bin/env python3

import re
import setuptools

long_description = open('README.md').read()

version = re.search(r'__version__\s*=\s*[\'"]([^\'"]*)[\'"]',
                    open('phylogenetics/__init__.py').read()).group(1)

setuptools.setup(
    name='VisProPhyl',
    version=version,
    author='Mathias Bockwoldt',
    author_email='mathias.bockwoldt@gmail.com',
    description='Map Blast results on a common-knowledge phylogenetic tree',
    long_description=long_description,
    long_description_content_type='text/markdown',
    url='https://github.com/MolecularBioinformatics/VisProPhyl',
    packages=setuptools.find_packages(),
	package_data = {'phylogenetics': ['templates/*']},
    entry_points={'console_scripts': [
                                    'phylogenetics = phylogenetics.cli:main',
                                    'phylotree = phylogenetics.phylotree:main',
                                    'blast2fasta = phylogenetics.blast2fasta:main'
                                    ]},
    install_requires=[
        'PyQt5>=5.11.3',
        'wheel>=0.33.0',
        'numpy>=1.15.1',
        'scipy>=1.11.0',
        'matplotlib>=3.1.1',
        'pandas>=1.0.0',
        'Pillow>=6.0.0',
        'biopython>=1.7.4',
        'ete3>=3.1.1',
	'taxfinder>=0.0.1'
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
