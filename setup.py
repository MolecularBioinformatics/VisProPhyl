from setuptools import setup

setup(
	name = 'Phylogenetics',
	version = '0.1.0',
	author = 'Mathias Bockwoldt',
	packages = ['venn', 'msa_blast', 'lineage_values'],
	scripts = ['bin/phylogenetics', 'bin/phylotree'],
	description = 'Running phylogenetic analyses on Blast results.',
	package_data = {'phylogenetics': ['templates/*']}
	)
