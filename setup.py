from setuptools import setup

setup(
	name = 'Phylogenetics',
	version = '0.1.0',
	author = 'Mathias Bockwoldt',
	author_email = 'mathias.bockwoldt@uit.no',
	packages = ['phylogenetics', 'taxfinder'],
	scripts = ['bin/phylogenetics', 'bin/phylotree'],
	description = 'Running phylogenetic analyses on Blast results.',
	package_data = {'phylogenetics': ['templates/*']}
	)
