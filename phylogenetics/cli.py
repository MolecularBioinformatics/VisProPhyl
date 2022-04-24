#!/usr/bin/env python3

'''
This is the command line interface for phylogenetics. It is not supposed
to be imported.
'''

import argparse
import os
import sys
import textwrap

import pkg_resources

import matplotlib.pyplot as plt
from taxfinder import TaxFinder

import phylogenetics as phylo


TF = None
CR = None
blastdb = ''


class ConfigReader():
	'''
	The ConfigReader reads the `limits.txt`, `proteinlist.txt`, and
	`heatmap_config.txt` and returns the contents in different formats
	upon request.
	'''

	def __init__(self, limits_file='limits.txt', \
				protein_list_file='proteinlist.txt', \
				heatmap_config_file='heatmap_config.txt'):
		'''
		Initiates the class.

		:param limits_file: A string to be prepended before the name
		:param protein_list_file: A string to be appended after the name
		:param heatmap_config_file: A string to be appended after the name
		'''

		self.limits_dict = {}
		self.proteins = []
		self.protein_files = []
		self.taxa_to_show = []
		self.hm_colors = {}
		self.clustering_method = ''

		self._read_proteinlist(protein_list_file)
		self._read_limits(limits_file)
		self._read_heatmap_values(heatmap_config_file)


	def _read_limits(self, filename):
		'''
		Reads `filename` and saves the content as dictionary to be used
		by the `get_limits()` method.

		:param filename: The filename to read the limits from.
		'''

		with open(filename) as limit_file:
			self.limits_dict = {'default': (1e-30, 50)}
			for line in limit_file:
				line = line.split('#', maxsplit=1)[0].rstrip()
				if not line:
					continue
				lline = line.split()
				try:
					evalue = float(lline[1])
				except ValueError:
					evalue = self.limits_dict['default'][0]
				try:
					length = int(lline[2])
				except ValueError:
					length = self.limits_dict['default'][1]
				self.limits_dict[lline[0]] = (evalue, length)


	def get_limits(self):
		'''
		Reads a limits file and returns a dict in the format:
		{'protein': (evalue, length), ...} where evalue is the negative
		exponent of the evalue (int) and length is the length limit for
		the protein (int).

		:returns: A dictionary with limits as described above
		'''

		return self.limits_dict


	def _read_proteinlist(self, filename):
		'''
		Reads a protein list file and saves the content as lists to be
		used by the get* methods.

		:param filename: The filename to read the protein list from.
		'''

		with open(filename) as protein_list_file:
			self.proteins = []
			self.protein_files = []
			for line in protein_list_file:
				line = line.split('#', maxsplit=1)[0].strip()
				if not line:
					continue
				elems = line.split()
				if elems[0] not in self.proteins:
					self.proteins.append(elems[0])
					self.protein_files.append([])
				pidx = self.proteins.index(elems[0])
				self.protein_files[pidx].append(os.path.splitext(elems[1])[0])


	def get_protein_names(self, prefix = '', suffix = ''):
		'''
		Returns a set of proteins in the format:
		{'protein1', 'protein2', ...}. A prefix (a path) can be added as
		well as a suffix. get_protein_names('abc/', '.txt') will turn
		'protein1' into 'abc/protein1.txt'.

		:param prefix: A string to be prepended before the name
		:param suffix: A string to be appended after the name
		:returns: A set with proteins
		'''

		return set((prefix + protein + suffix for protein in self.proteins))


	def get_protein_files(self, prefix = '', suffix = ''):
		'''
		Returns a set of filenames in the format: {'proteinfile1',
		'proteinfile2', ...}. The files are just the basenames without
		suffix. A prefix (a path) can be added as well as a suffix.
		get_protein_files('abc/', '.txt') will turn 'file1' into
		'abc/file1.txt'.

		:param prefix: A string to be prepended before the filename
		:param suffix: A string to be appended after the filename
		:returns: A set of filenames
		'''

		return set((prefix + fn + suffix for proteinlist in self.protein_files for fn in proteinlist))


	def get_protein_dict(self, prefix = '', suffix = ''):
		'''
		Returns a dict with proteins and their respective file names in
		the format: {'protein': set('file1', 'file2', ...), ...}. The
		files are just the basenames without suffix. A prefix (a path)
		can be added as well as a suffix. get_protein_dict('abc/', '.txt')
		will turn 'file1' into 'abc/file1.txt'.

		:param prefix: A string to be prepended before the filename
		:param suffix: A string to be appended after the filename
		:returns: A dictionary as described above
		'''

		ret = {}
		for i, prot in enumerate(self.proteins):
			ret[prot] = set((prefix + fn + suffix for fn in self.protein_files[i]))

		return ret


	def _read_heatmap_values(self, filename):
		'''
		Reads heatmap config file and saves the content as lists and
		dicts to be used by the get* methods.
		'''

		with open(filename) as heatmap_file:
			self.taxa_to_show = []
			self.hm_colors = {}
			self.clustering_method = ''

			mode = ''
			modes = {'TAXA_TO_SHOW': 'taxa_to_show', 'COLORS': 'colors', 'ALGO': 'clustering'}

			clustering_methods = {'single', 'complete', 'average', 'weighted', 'centroid', 'median', 'ward'}

			for line in heatmap_file:
				line = line.split('#', maxsplit=1)[0].strip()
				if not line:
					continue

				if line in modes:
					mode = modes[line]

				elif mode == 'taxa_to_show':
					self.taxa_to_show.append(line)

				elif mode == 'colors':
					lline = line.split(maxsplit=1)
					self.hm_colors[lline[0]] = f'#{lline[1]}'

				elif mode == 'clustering':
					if line in clustering_methods:
						self.clustering_method = line
					else:
						raise ValueError(f'Error in heatmap_config.txt. Unknown clustering method: {line}!')

				else:
					raise ValueError(f'Error in heatmap_config.txt. No mode selected and {line} found!')


	def get_heatmap_taxa(self):
		'''
		Returns a list with the taxa to show.

		:returns: A list of taxa
		'''

		return self.taxa_to_show


	def get_heatmap_colors(self):
		'''
		Returns a dict with the colors in the form {'letter': 'color_code', ...}

		:returns: A dict with colors
		'''

		return self.hm_colors


	def get_heatmap_clustering(self):
		'''
		Returns the clustering method of choice

		:returns: A string with the clustering method
		'''

		return self.clustering_method

	# end ConfigReader


def _get_basename(name):
	'''
	Strips the path and the extension from a filename.
	E.g. /a/path/to/file.png -> file

	:param name: The file and path to be processed
	:returns: The processed filename
	'''

	return os.path.splitext(os.path.basename(name))[0]


def init():
	'''
	Creates some files and a folder to start with a new project.
	Existing files will not be overwritten.

	:creates: .gitignore
	:creates: limits.txt
	:creates: proteinlist.txt
	:creates: tree_config.txt
	:creates: tree_to_prune.txt
	:creates: heatmap_config.txt
	:creates: heatmap_template.html
	:creates: fastas/
	:creates: blastresults/
	'''

	os.makedirs('fastas', exist_ok=True)
	os.makedirs('blastresults', exist_ok=True)

	if not os.path.isfile('.gitignore'):
		with open('.gitignore', 'wb') as out:
			out.write(pkg_resources.resource_string('phylogenetics', 'templates/gitignore'))

	if not os.path.isfile('limits.txt'):
		with open('limits.txt', 'wb') as out:
			out.write(pkg_resources.resource_string('phylogenetics', 'templates/limits.txt'))

	if not os.path.isfile('proteinlist.txt'):
		with open('proteinlist.txt', 'wb') as out:
			out.write(pkg_resources.resource_string('phylogenetics', 'templates/proteinlist.txt'))

	if not os.path.isfile('tree_config.txt'):
		with open('tree_config.txt', 'wb') as out:
			out.write(pkg_resources.resource_string('phylogenetics', 'templates/tree_config.txt'))

	if not os.path.isfile('tree_to_prune.txt'):
		with open('tree_to_prune.txt', 'wb') as out:
			out.write(pkg_resources.resource_string('phylogenetics', 'templates/tree_to_prune.txt'))

	if not os.path.isfile('heatmap_config.txt'):
		with open('heatmap_config.txt', 'wb') as out:
			out.write(pkg_resources.resource_string('phylogenetics', 'templates/heatmap_config.txt'))

	if not os.path.isfile('heatmap_template.html'):
		with open('heatmap_template.html', 'wb') as out:
			out.write(pkg_resources.resource_string('phylogenetics', 'templates/heatmap_template.html'))


def run_blast():
	'''
	Blasts a list of files with Blastp against the provided database.

	:creates: `blastresults/*.xml`
	'''

	os.makedirs('blastresults', exist_ok=True)

	if not blastdb:
		raise ValueError('blastdb is empty. Run this script with -d blastdb')

	file_list = CR.get_protein_files(prefix = 'fastas/', suffix = '.fasta')
	outstring = 'Now blasting {:' + str(len(str(len(file_list)))) + 'd}/{}: {}'

	for i, filename in enumerate(file_list):
		print(outstring.format(i+1, len(file_list), filename))
		outfilename = f'blastresults/{_get_basename(filename)}.xml'
		phylo.run_blastp(filename, outfilename, blastdb)


def parse_blast_results(to_parse = None, to_exclude=None):
	'''
	Parses Blast results to a table (tsv file).

	:param to_parse: Should be an iterable with filenames that shall be
		parsed. If all files in `blastresults` shall be parsed, this must be None.
	:param to_exclude: Set of taxonomy ids to exclude from the results.
		Leave empty (or set to None) to not exclude any species.
	:creates: `resulttables/*.tsv`
	'''

	if to_parse is None:
		to_parse = CR.get_protein_files(prefix = 'blastresults/', suffix = '.xml')

	if to_exclude is None:
		to_exclude = set()

	os.makedirs('resulttables', exist_ok=True)

	outstring = 'Now parsing {:' + str(len(str(len(to_parse)))) + 'd}/{}: {:<30}'

	header = '\t'.join((
		'Tax-ID',
		'Acc',
		'Species',
		'Rank',
		'e-value',
		'Length',
		'Lineage',
		'Prot-Name',
		'Query-Protein',
	))

	for i, filename in enumerate(to_parse):
		basename = _get_basename(filename)
		print(outstring.format(i+1, len(to_parse), basename), end='\r')
		parsed_result = phylo.parse_blast_result(filename, TF = TF, top = 0, exclude=to_exclude)

		with open(f'resulttables/{basename}.tsv', 'w') as out:
			out.write(header)
			out.write('\n')
			out.write('\n'.join(parsed_result))


def combine_parsed_results():
	'''
	Combines parsed Blast results according on proteinlist.txt.

	:creates: `combinedtables/*.tsv`
	'''

	os.makedirs('combinedtables', exist_ok=True)

	limits = CR.get_limits()
	proteins = CR.get_protein_dict(prefix = 'resulttables/', suffix = '.tsv')

	for protein in sorted(proteins):
		print(f'combining tsv files... {protein:<40}', end='\r')
		outfn = f'combinedtables/{protein}.tsv'
		try:
			max_evalue, min_length = limits[protein]
		except KeyError:
			max_evalue, min_length = limits['default']

		combined = phylo.combine_parsed_results(proteins[protein], max_evalue, min_length)

		if combined:
			open(outfn, 'w').write(combined)


def tables_for_interactive_heatmap():
	'''
	Combines parsed Blast results for use for interactive heatmaps.

	:creates: `interactivetables/*.tsv`
	'''

	proteins = CR.get_protein_dict(prefix = 'resulttables/', suffix = '.tsv')

	os.makedirs('interactivetables', exist_ok=True)

	for protein in sorted(proteins):
		print(protein.ljust(50), end='\r')

		entries = phylo.table_for_interactive_heatmaps(proteins[protein], TF)

		with open(f'interactivetables/{protein}.tsv', 'w') as outfile:
			for entry in sorted(entries):
				outfile.write(f'{entry}\t{entries[entry]}\n')


def unique_names():
	'''
	For each protein, creates lists with the species that possess this
	protein. Two lists are generated: One with unique species names and
	one with unique taxonomy ids. In addition a general.names and
	general.taxids are generated with all species that occur in the lists
	of all proteins.

	:creates: `names/*.names`
	:creates: `taxids/*.taxids`
	'''

	os.makedirs('names', exist_ok=True)
	os.makedirs('taxids', exist_ok=True)

	total_names = set()
	total_taxids = set()

	for filename in CR.get_protein_names():
		names, taxids = phylo.unique_names_and_taxids(f'combinedtables/{filename}.tsv')

		with open(f'names/{filename}.names', 'w') as outfile:
			for name in names:
				outfile.write(name + '\n')

		with open(f'taxids/{filename}.taxids', 'w') as outfile:
			for taxid in taxids:
				outfile.write(str(taxid) + '\n')

		total_names.update(names)
		total_taxids.update(taxids)

	with open('names/general.names', 'w') as outfile:
		for name in sorted(total_names):
			outfile.write(name + '\n')

	with open('taxids/general.taxids', 'w') as outfile:
		for taxid in sorted(total_taxids):
			outfile.write(str(taxid) + '\n')


def make_newick():
	'''
	Creates Newick trees for every taxonomy id list in the taxids folder.

	:creates: `trees/*.tre`
	'''

	todo = CR.get_protein_names(prefix = 'taxids/', suffix = '.taxids')
	todo.add('taxids/general.taxids')

	os.makedirs('trees', exist_ok=True)

	sanitizer = phylo.NodeSanitizer()

	for filename in todo:
		outfn = f'trees/{_get_basename(filename)}.tre'
		newick = phylo.make_newick(filename, sanitizer, TF)

		open(outfn, 'w').write(newick)

	sanitizer.print_bad_chars()


def tree_attributes():
	'''
	For each element in the general tree, determine the proteins
	(attributes) for this element. Each attribute will be a two-letter
	string. The meaning of each string will be written to `attributekeys.txt`.

	:uses: `trees/general.tre`
	:creates: `attributes.txt` (containing the keys for each protein)
	'''

	proteinlist = sorted(CR.get_protein_names())
	filenames = {name: f'trees/{name}.tre' for name in proteinlist}

	master_tree = open('trees/general.tre').read()

	keys, attributes = phylo.get_keys_and_attributes(proteinlist, filenames, master_tree)

	with open('attributes.txt', 'w') as out:
		for key in sorted(keys):
			out.write(f'{key} {keys[key]}\n')

		out.write('---\n')

		for key in sorted(attributes):
			out.write(f'{key}\t{"".join(attributes[key])}\n')


def make_histograms():
	'''
	For each seed protein, the distribution of hits is shown in a histogram.

	:creates: `histograms/*.png`
	'''

	protein_files = CR.get_protein_dict(prefix = 'fastas/', suffix = '.fasta')
	seed_lengths = {}

	for prot, files in protein_files.items():
		seed_lengths[prot] = []
		for filename in files:
			with open(filename) as fastafile:
				length = 0
				next(fastafile)
				for line in fastafile:
					length += len(line.rstrip())
			seed_lengths[prot].append(length)

	os.makedirs('histograms', exist_ok=True)

	colormap = plt.cm.get_cmap('rainbow')

	for protein in sorted(protein_files):
		infile = f'combinedtables/{protein}.tsv'
		print(f'Histogram: {protein:<50}', end='\r')

		fig = phylo.make_histogram(infile, seed_lengths[protein], colormap=colormap)

		fig.savefig(f'histograms/{protein}.png')


def show_blast_mapping():
	'''
	For each protein, create an overview over where the hits where mapped
	over the length of the protein.

	:creates: `blastmappings/*.png`
	'''

	os.makedirs('blastmappings', exist_ok=True)

	proteinnames = sorted(CR.get_protein_files())

	for proteinname in proteinnames:
		print(f'Mapping {proteinname:<50}', end='\r')

		query_length = 0
		with open(f'fastas/{proteinname}.fasta') as fastafile:
			next(fastafile)
			for line in fastafile:
				query_length += len(line.rstrip())

		filename = f'blastresults/{proteinname}.xml'

		image = phylo.show_blast_mapping(filename, query_length)

		image.save(f'blastmappings/{proteinname}.png')


def similarity_matrix(delimiter=','):
	'''
	Creates a similarity matrix of proteins. For each protein, the Blast
	hits are compared to each other protein. The lowest e-value of each
	protein-protein combination is saved. This gives an impression of
	how similar two given proteins are and how much the results overlap.

	:uses: combined parsed Blast results in `combinedtables`
	:creates: `matrix.csv`
	'''

	names = {name: f'combinedtables/{name}.tsv' for name in CR.get_protein_names()}
	sorted_names = sorted(names)

	res = phylo.similarity_matrix(names)

	with open('matrix.csv', 'w') as out:
		out.write(f'{delimiter}{delimiter.join(sorted_names)}\n')
		for i, line in enumerate(res):
			out.write(f'{sorted_names[i]}{delimiter}{delimiter.join(map(str, line))}\n')


def int_heatmap():
	'''
	Creates an interactive heatmap as html file with javascript.

	:uses: heatmap_template.html
	:creates: heatmap.html
	'''

	taxa_to_check = CR.get_heatmap_taxa()
	colors = CR.get_heatmap_colors()
	method = CR.get_heatmap_clustering()
	proteins_to_check = sorted(CR.get_protein_names())
	template = open('heatmap_template.html').read()

	taxids = {}
	tick_taxa = []
	for i, tax in enumerate(taxa_to_check):
		taxname, taxid = tax.split('^', maxsplit=1)
		tick_taxa.append(taxname.replace('_', ' '))
		taxids[int(taxid.rstrip())] = i

	tick_proteins = proteins_to_check[:]

	lineages = {}
	with open('taxids/general.taxids') as all_taxids:
		for taxid in all_taxids:
			lineages[int(taxid.rstrip())] = TF.get_lineage_fast(int(taxid))

	todos = {}
	for taxid in taxids:
		if TF.get_tax_info(taxid)['rank'] == 'species':
			todos[taxid] = [taxid]
			continue
		todos[taxid] = []
		for basic_taxid, lineage in lineages.items():
			if taxid in lineage:
				todos[taxid].append(basic_taxid)

	matrix = []
	for i, protein in enumerate(proteins_to_check):
		matrix.append([-100] * len(taxa_to_check))
		taxid2evalue = {}
		with open(f'interactivetables/{protein}.tsv') as tablefile:
			for line in tablefile:
				lline = line.split(maxsplit=1)
				taxid2evalue[int(lline[0])] = int(lline[1])

		for taxid, pos in taxids.items():
			best = 0
			for todo in todos[taxid]:
				try:
					evalue = taxid2evalue[todo]
				except KeyError:
					# If the key is not found, the species is not in
					# the results for this protein, i.e. no evalue
					continue
				if evalue > best:
					best = evalue
			matrix[i][pos] = best

	html = phylo.interactive_heatmap(matrix, tick_taxa, tick_proteins, colors, template, method)

	open('heatmap.html', 'w').write(html)


tasks = {
	'blast': ('run Blast', run_blast),
	'parse': ('parse Blast results', parse_blast_results),
	'combine': ('combine parsed results', combine_parsed_results),
	'comheat': ('combine parsed results for heatmaps', tables_for_interactive_heatmap),
	'unique': ('create unique lists of names and taxids', unique_names),
	'newick': ('create Newick tree for each protein', make_newick),
	'attrib': ('determine tree attributes', tree_attributes),
	'hist': ('create histograms with Blast hits for each protein', make_histograms),
	'map': ('create hit mapping diagrams for each protein', show_blast_mapping),
	'intheat': ('create an interactive heatmap (html)', int_heatmap),
	'matrix': ('create a similarity matrix of all proteins', similarity_matrix),
}

tasknames = list(tasks)


def run_workflow(start, end=''):
	'''
	Starts the workflow from `start` until `end`. If `end` is empty, the
	workflow is run until the last element of the workflow.

	:param start: The name of the first element to run.
	:param end: The name of the last element to run or a falsy value if
		the workflow shall be run until the end.
	'''

	if start not in tasknames:
		raise ValueError(f'{start} is no valid task.')

	if not end:
		endidx = len(tasknames)
	elif end not in tasknames:
		raise ValueError(f'{end} is no valid task.')
	else:
		endidx = tasknames.index(end) + 1

	startidx = tasknames.index(start)
	for taskname in tasknames[startidx:endidx]:
		print(f'{taskname}: "{tasks[taskname][0]}"')
		task = tasks[taskname][1]
		task()


def main():
	'''
	This is the start point for the command line tool.
	'''

	global TF, CR, blastdb

	def to_set(string):
		return {int(x) for x in string.split(',')}

	workflow = '\n'.join(textwrap.wrap('The following is a list of the '
	'workflow. The names or numbers can be used for the -s or -o '
	'arguments.', width = 80))

	workflow += '\n\n' + '\n'.join(
		(f'{i:>2}. {name:<8} {tasks[name][0]}' for i, name in enumerate(tasknames))
	)

	parser = argparse.ArgumentParser(
		description='This module provides you with tools to run '
		'phylogenetic analyses. Exactly one argument must be given.')

	parser.add_argument('-l', '--list', action='store_true',
		help='Shows the whole workflow with information and exits')

	parser.add_argument('-i', '--init', action='store_true',
		help='Initiates the working directory with necessary files and folders')

	parser.add_argument('-a', '--all', action='store_true',
		help='Run the full workflow without Blast')

	parser.add_argument('-b', '--blast', action='store_true',
		help='Run the full workflow including Blast')

	parser.add_argument('-s', '--startfrom', default='',
		help='Run from and including this step [e.g. 7 or hist]')

	parser.add_argument('-o', '--only', default='',
		help='Run only the given step [e.g. 4 or unique]')

	parser.add_argument('-d', '--database', default='',
		help='Path to the Blast database to use. Only needed when '
		'actually running Blast (-b or -[o|s] blast)')

	parser.add_argument('-e', '--exclude', type=to_set, default={118797, 59538, 7213},
		help='Comma-separated list of taxonomy ids to exclude. By default '
		'the species Lipotes vexillifer (118797), Pantholops hodgsonii '
		'(59538), and Ceratitis capitata (7213) are excluded due to '
		'known massive contamination of bacterial sequences in the genome.')

	args = parser.parse_args()

	if args.list:
		parser.print_help()
		print('')
		print(workflow)
		sys.exit()

	if args.init:
		init()
		sys.exit()

	num_arguments = args.all + args.blast + bool(args.startfrom) + bool(args.only)

	if not (num_arguments == 1 or (num_arguments == 2 and args.database != '')):
		parser.print_help()
		sys.exit()

	blastdb = args.database

	# We need objects of these two classes for most of the functions, so
	# we initialize them here already. TaxFinder takes some seconds to
	# load, so this is, what makes loading this module slow.
	TF = TaxFinder()
	try:
		CR = ConfigReader()
	except IOError:
		CR = None

	if args.all:
		run_workflow('parse')
	elif args.blast:
		run_workflow('blast')
	elif args.startfrom:
		try:
			todo_index = int(args.startfrom)
		except ValueError:
			try:
				todo_index = tasknames.index(args.startfrom)
			except ValueError:
				parser.print_help()
		if todo_index < len(tasknames):
			run_workflow(tasknames[todo_index])
		else:
			parser.print_help()
	elif args.only:
		try:
			todo_index = int(args.only)
		except ValueError:
			try:
				todo_index = tasknames.index(args.only)
			except ValueError:
				parser.print_help()
		if todo_index < len(tasknames):
			task = tasks[tasknames[todo_index]][1]
			task()
		else:
			parser.print_help()
	else:
		print('This should not happen!')
		parser.print_help()
