#!/usr/bin/env python3

import sys
import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy.cluster.hierarchy as hierarchy
from PIL import Image, ImageDraw, ImageFont
from Bio.Blast import NCBIXML
from Bio.Blast.Applications import NcbiblastpCommandline
from multiprocessing import cpu_count


class NodeSanitizer():
	'''
	The NodeSanitizer is needed to create Newick files. It filters out bad characters and replaces them with allowed, similar characters. If there were unknown characters, they can be printed out.
	'''

	def __init__(self):
		'''
		Initiates the class.
		'''

		self.bad_chars = set()
		self.good_chars = {'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z', ' ', 'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm', 'n', 'o', 'p', 'q', 'r', 's', 't', 'u', 'v', 'w', 'x', 'y', 'z', '0', '1', '2', '3', '4', '5', '6', '7', '8', '9', '^', '_', '=', '-', '/', '.', '*'}
		self.to_replace = {'(': '<', ')': '>', '[': '<', ']': '>', '#': '_', ':': '_', '+': '_', "'": '_', ',': '_'}


	def sanitize(self, nodes):
		'''
		Goes through a list of nodes and replaces bad characters.

		:param nodes: The list of nodes (list of strings) to be sanitized
		:returns: A list with sanitized nodes (strings)
		'''

		result = []

		for node in nodes:
			new_name = []
			for char in node:
				if char not in self.good_chars:
					if char in self.to_replace:
						new_name.append(self.to_replace[char])
					else:
						self.bad_chars.add(char)
						new_name.append('!')
				else:
					new_name.append(char)
			result.append(''.join(new_name))

		return result


	def print_bad_chars(self):
		'''
		If unknown characters were found in the last call(s) of `NodeSanitizer.sanitize()`, these characters are printed. Otherwise, nothing happens.
		'''

		if self.bad_chars:
			print('Unknown chars found:')
			for elem in sorted(self.bad_chars):
				print(elem)


def run_blastp(query, outfilename, db, evalue = 1, maxthreads = cpu_count()):
	'''
	Run Blastp. This may take a long time (several minutes up to half-an-hour depending on the computer power and the size of the query and the database.

	:param query: The filename of the query sequence (may include a path)
	:param outfilename: The filename (including path) to save the result to
	:param db: The database file
	:param evalue: The e-value cut-off (should usually be very high, as the results will be filtered out afterwards)
	:param maxthreads: The maximum number of threads to use by Blast
	:creates: `outfilename`
	'''

	stdout, stderr = NcbiblastpCommandline(query = query, db = db, evalue = evalue, outfmt = 5, out = outfilename, num_threads = maxthreads, max_target_seqs = 20000)

	if stdout:
		print(stdout, file=sys.stderr)
	if stderr:
		print(stderr, file=sys.stderr)


def parse_blast_result(blast_XML, TF, top = 0, exclude=None, new_header=True):
	'''
	Parses Blast result XML files and writes the best or all results with less information in a tsv file.

	:param blast_XML: The filename of the Blast output (must be Blast output type 5)
	:param TF: An instance of the TaxFinder class
	:param top: Return only the best `top` hits. If `top` is 0, all hits are returned.
	:param exclude: Set with taxids of species to exclude from the results
	:param new_header: Were the Blast results produced with new headers (database from 2016 and newer)?
	:returns: tsv table as string with the results
	'''

	if exclude is None:
		exclude = set()

	if top < 0:
		top = 0

	result = []

	with open(blast_XML) as f:
		records = NCBIXML.parse(f)

		result.append('\t'.join(('Tax-ID', 'Acc', 'Species', 'Rank', 'e-value', 'Length', 'Lineage', 'Prot-Name', 'Query-Protein')))

		for record in records:
			for i, alignment in enumerate(record.alignments):
				if top and i > top:
					break

				infos = TF.getInfoFromHitDef(alignment.hit_id, alignment.hit_def, newHeader = new_header)

				for info in infos:
					if info['taxid'] in exclude:
						continue

					lineage = ', '.join(TF.getLineage(info['taxid'], display = 'name'))

					for hsp in alignment.hsps:
						try:
							line = '\t'.join((str(info['taxid']), info['acc'], info['name'], info['rank'], str(hsp.expect), str(hsp.align_length), lineage, info['protname'], record.query.split('|')[1]))
						except IndexError:
							line = '\t'.join((str(info['taxid']), info['acc'], info['name'], info['rank'], str(hsp.expect), str(hsp.align_length), lineage, info['protname'], record.query))

						result.append(line)

	return '\n'.join(result)


def combine_parsed_results(parsed_results, max_evalue, min_length):
	'''
	Combine parsed results to a tsv table.

	:param parsed_results: List of filenames with parsed blast results to combine
	:param max_evalue: Highest e-value to include (float)
	:param min_length: Minimal length to include (int)
	:returns: String of a tsv with the combined table
	'''

	header = True
	results = []
	for filename in parsed_results:
		with open(filename) as f:
			# Only copy the header once, no matter how many files are parsed
			if header:
				results.append(next(f).rstrip())
				header = False
			else:
				next(f)

			for line in f:
				lline = line.split('\t')
				if float(lline[4]) <= max_evalue and int(lline[5]) >= min_length:
					results.append(line.rstrip())

	return '\n'.join(results)


def table_for_interactive_heatmaps(parsed_results, TF):
	'''
	Combine parsed results to a table for creation of interactive tables.

	:param parsed_results: List of filenames with parsed blast results to combine
	:param TF: Instance of taxfinder.TaxFinder
	:returns: dict with a mapping from taxonomy id to evalue
	'''

	entries = {}
	for filename in parsed_results:
		with open(filename) as f:
			next(f)
			for line in f:
				lline = line.split('\t')
				taxid = lline[0]
				rank = lline[3]
				evalue = lline[4]

				if rank != 'species':
					try:
						taxid = str(TF.getSpeciesFromSubspecies(taxid))
					except ValueError:
						continue

				if taxid not in entries:
					entries[taxid] = -1

				if evalue == '0.0':
					e = 200
				elif 'e' in evalue:
					e = -1 * int(evalue.split('e')[1])
				else:
					e = math.ceil(math.log10(float(evalue)))

				if e > entries[taxid]:
					entries[taxid] = e

	return entries


def unique_names_and_taxids(comb_table):
	'''
	Extract all names and taxonomy ids from a combined table and return them as lists with unique entries.

	:param comb_table: Filename of a combined table to extract names and taxids from
	:returns: Tuple of two lists; the unique names (strs) and taxonomy ids (ints)
	'''

	names = set()
	taxids = set()

	with open(comb_table) as f:
		next(f)
		for line in f:
			elems = line.split('\t')
			names.add(elems[2])
			taxids.add(int(elems[0]))

	names = sorted(names)
	taxids = sorted(taxids)

	return names, taxids


def make_newick(filename, sanitizer, TF):
	'''
	Creates a Newick tree from a list of taxonomy ids. The relationships of the given taxonomy ids are derived automatically using the taxfinder.TaxFinder.

	:param filename: Filename of a file with one taxonomy id per line of the taxa to include in the tree
	:param sanitizer: Instance of phylogenetics.NodeSanitizer
	:param TF: Instance of taxfinder.TaxFinder
	:returns: String with the Newick representation of the tree
	'''

	with open(filename) as f:
		lineages = []

		for line in f:
			lineages.append(TF.getLineage(int(line), display='both'))

		tree = {}

		for line in lineages:
			if line and line[0] not in tree:
				tree[line[0]] = set()

			for i in range(1, len(line)):
				if line[i] not in tree:
					tree[line[i]] = set()
				tree[line[i-1]].add(line[i])

		newick = '(root^1);'
		nodes = ['root^1']
		while nodes:
			node = nodes.pop(0)
			newnodes = tree.pop(node)
			if newnodes:
				sanitized_nodes = sanitizer.sanitize(newnodes)
				newick = newick.replace(node, '(' + ','.join(sanitized_nodes) + ')' + node)
				nodes.extend(newnodes)

	return newick


def _get_tree_elements(tree, return_set = False, splitting = True):
	'''
	Extract all nodes from a Newick tree.

	:param tree: The Newick tree as a string
	:param return_set: Return the tree elements as a set (True) or a dictionary where each element points to an empty list (False)?
	:param splitting: If the Newick tree elements are in the form `Name^Taxid`, this must be True. If they are in the form `taxid`, this must be False.
	:returns: Either a list with taxonomy ids or a dictionary with taxonomy ids as keys and empty lists as value, depending on `return_set`
	'''

	tree = tree.replace('(', '\n').replace(')', '\n').replace(',', '\n').replace(';', '')
	treelist = tree.split('\n')

	elements = {}

	for line in treelist:
		line = line.rstrip()
		if not line:
			continue
		if splitting:
			try:
				line = line.split('^')[1]
			except IndexError:
				print(line)
				raise

		elements[int(line)] = []

	if return_set:
		return set(elements)
	else:
		return elements


def _get_code(number):
	'''
	Returns a two-character code depending on the number. The form is: 1 -> aA, 2 -> aB, 26 -> aZ, 27 -> bA, ...

	:param number: The number to create the code from
	:returns: The two-character code
	'''

	# 97 = ord('a'); 65 = ord('A')
	return chr((int(number/26) % 26) + 97) + chr((number % 26) + 65)


def get_keys_and_attributes(proteinlist, treefiles, master_tree):
	'''
	Create a protein/two-letter-key table and a table with the presence of proteins in taxa.

	:param proteinlist: List of proteins to look for
	:param treefiles: Dict of `protein name`: `filenames of trees (Newick!)` two check for presence of taxonomy ids
	:param master_tree: The tree to use for extraction of taxonomy ids
	:returns: A tuple of two dicts; the two-letter-code to protein name and taxonomy id to two-letter-code
	'''

	attributes = _get_tree_elements(master_tree, return_set = False, splitting = True)

	keys = {}

	for i, name in enumerate(proteinlist):
		code = _get_code(i)
		keys[name] = code

	for protein in treefiles:
		taxids = _get_tree_elements(open(treefiles[protein]).read(), return_set = True, splitting = True)
		for taxid in taxids:
			if taxid in attributes:
				attributes[taxid].append(keys[protein])

	return keys, attributes


def make_histogram(combined_table, seed_length, width = 12, height = 6, colormap = None):
	'''
	Create a histogram showing the lengths of the blast hits for the query protein, helping to identify length cut-offs for Blast results.

	:param combined_table: Filename of the combined table of the query protein
	:param seed_length: Length of the query protein
	:param width: Width of the output figure (in inch)
	:param height: Height of the output figure (in inch)
	:param colormap: Instance of a matplotlib.cm.ScalarMappable (colormap). If this is None, the `rainbow` colormap will be used.
	:returns: Instance of matplotlib.pyplot.figure with the histogram
	'''

	if colormap is None:
		colormap = plt.cm.get_cmap('rainbow')

	values = np.loadtxt(combined_table, dtype=np.int, comments=None, delimiter='\t', skiprows=1, usecols=(5,))

	mi = np.amin(values)
	ma = np.amax(values)

	if ma - mi <= 1000:
		interval = 10
	elif ma - mi <= 2000:
		interval = 20
	elif ma - mi <= 2500:
		interval = 25
	elif ma - mi <= 5000:
		interval = 50
	else:
		interval = 100

	text = '''Distribution
min: {}
max: {}
average: {:.0f}
median: {:.0f}
total elements: {}

interval: {}'''.format(mi, ma, np.mean(values), np.median(values), values.size, interval)

	seeds = seed_length#s[protein]

	sizes = '''Seed protein(s)
min: {}
max: {}
average: {:.0f}
total seeds: {}'''.format(min(seeds), max(seeds), np.average(seeds), len(seeds))

	middle = int(ma/2 + mi/2)
	middle -= int(middle % interval)
	if middle - 50*interval < 0:
		middle = 50*interval

	bins = list(range(middle - 50*interval, middle + interval * 50 + 1, interval))

	plt.close()
	fig = plt.figure(1, figsize=(width, height))
	ax = fig.add_subplot(1,1,1)
	n, bins, patches = ax.hist(values, bins=bins)

	# The following block is to color the bars
	bin_centers = 0.5 * (bins[:-1] + bins[1:])
	col = bin_centers - min(bin_centers)
	col /= max(col)
	for c, p in zip(col, patches):
		plt.setp(p, 'facecolor', colormap(c))

	ax.text(0.05, 0.95, text, transform=ax.transAxes, horizontalalignment='left', verticalalignment='top')

	ax.text(0.95, 0.95, sizes, transform=ax.transAxes, horizontalalignment='right', verticalalignment='top')

	return fig


def show_blast_mapping(blast_result_file, query_length):
	'''
	Create an overview over where the Blast hits are mapped on the query protein.

	:param blast_result_file: The filename of a blast result (output 5!)
	:param query_length: Length of the query protein
	:returns: Instance of PIL.Image with the image
	'''

	fnt = ImageFont.load_default()

	counters = [np.zeros(query_length, np.int) for x in range(6)]
	num_hsps = [0] * 6

	with open(blast_result_file) as f:
		records = NCBIXML.parse(f)

		for record in records:
			for alignment in record.alignments:
				for hsp in alignment.hsps:
					if hsp.expect > 1e-15:
						n = 0
					elif hsp.expect > 1e-30:
						n = 1
					elif hsp.expect > 1e-60:
						n = 2
					elif hsp.expect > 1e-90:
						n = 3
					elif hsp.expect > 1e-120:
						n = 4
					else:
						n = 5
					counters[n][hsp.query_start - 1:hsp.query_end - 1] += 1
					num_hsps[n] += 1

	ma = [np.amax(counters[n]) * 0.01 for n in range(6)]

	counters = [counters[n] / ma[n] if ma[n] != 0 else np.ones(query_length, np.int) for n in range(6)]

	im = Image.new('RGB', (query_length + 60, 600), (255, 255, 255))
	dr = ImageDraw.Draw(im)

	dr.text((2, 40), '> 1e-15', (0, 0, 0), fnt)
	dr.text((2, 140), '> 1e-30', (0, 0, 0), fnt)
	dr.text((2, 240), '> 1e-60', (0, 0, 0), fnt)
	dr.text((2, 340), '> 1e-90', (0, 0, 0), fnt)
	dr.text((2, 440), '> 1e-120', (0, 0, 0), fnt)
	dr.text((2, 540), '<= 1e-120', (0, 0, 0), fnt)

	for n in range(6):
		dr.text((2, 60 + 100 * n), 'n = {}'.format(num_hsps[n]), (0, 0, 0), fnt)

	colors = [(0, 0, 0), (0, 0, 200), (0, 200, 0), (200, 0, 200), (200, 0, 0), (150, 150, 0)]

	for n in range(int(query_length / 100)):
		col = 160 + n*100
		dr.line([(col, 0), (col, 600)], fill=(125, 125, 125), width=1)

	for n in range(6):
		for col, thickness in enumerate(counters[n]):
			dr.line([(col + 60, n*100), (col + 60, thickness + n*100)], fill=colors[n], width=1)

	return im


def interactive_heatmap(matrix, tick_taxa, tick_proteins, colors, template, method):
	'''
	Create an interactive heatmap with HTML/JavaScript showing in which species proteins are found.

	:param matrix: List of lists with integers indication the -log(evalue) of a protein in a taxon. The first level (`matrix[x]`) should fit to `tick_proteins` and the second level (`matrix[…][x]`) should fit to `tick_taxa` and contain the -log(evalue) of that taxon to the protein.
	:param tick_taxa: List of taxa as strings
	:param tick_proteins: List of proteins as strings
	:param colors: Dict with five elements, mapping letters 'grcbm' (green, red, cyan, blue, magenta) to HTML color codes.
	:param template: HTML template to use for the output.
	:param method: Which clustering method to use (str). See scipy.cluster.hierarchy.linkage for options.
	:returns: HTML as string
	'''

	pdmatrix = pd.DataFrame(matrix, columns = tick_taxa, index = tick_proteins)

	linkage = hierarchy.linkage(pdmatrix, method=method)
	dendro = hierarchy.dendrogram(linkage, labels = tick_proteins, no_plot = True, distance_sort = True)

	data = []
	for num in dendro['leaves']:
		data.append(matrix[num])

	cl = [colors[c] for c in dendro['color_list']]

	xvalues = [x[:] for x in dendro['icoord']]
	yvalues = [y[:] for y in dendro['dcoord']]
	maxX = max((max(x) for x in xvalues))
	maxY = max((max(y) for y in yvalues))
	xvalues = [[x/maxX for x in a] for a in xvalues]
	yvalues = [[y/maxY for y in a] for a in yvalues]

	longest_spec_name = max(len(name) for name in tick_taxa)
	longest_prot_name = max(len(name) for name in tick_proteins)

	width = str(int(10 + 12 * len(tick_proteins) + 6.5 * longest_spec_name))
	height = str(int(75 + 10 + 12 * len(tick_taxa) + 7 * longest_prot_name))

	clusterp = [[cl[i]] + list(zip(*x)) for i, x in enumerate(zip(xvalues, yvalues))]

	cproteins_printable = repr(dendro['ivl'])
	cluster_printable = repr(clusterp).replace('(', '[').replace(')', ']')
	cdata_printable = repr(list(zip(*data))).replace('(', '[').replace(')', ']')
	taxa_printable = repr(tick_taxa)
	adata_printable = repr(list(zip(*matrix))).replace('(', '[').replace(')', ']')
	aproteins_printable = repr(tick_proteins)

	html = template.format(CWIDTH=width, CHEIGHT=height, CDATA=cdata_printable, TAXA=taxa_printable, CPROTEINS=cproteins_printable, CLUSTER=cluster_printable, ADATA=adata_printable, APROTEINS=aproteins_printable)

	return html


def similarity_matrix(names):
	'''
	Create a similarity matrix showing the best evalue of the best hit that is shared by two query proteins.

	:param names: Dict mapping the protein name to a filename
	:returns: The similarity matrix as tsv string
	'''

	sorted_names = sorted(names)

	values = {}
	for name in sorted_names:
		values[name] = [set() for _ in range(151)]
		with open(names[name]) as f:
			next(f)
			for line in f:
				lline = line.split('\t')
				evalue = lline[4]
				if 'e-' in evalue:
					evalue = int(evalue.split('e-')[1])
					if evalue > 150:
						evalue = 150
				elif evalue == '0.0':
					evalue = 150
				else:
					evalue = 0
				acc = lline[1]
				values[name][evalue].add(acc)

	res = []

	names = sorted(values)

	for i, name1 in enumerate(sorted_names):
		res.append([])
		for j, name2 in enumerate(sorted_names):
			if name1 == name2:	# this will give max anyway
				res[-1].append('-')
				continue
			if i > j:			# we save half of the calculation by using the symmetry
				res[-1].append(res[j][i])
				continue
			acc1 = set()
			acc2 = set()
			for evalue in range(150, -1, -1):
				acc1.update(values[name1][evalue])
				acc2.update(values[name2][evalue])
				inter = acc1.intersection(acc2)
				if len(inter) > 0.05 * (len(acc1) + len(acc2)):
					res[-1].append(evalue)
					break
			else:
				res[-1].append(0)

	return res
