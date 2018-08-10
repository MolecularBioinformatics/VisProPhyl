#!/usr/bin/env python3

import sys
import math
import pandas as pd
import scipy.cluster.hierarchy as hierarchy
from Bio.Blast import NCBIXML
from Bio.Blast.Applications import NcbiblastpCommandline
from multiprocessing import cpu_count


class NodeSanitizer():
	'''
	The NoteSanitizer is needed to create Newick files. It filters out bad characters and replaces them with allowed, similar characters. If there were unknown characters, they are printed out.
	'''

	def __init__(self):
		'''
		Initiates the class.
		'''

		self.bad_chars = set()
		self.good_chars = {'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z', ' ', 'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm', 'n', 'o', 'p', 'q', 'r', 's', 't', 'u', 'v', 'w', 'x', 'y', 'z', '0', '1', '2', '3', '4', '5', '6', '7', '8', '9', '^', '_', '=', '-', '/', '.', '*'}
		self.to_replace = {'(': '<', ')': '>', '[': '<', ']': '>', '#': '_', ':': '_', '+': '', "'": '', ',': '_'}


	def sanitize(self, nodes):
		'''
		Goes through a list of nodes and replace bad characters.

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

	:param blast_XML: The filename of the Blast output (must be output type 5)
	:param TF: An instance of the TaxFinder class
	:param top: Write only the best `top` hits to file. If `top` is 0, all hits are saved.
	:param contaminants: Set with taxids of species to exclude from the results
	:param new_header: Where the Blast results produced with new headers (database from 2016 and newer)?
	:returns: csv table as string with the results
	'''

	if contaminants is None:
		contaminants = set()

	if top < 0:
		top = 0

	result = []

	with open(blast_XML) as f, open(outfile, 'w') as out:
		records = NCBIXML.parse(f)

		result.append('\t'.join(('Tax-ID', 'Acc', 'Species', 'Rank', 'e-value', 'Length', 'Lineage', 'Prot-Name', 'Query-Protein')))

		for record in records:
			for i, alignment in enumerate(record.alignments):
				if top and i > top:
					break

				infos = TF.getInfoFromHitDef(alignment.hit_id, alignment.hit_def, newHeader = new_header)

				for info in infos:
					if info['taxid'] in contaminants:
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
	header = True
	results = []
	for fname in parsed_results:
		with open(fname) as f:
			# Only copy the header once, no matter how many files are parsed
			if header:
				results.append(next(f).rstrip())
				header = False
			else:
				next(f)

			for line in f:
				lline = line.split('\t')
				if float(lline[4]) <= maxevalue and int(lline[5]) >= minlength:
					results.append(line.rstrip())

	return '\n'.join(results)


def table_for_interactive_heatmaps(filenames, TF):
	entries = {}
	for filename in filenames:
		with open(fname) as f:
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


def unique_names_and_taxids(filename):
	names = set()
	taxids = set()

	with open('combinedtables/{}.tsv'.format(fn)) as f:
		next(f)
		for line in f:
			elems = line.split('\t')
			names.add(elems[2])
			taxids.add(int(elems[0]))

	names = sorted(names)
	taxids = sorted(taxids)

	return names, taxids


def make_newick(filename, sanitizer, TF):
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

	:param tree: The Newick tree as a string
	:param return_set: Shall the tree elements be returned as a set (True) or a dictionary where each element points to an empty list (False)?
	:param splitting: If the Newick tree elements are in the form `Name^Taxid`, this must be True. If it is in the form `taxid`, this must be False.
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

	return chr((int(number/26) % 26) + 97) + chr((number % 26) + 65)


def get_keys_and_attributes(proteinlist, filenames, master_tree):
	attributes = _get_tree_elements(master_tree, return_set = False, splitting = True)

	keys = []

	for i, name in enumerate(proteinlist):
		code = _get_code(i)
		keys[code] = name

	for filename in filenames:
		elements = _get_tree_elements(filename.read(), return_set = True, splitting = True)
		for element in elements:
			if element in attributes:
				attributes[element].append(code)

	return keys, attributes


def interactive_heatmap(matrix, tick_taxa, tick_proteins, colors, template, method):
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
	taxa_printable = repr(tick_taxons)
	adata_printable = repr(list(zip(*matrix))).replace('(', '[').replace(')', ']')
	aproteins_printable = repr(tick_proteins)

	html = template.format(CWIDTH=width, CHEIGHT=height, CDATA=cdata_printable, TAXA=taxa_printable, CPROTEINS=cproteins_printable, CLUSTER=cluster_printable, ADATA=adata_printable, APROTEINS=aproteins_printable)

	return html


def similarity_matrix(names):
	''' names is a dict: names[name] = filename '''

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
