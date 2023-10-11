#!/usr/bin/env python3

'''
This is the module for phylogenetics. This can be imported. If you want
to use it via the command line, see cli.py.
'''

import math
import os
import re
import sys

from functools import partial

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.cluster import hierarchy

from PIL import Image, ImageDraw, ImageFont
from Bio.Blast import NCBIXML
from Bio.Blast.Applications import NcbiblastpCommandline
from Bio import Entrez, SeqIO


class NodeSanitizer():
	'''
	The NodeSanitizer is needed to create Newick files. It filters out
	bad characters and replaces them with allowed, similar characters.
	If there were unknown characters, they can be printed out.
	'''

	def __init__(self):
		'''
		Initiates the class.
		'''

		self.bad_chars = set()
		self.good_chars = set('ABCDEFGHIJKLMNOPQRSTUVWXYZ abcdefghijklmnopqrstuvwxyz0123456789^_=-/.*')
		self.to_replace = {
			'(': '<',
			')': '>',
			'[': '<',
			']': '>',
			'#': '_',
			':': '_',
			'+': '_',
			"'": '_',
			',': '_',
		}


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
		If unknown characters were found in the last call(s) of
		`NodeSanitizer.sanitize()`, these characters are printed.
		Otherwise, nothing happens.
		'''

		if self.bad_chars:
			print('Unknown chars found:')
			for elem in sorted(self.bad_chars):
				print(elem)


def run_blastp(query, outfilename, blastdb, evalue=1, maxthreads=1, remote=True):
	'''
	Run Blastp. This may take a long time (several minutes up to
	half-an-hour depending on the computer power and the size of the
	query and the database.

	:param query: The filename of the query sequence (may include a path)
	:param outfilename: The filename (including path) to save the result to
	:param blastdb: The database file
	:param evalue: The e-value cut-off (should usually be very high, as
		the results will be filtered out afterwards)
	:param maxthreads: The maximum number of threads to use by Blast
	:param remote: If True, use remote Blast on the NCBI servers. Use local Blast otherwise.
	:creates: `outfilename`
	'''

	blast_cmd = partial(NcbiblastpCommandline, query=query, db=blastdb,
	evalue=evalue, outfmt=16, out=outfilename, max_target_seqs=20000)

	if remote:
		cmdline = blast_cmd(remote=True)
	else:
		cmdline = blast_cmd(num_threads=maxthreads)

	if cmdline:
		print(cmdline, file=sys.stderr)


def parse_blast_result(blast_xml, TF, top=0, exclude=None):
	'''
	Parses Blast result XML files and writes the best or all results with
	less information in a tsv file.

	:param blast_xml: The filename of the Blast output (must be Blast output type 16)
	:param TF: An instance of the TaxFinder class
	:param top: Return only the best `top` hits. If `top` is 0, all hits are returned.
	:param exclude: Set with taxids of species to exclude from the results
	:returns: tsv table as string with the results
	'''

	top = max(top, 0)

	result = []

	with open(blast_xml) as blast_result_file:
		records = NCBIXML.parse(blast_result_file)

		for record in records:

			try:
				query = record.query.split('|')[1]
			except IndexError:
				query = record.query

			for i, descr in enumerate(record.descriptions):
				if top and i > top:
					break

				evalue = descr.e

				alignment_length = record.alignments[i].hsps[0].align_length

				for hit in descr.items:
					taxid = hit.taxid
					lineage = TF.get_lineage(taxid, display = 'taxid')
					if exclude:
						stop = False
						for tid in exclude:
							if tid in lineage:
								stop = True
								break
						if stop:
							continue

					acc = hit.accession
					taxinfo = TF.get_tax_info(taxid)
					species = taxinfo['name']
					rank = taxinfo['rank']
					lineage_str = '>'.join(TF.get_lineage(taxid, display = 'name'))
					protname = hit.title.split('[', maxsplit=1)[0].rstrip()

					result.append(f'{taxid}\t{acc}\t{species}\t{rank}\t{evalue}\t{alignment_length}\t{lineage_str}\t{protname}\t{query}')

	return result


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
		with open(filename) as parsed_results_file:
			# Only copy the header once, no matter how many files are parsed
			if header:
				results.append(next(parsed_results_file).rstrip())
				header = False
			else:
				next(parsed_results_file)

			for line in parsed_results_file:
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
		with open(filename) as parsed_results_file:
			next(parsed_results_file)
			for line in parsed_results_file:
				lline = line.split('\t')
				taxid = lline[0]
				rank = lline[3]
				evalue = lline[4]

				if rank != 'species':
					try:
						taxid = str(TF.get_species_from_subspecies(taxid))
					except ValueError:
						continue

				if taxid not in entries:
					entries[taxid] = -1

				if evalue == '0.0':
					e_val = 200
				elif 'e' in evalue:
					e_val = -1 * int(evalue.split('e', maxsplit=1)[1])
				else:
					e_val = math.ceil(math.log10(float(evalue)))

				if e_val > entries[taxid]:
					entries[taxid] = e_val

	return entries


def unique_names_and_taxids(comb_table):
	'''
	Extract all names and taxonomy ids from a combined table and return
	them as lists with unique entries.

	:param comb_table: Filename of a combined table to extract names and taxids from
	:returns: Tuple of two lists; the unique names (strs) and taxonomy ids (ints)
	'''

	names = set()
	taxids = set()

	with open(comb_table) as comb_table_file:
		next(comb_table_file)
		for line in comb_table_file:
			elems = line.split('\t')
			names.add(elems[2])
			taxids.add(int(elems[0]))

	names = sorted(names)
	taxids = sorted(taxids)

	return names, taxids


def make_newick(filename, sanitizer, TF):
	'''
	Creates a Newick tree from a list of taxonomy ids. The relationships
	of the given taxonomy ids are derived automatically using the taxfinder.TaxFinder.

	:param filename: Filename of a file with one taxonomy id per line of
		the taxa to include in the tree
	:param sanitizer: Instance of phylogenetics.NodeSanitizer
	:param TF: Instance of taxfinder.TaxFinder
	:returns: String with the Newick representation of the tree
	'''

	with open(filename) as taxfile:
		lineages = []

		for line in taxfile:
			lineages.append(TF.get_lineage(int(line), display='both'))

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


def _get_tree_elements(tree, return_set=False, splitting=True):
	'''
	Extract all nodes from a Newick tree.

	:param tree: The Newick tree as a string
	:param return_set: Return the tree elements as a set (True) or a
		dictionary where each element points to an empty list (False)?
	:param splitting: If the Newick tree elements are in the form
		`Name^Taxid`, this must be True. If they are in the form `taxid`,
		this must be False.
	:returns: Either a list with taxonomy ids or a dictionary with
		taxonomy ids as keys and empty lists as value, depending on `return_set`
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
				line = line.split('^', maxsplit=1)[1]
			except IndexError:
				print(line)
				raise

		elements[int(line)] = []

	if return_set:
		return set(elements)

	return elements


def _get_code(number):
	'''
	Returns a two-character code depending on the number. The form is:
	1 -> aA, 2 -> aB, 26 -> aZ, 27 -> bA, ...

	:param number: The number to create the code from
	:returns: The two-character code
	'''

	# 97 = ord('a'); 65 = ord('A')
	return chr((int(number/26) % 26) + 97) + chr((number % 26) + 65)


def get_keys_and_attributes(proteinlist, treefiles, master_tree):
	'''
	Create a protein/two-letter-key table and a table with the presence
	of proteins in taxa.

	:param proteinlist: List of proteins to look for
	:param treefiles: Dict of `protein name`: `filenames of trees (Newick!)`
		two check for presence of taxonomy ids
	:param master_tree: The tree to use for extraction of taxonomy ids
	:returns: A tuple of two dicts; the two-letter-code to protein name
		and taxonomy id to two-letter-code
	'''

	attributes = _get_tree_elements(master_tree, return_set = False, splitting=True)

	keys = {}

	for num, name in enumerate(proteinlist):
		code = _get_code(num)
		keys[name] = code

	for protein in treefiles:
		taxids = _get_tree_elements(open(treefiles[protein]).read(), return_set=True, splitting=True)
		for taxid in taxids:
			if taxid in attributes:
				attributes[taxid].append(keys[protein])

	return keys, attributes


def make_histogram(combined_table, seed_length, width=12, height=6, colormap=None):
	'''
	Create a histogram showing the lengths of the blast hits for the
	query protein, helping to identify length cut-offs for Blast results.

	:param combined_table: Filename of the combined table of the query protein
	:param seed_length: Length of the query protein
	:param width: Width of the output figure (in inch)
	:param height: Height of the output figure (in inch)
	:param colormap: Instance of a matplotlib.cm.ScalarMappable (colormap).
		If this is None, the `rainbow` colormap will be used.
	:returns: Instance of matplotlib.pyplot.figure with the histogram
	'''

	if colormap is None:
		colormap = plt.cm.get_cmap('rainbow')

	values = np.loadtxt(combined_table, dtype=np.int, comments=None,
		delimiter='\t', skiprows=1, usecols=(5,))

	min_value = np.amin(values)
	max_value = np.amax(values)

	if max_value - min_value <= 1000:
		interval = 10
	elif max_value - min_value <= 2000:
		interval = 20
	elif max_value - min_value <= 2500:
		interval = 25
	elif max_value - min_value <= 5000:
		interval = 50
	else:
		interval = 100

	text = f'''Distribution
min: {min_value}
max: {max_value}
average: {np.mean(values):.0f}
median: {np.median(values):.0f}
total elements: {values.size}

interval: {interval}'''

	seeds = seed_length#s[protein]

	sizes = f'''Seed protein(s)
min: {min(seeds)}
max: {max(seeds)}
average: {np.average(seeds):.0f}
total seeds: {len(seeds)}'''

	middle = int(max_value/2 + min_value/2)
	middle -= int(middle % interval)
	if middle - 50*interval < 0:
		middle = 50*interval

	bins = list(range(middle - 50*interval, middle + interval * 50 + 1, interval))

	plt.close()
	fig = plt.figure(1, figsize=(width, height))
	axis = fig.add_subplot(1,1,1)
	_, bins, patches = axis.hist(values, bins=bins)

	# The following block is to color the bars
	bin_centers = 0.5 * (bins[:-1] + bins[1:])
	col = bin_centers - min(bin_centers)
	col /= max(col)
	for color, patch in zip(col, patches):
		plt.setp(patch, 'facecolor', colormap(color))

	axis.text(0.05, 0.95, text, transform=axis.transAxes,
		horizontalalignment='left', verticalalignment='top')

	axis.text(0.95, 0.95, sizes, transform=axis.transAxes,
		horizontalalignment='right', verticalalignment='top')

	return fig


def show_blast_mapping(blast_result_file, query_length):
	'''
	Create an overview over where the Blast hits are mapped on the query protein.

	:param blast_result_file: The filename of a blast result (output 16!)
	:param query_length: Length of the query protein
	:returns: Instance of PIL.Image with the image
	'''

	fnt = ImageFont.load_default()

	counters = [np.zeros(query_length, np.int) for x in range(6)]
	num_hsps = [0] * 6

	with open(blast_result_file) as blast_file_pointer:
		records = NCBIXML.parse(blast_file_pointer)

		for record in records:
			for alignment in record.alignments:
				for hsp in alignment.hsps:
					if hsp.expect > 1e-15:
						num = 0
					elif hsp.expect > 1e-30:
						num = 1
					elif hsp.expect > 1e-60:
						num = 2
					elif hsp.expect > 1e-90:
						num = 3
					elif hsp.expect > 1e-120:
						num = 4
					else:
						num = 5
					counters[num][hsp.query_start - 1:hsp.query_end - 1] += 1
					num_hsps[num] += 1

	max_e = [np.amax(counters[n]) * 0.01 for n in range(6)]

	counters = [
		counters[n] / max_e[n]
		if max_e[n] != 0
		else np.ones(query_length, np.int)
		for n in range(6)
	]

	image = Image.new('RGB', (query_length + 60, 600), (255, 255, 255))
	draw = ImageDraw.Draw(image)

	draw.text((2, 40), '> 1e-15', (0, 0, 0), fnt)
	draw.text((2, 140), '> 1e-30', (0, 0, 0), fnt)
	draw.text((2, 240), '> 1e-60', (0, 0, 0), fnt)
	draw.text((2, 340), '> 1e-90', (0, 0, 0), fnt)
	draw.text((2, 440), '> 1e-120', (0, 0, 0), fnt)
	draw.text((2, 540), '<= 1e-120', (0, 0, 0), fnt)

	for num in range(6):
		draw.text((2, 60 + 100 * num), f'n = {num_hsps[num]}', (0, 0, 0), fnt)

	colors = [(0, 0, 0), (0, 0, 200), (0, 200, 0), (200, 0, 200), (200, 0, 0), (150, 150, 0)]

	for num in range(int(query_length / 100)):
		col = 160 + num*100
		draw.line([(col, 0), (col, 600)], fill=(125, 125, 125), width=1)

	for num in range(6):
		for col, thickness in enumerate(counters[num]):
			draw.line([(col + 60, num*100), (col + 60, thickness + num*100)], fill=colors[num], width=1)

	return image


def interactive_heatmap(matrix, tick_taxa, tick_proteins, colors, template, method):
	'''
	Create an interactive heatmap with HTML/JavaScript showing in which
	species proteins are found.

	:param matrix: List of lists with integers indication the -log(evalue)
		of a protein in a taxon. The first level (`matrix[x]`) should fit
		to `tick_proteins` and the second level (`matrix[…][x]`) should
		fit to `tick_taxa` and contain the -log(evalue) ofthat taxon to
		the protein.
	:param tick_taxa: List of taxa as strings
	:param tick_proteins: List of proteins as strings
	:param colors: Dict with ten elements, mapping letters 'C0', 'C1', ... 'C9'
		to HTML color codes.
	:param template: HTML template to use for the output.
	:param method: Which clustering method to use (str). See
		scipy.cluster.hierarchy.linkage for options.
	:returns: HTML as string
	'''

	pdmatrix = pd.DataFrame(matrix, columns = tick_taxa, index = tick_proteins)

	linkage = hierarchy.linkage(pdmatrix, method=method)
	dendro = hierarchy.dendrogram(linkage, labels = tick_proteins,
		no_plot = True, distance_sort = True)

	data = []
	for num in dendro['leaves']:
		data.append(matrix[num])

	color_list = [colors[c] for c in dendro['color_list']]

	xvalues = [x[:] for x in dendro['icoord']]
	yvalues = [y[:] for y in dendro['dcoord']]
	max_x = max((max(x) for x in xvalues))
	max_y = max((max(y) for y in yvalues))
	xvalues = [[x/max_x for x in a] for a in xvalues]
	yvalues = [[y/max_y for y in a] for a in yvalues]

	longest_spec_name = max(len(name) for name in tick_taxa)
	longest_prot_name = max(len(name) for name in tick_proteins)

	width = str(int(10 + 12 * len(tick_proteins) + 6.5 * longest_spec_name))
	height = str(int(75 + 10 + 12 * len(tick_taxa) + 7 * longest_prot_name))

	clusterp = [[color_list[i]] + list(zip(*x)) for i, x in enumerate(zip(xvalues, yvalues))]

	cproteins_printable = repr(dendro['ivl'])
	cluster_printable = repr(clusterp).replace('(', '[').replace(')', ']')
	cdata_printable = repr(list(zip(*data))).replace('(', '[').replace(')', ']')
	taxa_printable = repr(tick_taxa)
	adata_printable = repr(list(zip(*matrix))).replace('(', '[').replace(')', ']')
	aproteins_printable = repr(tick_proteins)

	html = template.format(
		CWIDTH=width,
		CHEIGHT=height,
		CDATA=cdata_printable,
		TAXA=taxa_printable,
		CPROTEINS=cproteins_printable,
		CLUSTER=cluster_printable,
		ADATA=adata_printable,
		APROTEINS=aproteins_printable
	)

	return html


def similarity_matrix(names):
	'''
	Create a similarity matrix showing the best evalue of the best hit
	that is shared by two query proteins.

	:param names: Dict mapping the protein name to a filename
	:returns: The similarity matrix as tsv string
	'''

	sorted_names = sorted(names)

	values = {}
	for name in sorted_names:
		values[name] = [set() for _ in range(151)]
		with open(names[name]) as name_file:
			next(name_file)
			for line in name_file:
				lline = line.split('\t')
				evalue = lline[4]
				if 'e-' in evalue:
					evalue = min(int(evalue.split('e-', maxsplit=1)[1]), 150)
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
			if i > j:	# we save half of the calculation by using the symmetry
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


def get_entries_from_blast_result(filename, cutoff=1e-30, TF=None):
	'''
	Extract best hits from a Blast result

	:param filename: Filename of a Blast xml2 output (outfmt 16) to parse
	:param cutoff: e-value cutoff of hits to include
	:param TF: TaxFinder instance. Optional, but necessary to unify subspecies results
	:returns: Dict where taxids point to tuples with accession, e-value,
		scientific species name, start of the hit, end of the hit;
		entries[taxid] = (acc, evalue, sciname, hitfrom, hitto)
	'''

	entries = {}

	acc = None
	taxid = None
	sciname = None
	evalue = None
	hitfrom = None
	hitto = None

	# TODO: Isn't this something I should use Biopython for??

	with open(filename) as xml:
		for line in xml:
			if '<accession>' in line:
				acc = line.strip().replace('<accession>', '').replace('</accession>', '')

			elif not acc:
				continue

			elif '<taxid>' in line:
				taxid = int(line.strip().replace('<taxid>', '').replace('</taxid>', ''))

				try:
					taxid = TF.get_species_from_subspecies(taxid)
				except ValueError:
					print(f'Taxid {taxid} not found!', file=sys.stderr)
				except AttributeError:
					pass


			elif '<sciname>' in line:
				sciname = line.strip().replace('<sciname>', '').replace('</sciname>', '')

			elif '<evalue>' in line:
				evalue = float(line.strip().replace('<evalue>', '').replace('</evalue>', ''))

			elif '<hit-from>' in line:
				hitfrom = int(line.strip().replace('<hit-from>', '').replace('</hit-from>', ''))

			elif '<hit-to>' in line:
				hitto = int(line.strip().replace('<hit-to>', '').replace('</hit-to>', ''))

				if hitto < hitfrom:
					hitto, hitfrom = hitfrom, hitto

				if evalue < cutoff and (taxid not in entries or entries[taxid][1] > evalue):
					entries[taxid] = (acc, evalue, sciname, hitfrom, hitto)

				acc = None
				taxid = None
				sciname = None
				evalue = None
				hitfrom = None
				hitto = None

	return entries


def download_nucleotide_sequences(entries, mail, title=0, strip=False, logfile=None):
	'''
	Download nucleotide sequences from NCBI given a list of accessions.

	:param entries: List of strings that are accessions of the nucleotides to download
	:param mail: Your email address. Used by the NCBI in case something goes wrong
	:param title: (int) If greater than 0, shorten the title to this length
	:param strip: If True, trailing stop codons are stripped off from the sequenced
	:param logfile: Open file pointer stating where to write the log output.
		Can also be sys.stderr. If None, the logs will be discarded
	:returns: A string in fasta format ready to be written to disk
	'''

	if logfile is None:
		logfile = open(os.devnull, 'w')

	Entrez.email = mail
	out = []
	total = len(entries)

	for current, taxid in enumerate(entries, start=1):
		print('\r' + ' '*75, end='\r', flush=True, file=logfile)

		print(f'Running seq {current:3d} of {total:3d}: {entries[taxid][0]:<15}', end='', flush=True, file=logfile)

		print('dl:', end='', flush=True, file=logfile)
		try:
			handle = Entrez.efetch(db='nuccore', id=entries[taxid][0], rettype='gb', retmode='text')
		except IOError as err:
			print(f'\r{entries[taxid][0]}: {err}', file=logfile)
			continue

		print('ok parse:', end='', flush=True, file=logfile)
		seqobj = SeqIO.parse(handle, 'gb')
		record = next(seqobj)
		sequence = record._seq

		for feature in record.features:
			if feature.type == 'CDS':
				try:
					start = int(feature.location._start)
					end = int(feature.location._end)
				except (AttributeError, KeyError):
					continue

				if (not entries[taxid][3] and not entries[taxid][4]) or \
				(start <= entries[taxid][3] and end >= entries[taxid][4]):
					break
		else:
			print(f'\r{entries[taxid][0]}: No CDS{" "*40}', file=logfile)
			continue

		cds = str(sequence[start:end])

		if cds[:3] != 'ATG' or cds[-3:] not in ['TAA', 'TGA', 'TAG']:
			if cds[-3:] == 'CAT' and cds[:3] in ['CTA', 'TCA', 'TTA']:
				cds = str(sequence[start:end].reverse_complement())
			else:
				print(f'\r{entries[taxid][0]}: No ATG or Stop codon found! Sequence will be omitted{" "*30}', file=logfile)
				continue

		if len(cds) % 3:
			print(f'\r{entries[taxid][0]}: Possible frameshit! Sequence will be omitted{" "*40}', file=logfile)
			continue

		if re.search(r'[^ACGT]', cds):
			print(f'\r{entries[taxid][0]}: Non canonical DNA base! Sequence will be included in output.{" "*40}', file=logfile)

		if strip and cds[-3:] in ['TAA', 'TGA', 'TAG']:
			cds = cds[:-3]

		fasta_head = f'>{entries[taxid][0]}_{entries[taxid][2].replace(" ", "-")}'
		if title:
			fasta_head = fasta_head[:title]

		out.append(f'{fasta_head}\n{cds}')

	print('', file=logfile)

	return '\n'.join(out)


def download_protein_sequences(entries, mail, title=0, logfile=None):
	'''
	Download protein sequences from NCBI given a list of accessions.

	:param entries: List of strings that are accessions of the proteins to download
	:param mail: Your email address. Used by the NCBI in case something goes wrong
	:param title: (int) If greater than 0, shorten the title to this length
	:param logfile: Open file pointer stating where to write the log output.
		Can also be sys.stderr. If None, the logs will be discarded
	:returns: A string in fasta format ready to be written to disk
	'''

	if logfile is None:
		logfile = open(os.devnull, 'w')

	Entrez.email = mail
	out = []
	total = len(entries)

	for current, taxid in enumerate(entries, start=1):
		print('\r' + ' '*75, end='\r', flush=True, file=logfile)

		print(f'Running seq {current:3d} of {total:3d}: {entries[taxid][0]:<15}', end='', flush=True, file=logfile)

		print('dl:', end='', flush=True, file=logfile)
		try:
			handle = Entrez.efetch(db='protein', id=entries[taxid][0], rettype='fasta', retmode='text')
		except IOError as err:
			print(f'\r{entries[taxid][0]}: {err}', file=logfile)
			continue

		fasta = str(SeqIO.read(handle, 'fasta')._seq)
		handle.close()

		print('ok parse:', end='', flush=True, file=logfile)

		fasta_head = f'>{entries[taxid][0]}_{entries[taxid][2].replace(" ", "-")}'
		if title:
			fasta_head = fasta_head[:title]

		out.append(f'{fasta_head}\n{fasta}')

	print('', file=logfile)

	return '\n'.join(out)
