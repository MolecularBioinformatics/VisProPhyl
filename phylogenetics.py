#!/usr/bin/env python3

import sys
import os
import math
import numpy as np
import matplotlib.pyplot as plt
from glob import glob
from PIL import Image, ImageDraw, ImageFont
from Bio.Blast import NCBIXML
from Bio.Blast.Applications import NcbiblastpCommandline
from Bio import Entrez
from multiprocessing import cpu_count
from taxfinder import TaxFinder


class ConfigReader():
	'''
	The ConfigReader shall read the `limits.txt` and `proteinlist.txt` and return the contents in different formats upon request.
	'''

	def __init__(self):
		'''
		Initiates the class.
		'''

		self._readProteinlist()
		self._readLimits()


	def _readLimits(self):
		'''
		Reads `limits.txt` and saves the content as dictionary to be used by the `getLimits()` method.
		'''

		with open('limits.txt', 'r') as f:
			self.limitsDict = {}
			for line in f:
				line = line.split('#')[0].rstrip()
				if line == '':
					continue
				lline = line.split()
				try:
					evalue = float(lline[1])
				except ValueError:
					evalue = None
				try:
					length = int(lline[2])
				except ValueError:
					length = None
				self.limitsDict[lline[0]] = (evalue, length)


	def getLimits(self):
		'''
		Reads `limits.txt` and returns a dict in the format: {'protein': (evalue, length), ...} where evalue is the negative exponent of the evalue (int) and length is the length limit for the protein (int).

		:uses: `limits.txt`
		:returns: A dictionary with limits as described above
		'''

		return self.limitsDict


	def _readProteinlist(self):
		'''
		Reads `proteinlist.txt` and saves the content as lists to be used by the get* methods.
		'''

		with open('proteinlist.txt', 'r') as f:
			self.proteins = []
			self.proteinFiles = []
			for line in f:
				line = line.split('#')[0].strip()
				if not line:
					continue
				elems = line.split()
				if elems[0] not in self.proteins:
					self.proteins.append(elems[0])
					self.proteinFiles.append([])
				pidx = self.proteins.index(elem[0])
				self.proteinFiles[pidx].append(elems[1].replace('.fasta', ''))


	def getProteinNames(self):
		'''
		Returns a set of proteins in the format: {'protein1', 'protein2', ...}.

		:returns: A set with proteins
		'''

		return set(self.proteins)


	def getProteinFiles(self, prefix = '', suffix = ''):
		'''
		Returns a set of filenames in the format: {'proteinfile1', 'proteinfile2', ...}. The files are just the basenames without suffix. A prefix (a path) can be added as well as a suffix. getProteinFiles('abc/', '.txt') will turn 'file1' into 'abc/file1.txt'.

		:param prefix: A string to be prepended before the filename
		:param suffix: A string to be appended after the filename
		:returns: A set of filenames
		'''

		return set((prefix + fn + suffix for proteinlist in self.proteinFiles for fn in proteinlist))


	def getProteinDict(self, prefix = '', suffix = ''):
		'''
		Returns a dict with proteins and their respective file names in the format: {'protein': set('file1', 'file2', ...), ...}. The files are just the basenames without suffix. A prefix (a path) can be added as well as a suffix. getProteinDict('abc/', '.txt') will turn 'file1' into 'abc/file1.txt'.

		:param prefix: A string to be prepended before the filename
		:param suffix: A string to be appended after the filename
		:returns: A dictionary as described above
		'''

		ret = {}
		for i in range(len(self.proteins)):
			ret[self.proteins[i]] = set((prefix + fn + suffix for fn in self.proteinFiles[i]))

		return ret

	# end ConfigReader



# We need objects of these two classes for most of the functions, so we initialize them here already
# TaxFinder takes some seconds to load, so this is, what makes loading this module slow.
TF = TaxFinder()
CR = ConfigReader()

###### Here is a hardcoded path. That is not good!
def _runBlast(query, db = '/home/mathias/projects/nr/nr', evalue = 1, maxthreads = cpu_count()):
	'''
	Private method to run Blast. This may take a long time (several minutes up to half-an-hour depending on the computer power and the size of the query and the database.

	:param query: The filename of the query sequence (may include a path)
	:param db: The database file
	:param evalue: The e-value cut-off (should usually be very high, as the results will be filtered out afterwards)
	:param maxthreads: The maximum number of threads to use by Blast
	:creates: `blastresults/QUERY.xml`
	'''

	os.makedirs('blastresults', exist_ok=True)

	outfilename = 'blastresults/' + os.path.basename(query.replace('fasta', 'xml'))

	cmd = NcbiblastpCommandline(query = query, db = db, evalue = evalue, outfmt = 5, out = outfilename, num_threads = maxthreads, max_target_seqs = 20000)
	stdout, stderr = cmd()

	print('Blasting', fname, 'done')
	if stdout:
		print(stdout)
	if stderr:
		print(stderr, file=sys.stderr)
	else:
		print('No errors')


def runBlast():
	'''
	Blasts a list of files.
	'''

	li = glob('fastas/*.fasta')		### maybe get from proteinlist.txt
	outstring = 'Now blasting {:' + str(len(str(len(li)))) + 'd}/{}: {}'

	for i, fname in enumerate(li):
		print(outstring.format(i+1, len(li), fname))
		_runBlast(fname)


def _getTopFromBlast(blastXML, TF, top = 0, exContaminSpecies=True, outfile=None, newHeader=True):
	'''
	Parses Blast result XML files and writes the best or all results with less information in a tsv file.

	:param blastXML: The filename of the Blast output (must be output type 5)
	:param TF: An instance of the TaxFinder class
	:param top: Write only the best `top` hits to file. If `top` is 0, all hits are saved.
	:param exContaminSpecies: Shall hits of known contaminated species be excluded?
	:param outfile: The file to write the results to (including path). If it is None, use the basename of `blastXML`
	:param newHeader: Where the Blast results produced with new headers (database from 2016 and newer)?
	:creates: `resulttables/FILE.xml.tsv`
	'''

	contaminantSpecies = set((118797, 59538, 7213))	# Lipotes vexillifer, Pantholops hodgsonii, Ceratitis capitata

	if outfile is None:
		outfile = 'resulttables/{}.tsv'.format(os.path.basename(blastXML))

	if top < 0:
		top = 0

	with open(blastXML, 'r') as f, open(outfile, 'w') as out:
		records = NCBIXML.parse(f)

		out.write('\t'.join(('Tax-ID', 'Acc', 'Species', 'Rank', 'e-value', 'Length', 'Lineage', 'Prot-Name', 'Query-Protein')) + '\n')

		for record in records:
			for i, alignment in enumerate(record.alignments):
				if top and i > top:
					break

				infos = TF.getInfoFromHitDef(alignment.hit_id, alignment.hit_def, newHeader = newHeader)

				for info in infos:
					if exContaminSpecies and info['taxid'] in contaminantSpecies:
						continue

					lineage = ', '.join(TF.getLineage(info['taxid'], display = 'name'))

					for hsp in alignment.hsps:
						line = '\t'.join((str(info['taxid']), info['acc'], info['name'], info['rank'], str(hsp.expect), str(hsp.align_length), lineage, info['protname'], record.query.split('|')[1]))

						out.write(line + '\n')


#### This should get some argument to allow only certain files to be reparsed (as this takes time!)
def parseBlastResults():
	'''
	Parses Blast results.

	:creates: `resulttables/*.tsv`
	'''

	li = glob('blastresults/*.xml')	### This may be taken from the proteinlist.txt

	#li = ['blastresults/Athaliana_nadA_QS.xml']

	os.makedirs('resulttables', exist_ok=True)

	outstring = 'Now parsing {:' + str(len(str(len(li)))) + 'd}/{}: {:<30}'

	for i, fname in enumerate(li):
		basename = os.path.basename(fname)
		print(outstring.format(i+1, len(li), basename), end='\r')
		getTopFromBlast(fname, TF = TF, top = 0, outfile='resulttables/{}.tsv'.format(basename), exContaminSpecies=True)

	print('')


def combineParsedResults():
	'''
	Combines parsed Blast results.

	:creates: `combinedtables/*.tsv`
	'''

	os.makedirs('combinedtables', exist_ok=True)

	limits = CR.getLimits()
	proteins = CR.getProteinDict(prefix = 'resulttables/', suffix = '.tsv')

	for k in sorted(proteins):
		print('combining tsv files... {:<40}'.format(k), end='\r')
		outfn = 'combinedtables/{}.tsv'.format(k)
		try:
			maxevalue, minlength = limits[k]
			if maxevalue is None:
				maxevalue = limits['default'][0]
			if minlength is None:
				minlength = limits['default'][1]
		except KeyError:
			maxevalue, minlength = limits['default']
		header = True
		with open(outfn, 'w') as outfile:
			for fname in proteins[k]:
				with open(fname, 'r') as f:
					if header:
						outfile.write(next(f).rstrip() + '\n')
						header = False
					else:
						next(f)
					for line in f:
						lline = line.split('\t')
						if float(lline[4]) <= maxevalue and int(lline[5]) >= minlength:
							outfile.write(line.rstrip() + '\n')

		if not os.path.getsize(outfn):
			os.remove(outfn)

	print('')


def tablesForInteractiveHeatmap():
	'''
	Combines parsed Blast results for use for interactive heatmaps.

	:creates: `interactivetables/*.tsv`
	'''

	proteins = CR.getProteinDict(prefix = 'resulttables/', suffix = '.tsv')

	os.makedirs('interactivetables', exist_ok=True)

	totalNames = set()
	totalTaxids = set()

	for k in sorted(proteins):
		print(k.ljust(50), end='\r')

		entries = {}
		for fname in proteins[k]:
			with open(fname, 'r') as f:
				next(f)
				for line in f:
					lline = line.split('\t')
					taxid = lline[0]
					rank = lline[3]
					evalue = lline[4]

					if rank != 'species':
						res = TF.getSpeciesFromSubspecies(taxid)
						if res is None:
							continue
						taxid = str(res)

					if taxid not in entries:
						entries[taxid] = -1

					if evalue == '0.0':
						tid = 200
					elif 'e' in evalue:
						tid = -1 * int(evalue.split('e')[1])
					else:
						tid = math.ceil(math.log10(float(evalue)))

					if tid > entries[taxid]:
						entries[taxid] = tid

		with open('interactivetables/{}.tsv'.format(k), 'w') as outfile:
			for k in sorted(entries.keys()):
				outfile.write('{}\t{}\n'.format(k, entries[k]))


def uniqueNames():
	'''
	For each protein, creates lists with the species that possess this protein. Two lists are generated: One with unique Species names and one with unique taxonomy ids. In addition a general.names and general.taxids are generated with all Species that occur in the lists of all proteins.

	:creates: `names/*.names`
	:creates: `taxids/*.taxids`
	'''

	li = glob('combinedtables/*.tsv')	##### Use proteinlist.txt

	os.makedirs('names', exist_ok=True)
	os.makedirs('taxids', exist_ok=True)

	totalNames = set()
	totalTaxids = set()

	for fn in li:
		basename = os.path.basename(fn).replace('.tsv', '')
		with open(fn, 'r') as f:
			names = set()
			taxids = set()
			next(f)	# skip header
			for line in f:
				elems = line.split('\t')
				names.add(elems[2])
				taxids.add(int(elems[0]))

		names = sorted(list(names))
		taxids = sorted(list(taxids))

		with open('names/{}.names'.format(basename), 'w') as f:
			for name in names:
				f.write(name + '\n')

		with open('taxids/{}.taxids'.format(basename), 'w') as f:
			for taxid in taxids:
				f.write(str(taxid) + '\n')

		totalNames.update(names)
		totalTaxids.update(taxids)

	with open('names/general.names', 'w') as f:
		for name in sorted(list(totalNames)):
			f.write(name + '\n')

	with open('taxids/general.taxids', 'w') as f:
		for taxid in sorted(list(totalTaxids)):
			f.write(str(taxid) + '\n')


class NodeSanitizer():
	'''
	The NoteSanitizer is needed to create Newick files. It filters out bad characters and replaces them with allowed, similar characters. If there were unknown characters, they are printed out.
	'''

	def __init__(self):
		'''
		Initiates the class.
		'''

		self.badChars = set()
		self.goodChars = set(('A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z', ' ', 'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm', 'n', 'o', 'p', 'q', 'r', 's', 't', 'u', 'v', 'w', 'x', 'y', 'z', '0', '1', '2', '3', '4', '5', '6', '7', '8', '9', '^', '_', '=', '-', '/', '.', '*'))
		self.toReplace = {'(': '<', ')': '>', '[': '<', ']': '>', '#': '_', ':': '_', '+': '', "'": '', ',': '_'}


	def sanitize(self, nodes):
		'''
		Goes through a list of nodes and replace bad characters.

		:param nodes: The list of nodes (list of strings) to be sanitized
		:returns: A list with sanitized nodes (strings)
		'''

		result = []

		for node in nodes:
			newName = []
			for char in node:
				if char not in self.goodChars:
					if char in self.toReplace:
						newName.append(self.toReplace[char])
					else:
						self.badChars.add(char)
						newName.append('!')
				else:
					newName.append(char)
			result.append(''.join(newName))

		return result


	def printBadChars(self):
		'''
		If unknown characters were found in the last call(s) of `NodeSanitizer.sanitize()`, these characters are printed. Otherwise, nothing happens.
		'''

		if self.badChars:
			print('Unknown chars found:')
			for elem in sorted(list(self.badChars)):
				print(elem)


def makeNewick():
	'''
	Creates Newick trees for every taxonomy-id list in the taxids folder.

	:creates: `trees/*.tre`
	'''

	li = glob('taxids/*.taxids')	#### Read from proteinlist.txt? But general.taxid must also be taken.

	os.makedirs('trees', exist_ok=True)

	Sani = NodeSanitizer()

	for fname in li:
		outfn = 'trees/{}.tre'.format(os.path.basename(fname).replace('.taxids', ''))
		with open(fname, 'r') as infile:
			lineages = []

			for line in infile:
				lineages.append(TF.getLineage(int(line), display='both'))

			tree = {}

			for line in lineages:
				if line[0] not in tree:
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
					sanitizedNodes = Sani.sanitize(newnodes)
					newick = newick.replace(node, '(' + ','.join(sanitizedNodes) + ')' + node)
					nodes.extend(newnodes)

		with open(outfn, 'w') as outfile:
			outfile.write(newick)

	Sani.printBadChars()


def _getTreeElements(tree, returnSet = False, splitting = True):
	'''

	:param tree: The Newick tree as a string
	:param returnSet: Shall the tree elements be returned as a set (True) or a dictionary where each element points to an empty list (False)?
	:param splitting: If the Newick tree elements are in the form `Name^Taxid`, this must be True. If it is in the form `taxid`, this must be False.
	:returns: Either a list with taxonomy ids or a dictionary with taxonomy ids as keys and empty lists as value, depending on `returnSet`
	'''

	tree = tree.replace('(', '\n').replace(')', '\n').replace(',', '\n').replace(';', '')
	treelist = tree.split('\n')

	if returnSet:
		elements = set()
	else:
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

		if returnSet:
			elements.add(int(line))
		else:
			elements[int(line)] = []

	return elements


def _getCode(number):
	'''
	Returns a two-character code depending on the number. The form is: 1 -> aA, 2 -> aB, 26 -> aZ, 27 -> bA, ...

	:param number: The number to create the code from
	:returns: The code
	'''

	return chr((int(number/26) % 26) + 97) + chr((number % 26) + 65)


def treeAttributes():
	'''
	For each element in the general tree, determine the proteins (attributes) for this element. Each attribute will be a two-letter string. The meaning of each string will be written to `attributekeys.txt`.

	:uses: `trees/general.tre`
	:creates: `attributekeys.txt` (containing an explanation of the keys)
	:creates: `attributes.txt` (containing the keys for each protein)
	'''

	proteinlist = CR.getProteinNames()

	with open('trees/general.tre', 'r') as f:
		tree = f.read()

	allElements = _getTreeElements(tree, returnSet = False, splitting = True)

	with open('attributekeys.txt', 'w') as keysout:
		for i in range(len(proteinlist)):
			code = _getCode(i)
			keysout.write('{} -> {}\n'.format(code, proteinlist[i]))
			try:
				with open('trees/{}.tre'.format(proteinlist[i]), 'r') as f:
					specificElements = _getTreeElements(f.read(), returnSet = True, splitting = True)
					for k in specificElements:
						if k in allElements:
							allElements[k].append(code)
			except IOError:
				print('File {} not found.'.format(proteinlist[i]))

	with open('attributes.txt', 'w') as f:
		for k in sorted(list(allElements.keys())):
			f.write('{}\t{}\n'.format(k, ''.join(allElements[k])))


def makeHistograms():
	'''
	For each seed protein, the distribution of hits is shown in a histogram.

	:creates: `histograms/*.png`
	'''

	proteins = {}

	### use the ConfigReader
	with open('proteinlist.txt', 'r') as f:
		for line in f:
			line = line.split('#')[0].rstrip()
			protein, fasta = line.split()
			if protein not in proteins:
				proteins[protein] = []
			length = 0
			with open('fastas/'+fasta, 'r') as fastafile:
				next(fastafile)
				for line in fastafile:
					length += len(line.rstrip())
			proteins[protein].append(length)

	li = sorted(proteins.keys())

	os.makedirs('histograms', exist_ok=True)

	cm = plt.cm.get_cmap('rainbow')

	for protein in li:
		fn = 'combinedtables/{}.tsv'.format(protein)
		print('Histogram: {:<50}'.format(protein), end='\r')
		values = np.loadtxt(fn, dtype=np.int, comments=None, delimiter='\t', skiprows=1, usecols=(5,))

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

		seeds = proteins[protein]

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

		fig = plt.figure(1, figsize=(12,6))
		ax = fig.add_subplot(1,1,1)
		n, bins, patches = ax.hist(values, bins=bins)

		# The following block is to color the bars
		bin_centers = 0.5 * (bins[:-1] + bins[1:])
		col = bin_centers - min(bin_centers)
		col /= max(col)
		for c, p in zip(col, patches):
			plt.setp(p, 'facecolor', cm(c))

		ax.text(0.05, 0.95, text, transform=ax.transAxes, horizontalalignment='left', verticalalignment='top')

		ax.text(0.95, 0.95, sizes, transform=ax.transAxes, horizontalalignment='right', verticalalignment='top')

		fig.savefig('histograms/{}.png'.format(os.path.basename(fn).replace('.tsv', '')))

		fig.clear()


def showBlastMapping():
	'''
	For each protein, create an overview over where the hits where mapped over the length of the protein.

	:creates: `blastmappings/*.png`
	'''

	os.makedirs('blastmappings', exist_ok=True)

	fnames = sorted(list(CR.getProteinFiles(suffix = '.fasta')))

	fnt = ImageFont.load_default()

	for fname in fnames:
		print('Mapping {:<50}'.format(fname), end='\r')

		query_length = 0
		with open('fastas/' + fname, 'r') as f:
			next(f)
			for line in f:
				query_length += len(line.rstrip())

		counters = [np.zeros(query_length, np.int) for x in range(6)]
		numHsps = [0] * 6

		with open('blastresults/' + fname.replace('.fasta', '.xml'), 'r') as f:
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
						numHsps[n] += 1

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
			dr.text((2, 60 + 100 * n), 'n = {}'.format(numHsps[n]), (0, 0, 0), fnt)

		colors = [(0, 0, 0), (0, 0, 200), (0, 200, 0), (200, 0, 200), (200, 0, 0), (150, 150, 0)]

		for n in range(int(query_length / 100)):
			col = 160 + n*100
			dr.line([(col, 0), (col, 600)], fill=(125, 125, 125), width=1)

		for n in range(6):
			for col, thickness in enumerate(counters[n]):
				dr.line([(col + 60, n*100), (col + 60, thickness + n*100)], fill=colors[n], width=1)

		#im.show()
		im.save('blastmappings/' + fname.replace('.fasta', '.png'))


def similarityMatrix():
	'''
	Creates a similarity matrix of proteins. For each protein, the Blast hits are compared to each other protein. The lowest e-value of each protein-protein combination is saved. This gives an impression of how similar two given proteins are and how much the results overlap.

	:uses: combined parsed Blast results in `combinedtables`
	:creates: `matrix.csv`
	'''

	values = {}
	for fn in glob('combinedtables/*.combined.tsv'):	##### Maybe use proteinlist.txt
		name = fn.replace('combinedtables/', '').replace('.combined.tsv', '')
		values[name] = [set() for _ in range(151)]
		with open(fn) as f:
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

	names = sorted(values.keys())

	for i, name1 in enumerate(names):
		res.append([])
		for j, name2 in enumerate(names):
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

	with open('matrix.csv', 'w') as out:
		out.write('\t' + '\t'.join(names) + '\n')
		for i, line in enumerate(res):
			out.write(names[i] + '\t' + '\t'.join(map(str, line)) + '\n')


def heatmap():
	raise NotImplementedError('The heatmap is not implemented, yet.')


def intHeatmap():
	raise NotImplementedError('The interactive heatmap is not implemented, yet.')


tasknames = ['blast', 'parse', 'combine', 'comheat', 'unique', 'newick', 'attrib', 'hist', 'map', 'heatmap', 'intheat', 'matrix']

tasks = {'blast': ('run Blast', runBlast),
'parse': ('parse Blast results', parseBlastResults),
'combine': ('combine parsed results', combineParsedResults),
'comheat': ('combine parsed results for heatmaps', tablesForInteractiveHeatmap),
'unique': ('create unique lists of names and taxids', uniqueNames),
'newick': ('create Newick tree for each protein', makeNewick),
'attrib': ('determine tree attributes', treeAttributes),
'hist': ('create histograms with Blast hits for each protein', makeHistograms),
'map': ('create hit mapping diagrams for each protein', showBlastMapping),
'heatmap': ('create a heatmap (image)', heatmap),
'intheat': ('create an interactive heatmap (website)', dynHeatmap),
'matrix': ('create a similarity matrix of all proteins', similarityMatrix)}


def runWorkflow(start, end=''):
	'''
	Starts the workflow from `start` until `end`. If `end` is empty, the workflow is run until the last element of the workflow.

	:param start: The name of the first element to run.
	:param end: The name of the last element to run or a falsy value if the workflow shall be run until the end.
	'''

	if start not in tasknames:
		raise ValueError('{} is no valid task.'.format(start))

	if not end:
		endidx = len(tasknames)
	elif end not in tasknames:
		raise ValueError('{} is no valid task.'.format(end))
	else:
		endidx = tasknames.index(end) + 1

	startidx = tasknames.index(start)
	for taskname in tasknames[startidx:endidx]:
		print('{}: "{}"'.format(taskname, tasks[taskname][0]))
		task = tasks[taskname][1]
		task()


if __name__ == '__main__':
	import argparse
	import textwrap

	workflow = '\n'.join(textwrap.wrap("""The following is a list of the workflow. The names or numbers can be used for the -s or -o arguments.""", width = 80))

	workflow += '\n\n' + '\n'.join(('{:>2}. {:<8} {}'.format(i, name, tasks[name][0]) for i, name in enumerate(tasknames)))

	parser = argparse.ArgumentParser(description='This module provides you with tools to run phylogenetic analyses. Exactly one argument must be given.')

	parser.add_argument('-l', '--list', action='store_true', help='Shows the whole workflow with information and exits')
	parser.add_argument('-a', '--all', action='store_true', help='Run the full workflow without Blast')
	parser.add_argument('-b', '--blast', action='store_true', help='Run the full workflow including Blast')
	parser.add_argument('-s', '--startfrom', default='', help='Run from and including this step [e.g. 7 or hist]')
	parser.add_argument('-o', '--only', default='', help='Run only the given step [e.g. 4 or unique]')

	args = parser.parse_args()

	if args.list:
		parser.print_help()
		print('')
		print(workflow)
		sys.exit()

	numArguments = args.all + args.blast + bool(args.startfrom) + bool(args.only)

	if numArguments != 1:
		parser.print_help()
		sys.exit()

	if args.all:
		runWorkflow('parse')
	elif args.blast:
		runWorkflow('blast')
	elif args.startfrom:
		try:
			a = int(args.startfrom)
		except ValueError:
			if args.startfrom in tasknames:
				runWorkflow(args.startfrom)
			else:
				parser.print_help()
		else:
			if a < len(tasknames):
				runWorkflow(tasknames[a])
			else:
				parser.print_help()
	elif args.only:
		try:
			a = int(args.only)
		except ValueError:
			if args.only in tasknames:
				task = tasks[args.only][1]
				task()
			else:
				parser.print_help()
		else:
			if a < len(tasknames):
				task = tasks[tasknames[a]][1]
				task()
			else:
				parser.print_help()
	else:
		print('This should not happen!')
		parser.print_help()







