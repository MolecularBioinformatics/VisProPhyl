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


class General():

	def __init__(self):
		self.proteinSet = None
		self.limitsDict = None


	def readLimits(self):
		'''
		Reads limits.txt and returns a dict in the format: {'protein': (evalue, length), ...} where evalue is the negative exponent of the evalue (int) and length is the length limit for the protein (int).
		'''

		if self.limitsDict is not None:
			return self.limitsDict

		with open('limits.txt', 'r') as f:
			limits = {}
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
				limits[lline[0]] = (evalue, length)

		self.limitsDict = limits

		return limits


	def readProteinsWithFN(self, prefix = '', suffix = ''):
		'''
		Reads proteinlist.txt and returns a dict in the format: {'protein': set('file1', 'file2', ...), ...}. The files are just the basenames without suffix. A prefix (a path) can be added as well as a suffix. readProteins('abc/', '.txt') will turn 'file1' into 'abc/file1.txt'.
		'''

		with open('proteinlist.txt', 'r') as f:
			proteins = {}
			for line in f:
				line = line.split('#')[0].strip()
				if not line:
					continue
				elems = line.split()
				if elems[0] not in proteinlist:
					proteins[elems[0]] = set()
				proteins[elems[0]].add(prefix + elems[1].replace('.fasta', suffix))

		return proteins


	def readProteins(self):
		'''
		Reads proteinlist.txt and returns a set in the format: {'protein1', 'protein2', ...}.
		'''

		if self.proteinSet is not None:
			return self.proteinSet

		self.proteinSet = set(self.readProteinsWithFN().keys())

		return self.proteinSet


	def readProteinFNs(self, prefix = '', suffix = ''):
		'''
		Reads proteinlist.txt and returns a set in the format: {'proteinfile1', 'proteinfile2', ...}. The files are just the basenames without suffix. A prefix (a path) can be added as well as a suffix. readProteins('abc/', '.txt') will turn 'proteinfile1' into 'abc/proteinfile1.txt'.
		'''

		p = self.readProteinsWithFN(prefix, suffix)
		fns = set()
		for k in p:
			fns.update(p[k])

		return fns


TF = TaxFinder()

g = General()


"""
def readLimits():
	return g.readLimits()

def readProteins():
	return g.readProteins()

def readProteinsWithFN(prefix = '', suffix = ''):
	return g.readProteinsWithFN(prefix, suffix)

def readProteinFNs(prefix = '', suffix = ''):
	return g.readProteinFNs(prefix, suffix)
"""


def _runBlast(query, db = '/home/mathias/projects/nr/nr', evalue = 0.001, maxthreads = cpu_count()):
	cmd = NcbiblastpCommandline(query = query, db = db, evalue = evalue, outfmt = 5, out = query.replace('.fasta', '.xml').replace('fastas/', 'blastresults/'), num_threads = maxthreads, max_target_seqs = 20000)
	stdout, stderr = cmd()
	print('Blasting', fname, 'done')
	if stdout:
		print(stdout)
	if stderr:
		print(stderr, file=sys.stderr)
	else:
		print('No errors')


def runBlast():
	li = glob('fastas/*.fasta')		### maybe get from proteinlist.txt
	os.makedirs('blastresults', exist_ok=True)
	for i, fname in enumerate(li):
		outstring = 'Now doing {:' + str(len(str(len(li)))) + 'd}/{}: {:<30}'
		print(outstring.format(i+1, len(li), fname), end='\r')
		_runBlast(fname)



def _getTopFromBlast(blastXML, TF, top = 0, exContaminSpecies=True, outfile=None, newHeader=True):

	contaminantSpecies = set((118797, 59538, 7213))	# Lipotes vexillifer, Pantholops hodgsonii, Ceratitis capitata

	if outfile is None:
		outfile = 'resulttables/{}.tsv'.format(os.path.basename(blastXML))

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


def parseBlastResults():
	li = glob('blastresults/*.xml')	### This may be taken from the proteinlist.txt

	#li = ['blastresults/Athaliana_nadA_QS.xml']

	os.makedirs('resulttables', exist_ok=True)

	outstring = 'Now doing {:' + str(len(str(len(li)))) + 'd}/{}: {:<30}'

	for i, fname in enumerate(li):
		basename = os.path.basename(fname)
		print(outstring.format(i+1, len(li), basename), end='\r')
		getTopFromBlast(fname, TF = TF, top = 0, outfile='resulttables/{}.tsv'.format(basename), exContaminSpecies=True)

	print('')


def combineTSV():
	os.makedirs('combinedtables', exist_ok=True)

	limits = g.readLimits()
	proteins = g.readProteinsWithFN(prefix = 'resulttables/', suffix = '.xml.tsv')

	for k in sorted(proteins):
		print('combining tsv files... {:<40}'.format(k), end='\r')
		outfn = 'combinedtables/{}.combined.tsv'.format(k)
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


def tablesForDynamicHeatmap():
	proteins = g.readProteinsWithFN(prefix = 'resulttables/', suffix = '.xml.tsv')

	os.makedirs('dynamictables', exist_ok=True)

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

		with open('dynamictables/{}.dyn.tsv'.format(k), 'w') as outfile:
			for k in sorted(entries.keys()):
				outfile.write('{}\t{}\n'.format(k, entries[k]))


def uniqueNames():
	li = glob('combinedtables/*.combined.tsv')

	os.makedirs('names', exist_ok=True)
	os.makedirs('taxids', exist_ok=True)

	totalNames = set()
	totalTaxids = set()

	for fn in li:
		basename = os.path.basename(fn).replace('.combined.tsv', '')
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

	def __init__(self):
		self.badChars = set()
		self.goodChars = set(('A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z', ' ', 'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm', 'n', 'o', 'p', 'q', 'r', 's', 't', 'u', 'v', 'w', 'x', 'y', 'z', '0', '1', '2', '3', '4', '5', '6', '7', '8', '9', '^', '_', '=', '-', '/', '.', '*'))
		self.toReplace = {'(': '<', ')': '>', '[': '<', ']': '>', '#': '_', ':': '_', '+': '', "'": '', ',': '_'}


	def sanitize(self, nodes):
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
		if self.badChars:
			print('Unknown chars found:')
			for elem in sorted(list(self.badChars)):
				print(elem)


def makeNewick():
	li = glob('taxids/*.taxids')	#### Read from proteinlist.txt?

	os.makedirs('trees', exist_ok=True)

	TF = TaxFinder()

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


def _getTreeElements(tree, giveSet = False, splitting = True):
	tree = tree.replace('(', '\n').replace(')', '\n').replace(',', '\n').replace(';', '')
	treelist = tree.split('\n')

	if giveSet:
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

		if giveSet:
			elements.add(int(line))
		else:
			elements[int(line)] = ''

	return elements


def _getCode(number, easy = False):
	if easy:
		return chr((i % 26) + 65)
	else:
		return chr((int(i/26) % 26) + 97) + chr((i % 26) + 65)


def treeAttributes():
	proteinlist = g.readProteins()

	with open('trees/general.tre', 'r') as f:
		tree = f.read()

	allElements = _getTreeElements(tree, giveSet = False, splitting = True)

	with open('attributekeys.txt', 'w') as keysout:
		for i in range(len(proteinlist)):
			code = _getCode(i)
			keysout.write('{} -> {}\n'.format(code, proteinlist[i]))
			try:
				with open('trees/{}.tre'.format(proteinlist[i]), 'r') as f:
					specificElements = _getTreeElements(f.read(), giveSet = True, splitting = True)
					for k in specificElements:
						if k in allElements:
							allElements[k] += code
			except IOError:
				print('File {} not found.'.format(proteinlist[i]))

	with open('attributes.txt', 'w') as f:
		for k in sorted(list(allElements.keys())):
			f.write(str(k) + '\t' + allElements[k] + '\n')


def makeHistograms():
	proteins = {}

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
		fn = 'combinedtables/{}.combined.tsv'.format(protein)
		print(protein.ljust(50), end='\r')
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

		fig.savefig('histograms/{}.png'.format(os.path.basename(fn).replace('.combined.tsv', '')))

		fig.clear()


def showBlastMapping():
	os.makedirs('blastmappings', exist_ok=True)

	fnames = sorted(list(g.readProteinFNs(suffix = '.fasta')))

	fnt = ImageFont.load_default()

	for fname in fnames:
		print(fname.ljust(50), end='\r')

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
	values = {}
	for fn in glob('combinedtables/*.combined.tsv'):
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


def dynHeatmap():
	raise NotImplementedError('The interactive heatmap is not implemented, yet.')


tasknames = ['blast', 'parse', 'combine', 'comheat', 'unique', 'newick', 'attrib', 'hist', 'map', 'heatmap', 'dynheat', 'matrix']

tasks = {'blast': ('run Blast', runBlast),
'parse': ('parse Blast results', parseBlastResults),
'combine': ('combine parsed results', combineTSV),
'comheat': ('combine parsed results for heatmaps', tablesForDynamicHeatmap),
'unique': ('create unique lists of names and taxids', uniqueNames),
'newick': ('create Newick tree for each protein', makeNewick),
'attrib': ('determine tree attributes', treeAttributes),
'hist': ('create histograms with Blast hits for each protein', makeHistograms),
'map': ('create hit mapping diagrams for each protein', showBlastMapping),
'heatmap': ('create a heatmap (image)', heatmap),
'dynheat': ('create an interactive heatmap (website)', dynHeatmap),
'matrix': ('create a similarity matrix of all proteins', similarityMatrix)}


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







