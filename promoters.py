#!/usr/bin/env python3
import os
import sys
from ete3 import Tree, TreeStyle, faces, AttrFace
import pandas as pd
import numpy as np
from scipy.stats.mstats import gmean
from time import strftime
import subprocess
from Bio import Entrez, SeqIO, motifs
from Bio.Alphabet import DNAAlphabet, IUPAC
from copy import deepcopy
import re
from glob import glob
from PIL import Image
from taxfinder import TaxFinder
from phylogenetics import NodeSanitizer
from phylotree import Clusters

VERBOSITY = False
Entrez.email = 'nicolai.von-kuegelgen@student.uni-tuebingen.de'
Entrez.tool = 'PromoterAnalysis-NvK'

# ToDo: extras for config file
#	- functions (& threshhold) used in treeClusters
#	- numbers & normalisation used in _similarity
#	- numbers in motif analysis function


# ToDo: clean state & file checking
# - to avoid wierd bugs always clean folders into which files will be downloaded
#  possible exceptions: featuretables (should be cleaned/updated regularly but not with every run)
# - see that all functions give (nice) errors if necessary files are missing

def _myGetBasename(name):
	'''
	Strips the path and the extension from a filename. E.g. /a/path/to/file.png -> file

	:param name: The file and path to be processed
	:returns: The processed filename
	'''

	return os.path.splitext(os.path.basename(name))[0]


# ToDo: forced call of this should update all files from scratch / 'reinitialise'
def init():
	'''
	Creates some files to start with a new project. Existing files will not be overwritten.

	:creates: infer_NJ_nucleotide_pair.mao
	:creates: assembly_summary_refseq.txt
	'''
	tx = taxids()

	origDir = os.path.dirname(os.path.realpath(__file__))

	if not os.path.isfile('infer_NJ_nucleotide_pair.mao'):
		with open('infer_NJ_nucleotide_pair.mao', 'w') as out, open(os.path.join(origDir, 'infer_NJ_nucleotide_pair.mao'), 'r') as f:
			out.write(f.read())

	if not os.path.isfile('promotersConfig.txt'):
		with open('promotersConfig.txt', 'w') as out, open(os.path.join(origDir, 'promotersConfig.txt'), 'r') as f:
			out.write(f.read())

		if not os.path.isfile('tree_tp_prune.txt'):
			with open('tree_tp_prune.txt', 'w') as out, open(os.path.join(origDir, 'tree_tp_prune.txt'), 'r') as f:
				out.write(f.read())


# opt: also use this class to provide/make a base tree for the given g-quality
class taxids:
	def __init__(self, qlevel=1):
		"""Provide an instance with all taxids of genomes in refseq for a given quality level, their accession numbers and ftp-download locations.
		Neccessary files will be downloaded/generated automatically.
		qlevel 0: all gemones, qlevel 1: assembly of at least one chromosome, qlevel2: representative/reference genome
		"""
		self.ids = set()
		self.map = dict()
		#not really needed, but maybe someday
		self.qlevel = qlevel

		if not os.path.isfile('refseq_taxids_l{}.txt'.format(qlevel)):
			self.generatefiles(qlevel)

		with open('refseq_taxids_l{}.txt'.format(qlevel), 'r') as f:
			next(f)
			for line in f:
				tax, acc, ftp = line.rstrip().split('\t')
				self.ids.add(tax)
				self.map[tax] = (acc, ftp)


	def generatefiles(self, qlevel):
		"""Creates necessary files for the class to work."""

		if VERBOSITY:
			print("Initialising refseq_genome taxid files.")

		if not os.path.isfile('assembly_summary_refseq.txt'):
			subprocess.call('rsync -q --copy-links --times rsync://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_refseq.txt .'.split())

		taxid_cache = {} #tax: (acc, quality, category)

		with open('assembly_summary_refseq.txt', 'r') as f, open('refseq_taxids_l{}.txt'.format(qlevel), 'w') as out:
			next(f)
			next(f)
			out.write('#Taxid\tGenome-Accseion\tftp-folder-location\n')
			for line in f:
				lline = line.rstrip().split('\t')

				completeness = ('Contig', 'Scaffold', 'Chromosome', 'Complete Genome')
				category = ('na', 'representative genome', 'reference genome')

				if qlevel == 2 and lline[4] == 'na':  # refseq categeory unspecified (= not a representative/reference genome)
					continue
				if qlevel and lline[11] in ('Scaffold', 'Contig'):  # completeness of genome at least on one chromosome
					continue

				# Only write a second entry for a given species if the second one is 'better' then the first
				# Thereby the second will overwrite the first one when reading the file
				if lline[5] in taxid_cache:
					if completeness.index(lline[11]) > completeness.index(taxid_cache[lline[5]][1]):
						pass
					elif category.index(lline[4]) > category.index(taxid_cache[lline[5]][2]):
						pass
					else:
						continue

				taxid_cache[lline[5]] = [lline[0], lline[11], lline[4]]

				out.write("\t".join((lline[5], lline[0], lline[19])) + '\n')  # taxid, acc, ftp-link


def readConfig(defaultargs):
	'''
	Read the config file and run the workflow accordingly.

	:uses: promotersConfig.txt
	:param defaultargs: argparser instance with default flags
	'''

	if not os.path.isfile('promotersConfig.txt'):
		print('Could not find config file, please run init (-i) and edit the config base file.')
		sys.exit()

	def checkbool(x):
		if x.lower() == 'true':
			return True
		elif x.lower() == 'false':
			return False
		else:
			raise ValueError

	vartransform = {'verbose': checkbool, 'startfrom': lambda x: x, 'only': lambda x: x, 'end': lambda x: x, 'path': lambda x: x,
					'evalue': int, 'genomequality': int, 'phylopath': lambda x: x, 'length': int,
					'drop': lambda x: x.split(), 'threshold': float, 'skipclusters': lambda x: x.split()}

	runargs = []
	current = None

	# Read Config file, flags before a block (current == None) are set as default
	# Specific exceptions should catch the most common problems (not excessively tested)
	with open('promotersConfig.txt', 'r') as f:
		try:
			for line in f:
				line = line.split('#')[0].rstrip()

				if not line:
					continue

				if line.startswith('!!'):
					setattr(currentargs, 'protein', current)
					runargs.append(currentargs)
					continue

				elif line.startswith('!'):
					current = line[1:]
					currentargs = deepcopy(defaultargs)
					continue

				key, value = line.split(':')
				value = vartransform[key](value)

				if current is None:
					setattr(defaultargs, key, value)
				else:
					setattr(currentargs, key, value)

		except ValueError:
			print('Error in Config file: Error with flag values (empty value, too many `:`, wrong bool value)')
			sys.exit()

		except KeyError:
			print('Error in Config file: Unknown flag given')
			sys.exit()

		except NameError:
			print('Error in Config file: Block Syntax is wrong (!Protein .... !!Protein)')
			sys.exit()

		except:
			print('Unknown Error in Config file')
			sys.exit()


	for args in runargs:
		global PATH, VERBOSITY
		PATH = args.path
		VERBOSITY = args.verbose

		if VERBOSITY:
			print('Running workflow for: {}'.format(args.protein))

		if PATH:
			os.makedirs(PATH, exist_ok=True)

		if args.only:
			runWorkflow(args, args.only, args.only)
		else:
			runWorkflow(args, args.startfrom, args.end)


#  If this script is to be used without the phylogenetics.py workflow (up to step 2) then this function will have to be
#  substituted with/adapted to a parser for a list of protein Accession-IDs and probably also the TaxFinder module/class
def _getUniqueAccs(protein, taxid, cutoff=50, phylopath=''):
	"""Read combinedtables from Phylogenetics workflow.
	 Extract Accesion numbers and taxids for hits, apply refseq as filter, save to file

	 :param protein: a protein name from the phylogenetics workflow
	 :param taxid: an instance of the taxids class
	 :param cutoff: worst evalue (e^-{cutoff} accepted for hits
	 :creates: {prortein}_hits.txt
	 """

	# add option to read from config file
	# search for phylogenetics workflow
	if os.path.isfile("combinedtables/{}.tsv".format(protein)):
		p = "combinedtables/{}.tsv".format(protein)
	elif os.path.isfile(os.path.join(PATH, "combinedtables/{}.tsv".format(protein))):
		p = os.path.join(PATH, "combinedtables/{}.tsv".format(protein))
	else:
		if not phylopath and VERBOSITY:
			print('Could not locate combined tables from phylogenetics workflow.')
			phylopath = input('Enter path to file(s):')

		if os.path.isfile(phylopath):
			func = os.path.dirname()
		else:
			func = lambda x: x

		if os.path.isfile(os.path.join(func(phylopath), "{}.tsv".format(protein))):
			p = os.path.join(func(phylopath), "{}.tsv".format(protein))
		elif os.path.isfile(os.path.join(func(phylopath), "combinedtables/{}.tsv".format(protein))):
			p = os.path.join(func(phylopath), "combinedtables/{}.tsv".format(protein))
		else:
			print('Path to phylogenetics workflow seems incorrect.\nExiting script.')
			sys.exit()

	if VERBOSITY:
		print('Getting usable blast hits for {:<10}'.format(protein), end='\r')

	with open(p, "r") as f:
		next(f)
		accs = set()
		mapping = dict()
		for line in f:
			lline = line.split('\t')
			tax = lline[0]
			acc = lline[1]
			evalue = lline[4]

			#Only work with hits, that have a (good) genome in the refseq db
			if tax not in taxid.ids:
				continue

			#Only accept refseq protein hits, since others will not be found in the refseq genomes anyway (probably? not 100% tested)
			if acc[2] != '_':
				continue

			#reduce hits by evalue
			if evalue == '0.0':
				evalue = 201
			elif 'e-' in evalue:
				evalue = int(evalue.split('e-')[1])
			else:
				evalue = 0

			if evalue < cutoff:
				continue

			accs.add(acc)
			mapping[acc] = tax

	with open(os.path.join(PATH, '{}_hits.tsv'.format(protein)), 'w') as out:
		for acc in accs:
			out.write(acc+'\t'+mapping[acc]+'\n')


def _loadGenomeFeatures(protein, taxid):
	"""Download genome annotation tables for all hits of given protein from NCBI ftp-server

	:uses: Promoters/{protein}_hits.tsv
	:param protein: a protein name from the phylogenetics workflow
	:param taxid: an instance of the taxids class
	:creates: featuretables/*.txt
	"""

	if VERBOSITY:
		print('Downloading Genome annotation tables for hits of {:<10}'.format(protein), end='\r')

	accs = set()
	os.makedirs('featuretables', exist_ok=True)

	with open(os.path.join(PATH, '{}_hits.tsv').format(protein), 'r') as f:
		for line in f:
			acc, tax = line.rstrip().split('\t')
			accs.add((acc, tax))

	cmd = 'rsync --copy-links -q --times rsync{}_feature_table.txt.gz featuretables/'

	for acc, tax in accs:
		if os.path.isfile('featuretables/{}_feature_table.txt'.format(tax)):
			continue

		if VERBOSITY:
			print('Downloading Genome annotation table for species: {:<10}'.format(tax), end='\r')

		ftp = taxid.map[tax][1]
		fn = ftp.split('/')[-1]

		ret = subprocess.call(cmd.format(ftp[3:]+'/'+fn).split()) #downlaod from ncbi
		p = 'featuretables/{}_feature_table.txt.gz'.format(fn)
		ret += subprocess.call('gzip -d {}'.format(p).split()) #unpack
		p = 'featuretables/{}_feature_table.txt'.format(fn)
		p2 = 'featuretables/{}_feature_table.txt'.format(tax)
		ret += subprocess.call('mv {} {}'.format(p, p2).split()) # rename
		if ret:
			print(ret) #should print errors


# Searching through genome annotation files may be faster with reg-ex (probably only when applying regex to whole file)
def _extractPromoterPos(protein):
	"""Extract start & strand of Promoter for all (usable) hits a of protein
	For each write Prot-Acc, TaxID, Genome-Acc, StartPos, Strand, GeneSym, Prot-Name into File

	:uses: hits_{protein}.tsv, featuretables/*.txt
	:param protein: a protein name from the phylogenetics workflow
	:creates: {protein}_hitfeatures.tsv
	"""

	#sort hits by species, so feature tables have to be accessed only once per species
	tax_accs = dict()
	with open(os.path.join(PATH, '{}_hits.tsv'.format(protein)), 'r') as f:
		for line in f:
			acc, tax = line.rstrip().split('\t')
			if tax in tax_accs:
				tax_accs[tax] += [acc]
			else:
				tax_accs[tax] = [acc]

	with open(os.path.join(PATH, '{}_hitfeatures.tsv'.format(protein)), 'w') as out:
		out.write('#Prot-Acc\tTaxID\tGenome-Acc\tStartPos\tEndPos\tStrand\tGeneSym\tProt-Name\n')
		for tax, accs in tax_accs.items():
			if VERBOSITY:
				print('Extracting promoter sequences for hits of {}, species: {:<10}'.format(protein, tax), end='\r')

			with open('featuretables/{}_feature_table.txt'.format(tax), 'r') as f:
				for line in f:
					#safe some time / only look for genes
					if not line.startswith('CDS'):
						continue

					#without splitting first might also find 'related_accession', maybe not though after sorting for CDS
					lline = line.rstrip().split('\t')
	#Columns in *_feature_table.txt
	#feature0	class1	assembly2	assembly_unit3	seq_type4	chromosome5	genomic_accession6	start7	end8	strand9
	#product_accession10	non-redundant_refseq11	related_accession12	name13	symbol14	GeneID15	locus_tag16 feature_interval_length17	product_length18	attributes19

					if lline[10] not in accs:
						continue

					accs.remove(lline[10])

					outl = [lline[10], tax, lline[6], lline[7], lline[8], lline[9], lline[14], lline[13]]
					out.write('\t'.join(outl)+'\n')

					if len(accs) ==0:
						break
				else:
					#only with verbosity ?
					print('Extracting promoter sequences for hits of {}, species: {}, features not found: {}'.format(protein, tax, accs))


def _loadPromoterSeq(protein, length=1500):
	"""Using the information of {protein}_hitfeatures.tsv download the associated promoter sequence from NCBI

	:uses: {protein}_hitfeatures.tsv
	:param protein: a protein name from the phylogenetics workflow
	:param length: The number of NT in front of CDS start to pull
	:creates: PromoterSeqs_{Protein}/*.fasta, PromoterSeqs_{Protein}/annotation.tsv
	"""

	entries = set()
	os.makedirs(os.path.join(PATH, 'PromoterSeqs_{}'.format(protein)), exist_ok=True)

	with open(os.path.join(PATH, '{}_hitfeatures.tsv'.format(protein)), 'r') as f, open(os.path.join(PATH, 'PromoterSeqs_{}/annotation.tsv'.format(protein)), 'w') as annotation:
		next(f)
		annotation.write('#Protein-Acc\tTaxID\tGenomic-Acc\tStart-Pos\tEnd-Pos\tstrand\tGene symbol\tGene name\n')
		for i, line in enumerate(f):
			lline = line.rstrip().split('\t')
			#lline: protein, taxid, genomicacc, start, end, strand, symbol, name

			tax, genome_acc, start, end, strand = lline[1:6]

			# Multiple proteins may come from one (unique) genomic locus
			# end of entry (e.g. from multiple isoforms) is irrelevant since we only look at start
			entry = (tax, genome_acc, start)
			if entry in entries:
				continue
			else:
				entries.add(entry)

			#Dont load same sequence again - mostly for test use
			#if os.path.isfile('PromoterSeqs_{}/{}.fasta'.format(protein, lline[6]+'_'+tax)): # / lline[0]
				#continue
			if VERBOSITY:
				print('Downloading promoter sequences for hits of {}. Found: {} ({} multiple){:<6}'.format(protein, len(entries), i-len(entries)+1, ''), end='\r')

			#start from the *_hitfeatures file is the CDS start, we want region downstream of that
			if strand == '+':
				end = int(start)
				start = int(start) - length
				strand = 1
			#If the gene is on the reverse strand, the Promoter is upstream of the gene
			else:
				start = int(end)
				end = int(end) + length
				strand = 2


			handle = Entrez.efetch(db="nucleotide",
								   id=genome_acc,
								   rettype="fasta",
								   strand=strand,
								   seq_start=start,
								   seq_stop=end)
			record = SeqIO.read(handle, "fasta")

			#Only short header in fasta file
			record.id = '{}^{}'.format(lline[0], tax)
			record.description = ''

			with open(os.path.join(PATH, 'PromoterSeqs_{}/{}.fasta'.format(protein, lline[0])), 'w') as out:
				SeqIO.write(record, out, 'fasta')

			#Complete Information in annotation file
			annotation.write('\t'.join(lline)+'\n')


def getPromoterSeqs(protein, evalue, length, genomequality, phylopath):
	'''wrapper for all functions needed to download Promoter Seqs

	:param protein: protein name from phylogenetics workflow
	:param evalue: maximum evalue for collection of hits from workflow
	:param length: number of bases before CDS start to pull as promoter sequence
	:creates: {prortein}_hits.txt, featuretables/*.txt, {protein}_hitfeatures.tsv, PromoterSeqs_{Protein}/*.fasta
	'''

	# determine/flag if protein is from phylogenetics workflow or somewhere else

	taxid = taxids(genomequality)
	#if phylogenetics
	_getUniqueAccs(protein, taxid, evalue, phylopath)
	#else: new parser for tsv or similar
	_loadGenomeFeatures(protein, taxid)
	_extractPromoterPos(protein)
	_loadPromoterSeq(protein, length)


def runClustal(protein):
	'''
	Wrapper to call run clustalo over all downloaded fasta files

	:uses: PromoterSeqs_{proteins}/*.fasta
	:param protein: protein name from phylogenetics workflow
	:creates: {protein}_promoters_aligned.fasta
	'''

	# TODO: check there are some fasta files there, Popen.wait() will otherwise deadlock

	# maybe only with 2nd level verbosity ?
	v = ''
	if VERBOSITY:
		v = '-v '

	pin = os.path.join(PATH, 'PromoterSeqs_{}/*.fasta'.format(protein))
	pout = os.path.join(PATH, '{}_promoters_aligned.fasta'.format(protein))

	files = glob(pin)

	cmd = 'clustalo -i - {}-o {} --force'.format(v, pout)

	# shell=True (& strcmd instead of list) required for the pipe to work with subprocess.call
	# It is however not recommended and probably bad, a different solution would be better
	# call(cmd, shell=True)

	p1 = subprocess.Popen(['cat']+files, stdout=subprocess.PIPE)
	p2 = subprocess.Popen(cmd.split(), stdin=p1.stdout)
	p2.wait()
	p1.stdout.close()


def _pw_similarity(seqs, names):
	'''
	Calculates a pairwise similarity matrix (using _similarity)
	:param seqs: list/tuple/... of sequences
	:param names: list/tuple/... of names of sequences (same lenght as seqs)
	:return: pd.Dataframe with pw-similarity matrix ('-' on diagonal)
	'''
	mat = pd.DataFrame(None, index=range(len(seqs)), columns=range(len(seqs)))
	for i, seq1 in enumerate(seqs):
		for j, seq2 in enumerate(seqs):
			if j > i:
				continue

			if i == j:
				mat[i][j] = '-'
			else:
				d = _similarity(seq1, seq2)
				mat[i][j] = d
				mat[j][i] = d
	mat.index = names
	mat.columns = names
	return mat


# The specific numbers here will have an impact on downstream results
# In the context of an MSA the normalisation of the score to len(seq) is only good to keep the number between 0 & 1
# ACTIVE: Normalisation to overlapping lenght (w/o trailing/leading gaps)
def _similarity(seq1, seq2):
	'''
	Calculate similarity score for seq1 & seq2 (normalised to 0-1). Sequences should be aligned (otherwise stops at end of shorter seq).
	Equal bases: score +1, unequal bases: score +0.1, gaps: score +0
	Normalisation: overlapping length (w/o trailing/leading gaps)

	:param seq1: str
	:param seq2: str
	:return: float or np.NaN if score is 0
	'''
	alphabet = ['A', 'C', 'T', 'G']
	score = 0

	for s1, s2 in zip(seq1, seq2):
		if s1 not in alphabet or s2 not in alphabet:
			continue

		if s1 == s2:
			score += 1.0
		else:
			score += 0.1

	if score:
		# raw score (normalise to full lenght, identical for all seqs because of gaps, length w/o gaps is also identical for all)
		# return score/len(seq1)

		# Normalise to overlapping lenght (w/ gaps): diff(min(start(seq1, seq2)), max(end(seq1, seq2)))
		first = lambda x: min(x.find(i) for i in alphabet if x.find(i)!=-1)
		last = lambda x: len(x) - first(x[::-1])
		return score / (max(last(seq1), last(seq2)) - min(first(seq1), first(seq2)))


	else:
		return np.NaN


def autosort(protein, drop):
	'''
	Takes a fasta file with multiple aligned sequences, calculates pw-similarity matrix & removes entries that have no similarity score with other entries (greedy alg.)

	:uses: {protein}_promoters_aligned.fasta
	:param protein: protein name from phylogenetics workflow
	:param drop: list (or similar) of entries that are to be dropped regardless of score
	:creates: {protein}_promoters_matrix.tsv, {protein}_promoters_kept.fasta, {protein}_promoters_dropped.fasta
	'''
	seqs = []
	names = []
	fname = os.path.join(PATH, protein + '_promoters_aligned.fasta')


	# this may be faster/easier with using SeqIO as parser
	with open(fname, 'r') as f:
		current = ''
		for line in f:
			line = line.rstrip()
			if line.startswith('>'):
				#first entry is empty
				if current:
					seqs.append(current)
				# #OLD/current fastas: >promoter|{acc}
				# names.append(line.split('|')[1])
				# NEW: >{acc}^{tax}
				names.append(line[1:])
				current = ''
				continue

			if line:
				current += line
		#last seqence
		seqs.append(current)

	if VERBOSITY:
		print('Calculating pairwise similarity matrix for {:<10}.'.format(protein))

	matrix = _pw_similarity(seqs, names)

	with open(os.path.join(PATH, protein+'_promoters_matrix.tsv'), 'w') as out:
		out.write('# Pairwise similarity scores for all sequences\n')
		matrix.to_csv(path_or_buf=out, sep='\t', na_rep='N/A')

	dropped = []
	# (true/false NaN for each cell).(count per col).(count total)
	while matrix.isnull().sum().sum():
		nans = matrix.isnull().sum()
		#Dropping all columns with the highest number of NaNs + the respective rows
		worst = nans[nans == max(nans)].index
		matrix.drop(worst, axis=0, inplace=True)
		matrix.drop(worst, axis=1, inplace=True)
		dropped += list(worst)

	if VERBOSITY:
		print('Following IDs were dropped to remove N/A from the similarity matrix:')
		print(', '.join(dropped))

	keep = list(matrix.index)

	i = 0

	with open(fname, 'r') as f, open(os.path.join(PATH, protein+'_promoters_kept.fasta'), 'w') as out, open(os.path.join(PATH, protein+'_promoters_dropped.fasta'), 'w') as out2:
		for line in f:
			if line.startswith('>'):
				# # OLD/current fastas: >promoter|{acc}
				# name = line.rstrip().split('|')[1]

				# NEW: >{acc}^{tax}
				name = line[1:].rstrip()
				if name in drop or name.split('^')[0] in drop or str(i) in drop:
					keeping = False
				else:
					keeping = name in keep
				i += 1

			if keeping:
				out.write(line)
			else:
				out2.write(line)


def runMegacc(protein):
	'''
	Wrapper to call run Megacc over sorted MSA

	:uses: {protein}_promoters_kept.fasta
	:param protein: protein name from phylogenetics workflow
	:creates: {protein}_promoters_megacc.nwk
	'''

	if not os.path.isfile('infer_NJ_nucleotide_pair.mao'):
		print('MegaCC config file is missing, please rerun initialisation.\nExiting script.')
		sys.exit()

	# only with 2nd level verbosity ??
	v = '-s '
	if VERBOSITY:
		v = ''

	pin = os.path.join(PATH, '{}_promoters_kept.fasta'.format(protein))
	pout = os.path.join(PATH, '{}_promoters_megacc.txt'.format(protein))
	treefn = pout.replace('.txt', '.nwk')

	# Megacc does NOT overwrite its own results, do it manually
	if os.path.isfile(treefn):
		os.remove(treefn)

	cmd = 'megacc -n {}-a infer_NJ_nucleotide_pair.mao -d {} -o {}'.format(v, pin, pout)
	subprocess.call(cmd.split())

	# Megacc replaces '^' with an underscore (unlike other characters that do cause problems ...)
	# use RE to replace the _ introduced by megacc & overwrite file
	with open(treefn, 'r') as f:
		filestr = f.read()

	reg = '_(?=[0-9]+:[0-9].[0-9]+)'
	filestr = re.sub(reg, '^', filestr)

	with open(treefn, 'w') as out:
		out.write(filestr)


# optional/possible: return Tree object (for msa viewer ?)
#(opt.?) ToDo: addition of removed files as extra cluster group (+ [recursive?] grouping of dropped files)
def treeClusters(protein, threshold):
	'''
	Use pw-similarity to divide a clustered MSA in different groups

	:uses: {protein}_promoter_matrix.tsv, {protein}_promoters_megacc.nwk
	:param protein: protein name from phylogenetics workflow
	:param threshold: similarity threshhold for generation of groups
	:creates: {protein}_promoters_groups.csv
	'''

	treefile = os.path.join(PATH, protein + '_promoters_megacc.nwk')
	matrixfile = os.path.join(PATH, protein + '_promoters_matrix.tsv')

	if not os.path.isfile(treefile):
		print('Newick-tree file not found. Run step 3 for this file to be created.')
		sys.exit()
	if not os.path.isfile(matrixfile):
		print('pws matrix file not found. Run step 2 for this file to be created.')
		sys.exit()

	#These should be specified by a config file or flag or whatever
	meanfunc1 = gmean
	meanfunc2 = gmean

	t = Tree(treefile)

	# matrix conatins also the dropped elements, even though probably/not yet needed
	# if VERBOSITY:
	# 	print()

	pws = pd.read_csv(matrixfile, sep='\t', header=1, index_col=0, na_values='-')

	for node in t.traverse('postorder'):
		# Only leaves are named in a newick file!
		node.add_features(sim=None)

		if node.is_leaf():
			# all node.names are acc^tax

			node.add_features(cl=None)
			node.add_features(distances = pws[node.name])

		else:
			# differentiate between children that are leaves (get mean of values for all leaves under current node from pws-matrix-slice stored on leaf)
			# and children that are not (they have [recursive mean of children] similarity values )
			leaves = node.get_leaf_names()

			values = []

			for n in node.children:
				if n.name in leaves:
					meansim = meanfunc1(n.distances[leaves].dropna())
					if not np.isnan(meansim):
						values.append(meansim)
				else:
					values.append(n.sim)

			node.sim = meanfunc2(values)

	# Cluster nodes based on similarity threshold
	# From bottom up (meaning a cluster will contain all leaves under the point where it was created)
	# If a node has sim lower than th, put each child-branch in a separate cluster
	clusters = []
	for node in t.traverse('postorder'):
		if node.is_leaf():
			continue

		if node.sim < threshold or node.is_root():
			for x in node.iter_descendants():
				# checking just for l.cl will/may overwrite the first generated cluster (l.cl = 0)
				if any(isinstance(l.cl, int) for l in x.iter_leaves()):
					continue

				# add a list with all node instances to clusters
				clusters.append(x.get_leaves())
				for leaf in clusters[-1]:
					leaf.cl = len(clusters)-1

	# go through clusters, >=2 neighbors have len <= 2, fuse them
	cleanclusters = []
	current = []
	merged = []
	for cluster in clusters:
		# Number by variable ??
		if len(cluster) <= 2:
			current += cluster
		else:
			if current:
				merged.append(len(cleanclusters)) #-> index of new (ex-)small (merged) cluster
				cleanclusters.append(current)
				current = []
			cleanclusters.append(cluster)
	# if there's still an open merging cluster left 'close it'
	if current:
		merged.append(len(cleanclusters))  # -> index of new (ex-)small (merged) cluster
		cleanclusters.append(current)

	with open(os.path.join(PATH, protein+'_promoters_groups.csv'), 'w') as clout:
		clout.write('# Clusters for {} determined with similarity threshhold {} on {}.\n'.format(os.path.basename(matrixfile), threshold, strftime('%x')))

		for i, cluster in enumerate(cleanclusters):
			if i in merged:
				note = ' (remerged)'
			else:
				note = ''

			clout.write('!Cluster {}{}\n'.format(i, note))
			clout.write(','.join([n.name for n in cluster]) + '\n')

			#For tests / returning Tree
			# for leaf in cluster:
			# 	leaf.cl = i # updating to clean cluster numbers
			# 	if i in merged:
			# 		leaf.add_features(remerged=True)
			# 	else:
			# 		leaf.add_features(remerged=False)


	# #simple visualisation for tests
	# def layout(node):
	# 	if node.is_leaf():
	# 		faces.add_face_to_node(AttrFace("cl"), node, column=0)
	# 	else:
	# 		faces.add_face_to_node(AttrFace("sim"), node, column=0)
	# ts = TreeStyle()
	# ts.layout_fn = layout
	# t.show(tree_style=ts)
	# # test vis end

	# for passing to next function (msa viewer) - maybe not needed
	# return t


def msaview(protein):

	'''
	Visualise the MSA & the groups of the megaCC clustering.
	Adapted from msa_viewer.py

	:uses: {protein}_promoters_groups.csv, {protein}_promoters_kept.fasta
	:param protein: protein name from phylogenetics workflow
	:creates: {protein}_grouped_msa.png
	'''

	alignmentFile = os.path.join(PATH, protein + '_promoters_kept.fasta')
	groupfn = os.path.join(PATH, protein + '_promoters_groups.csv')
	outfn = os.path.join(PATH, protein + '_grouped_msa.png')

	if not os.path.isfile(alignmentFile):
		print('Sorted MSA fasta file not found. Run step 2 for this file to be created.')
		sys.exit()
	if not os.path.isfile(groupfn):
		print('Grouping file not found. Run step 4 for this file to be created.')
		sys.exit()

	dnacolors = {
		'A': (203, 0, 0),
		'T': (0, 0, 203),
		'U': (0, 0, 203),
		'G': (0, 203, 0),
		'C': (203, 0, 191),
		'N': (0, 0, 0),
		'-': (255, 255, 255),
		# these are sometimes added by clustalo aswell
		'R': (0, 0, 0),  # puRine, A/G
		'Y': (0, 0, 0),  # pYrimidine, C/T
		'W': (0, 0, 0),
		'K': (0, 0, 0),
		'M': (0, 0, 0),
		'S': (0, 0, 0)	}

	clustercolors = [(166, 206, 227), (31, 120, 180), (178, 223, 138), (51, 160, 44), (251, 154, 153), (227, 26, 28),
					 (253, 191, 111), (255, 127, 0), (202, 178, 214), (106, 61, 154), (255, 255, 179), (255, 234, 77), (177, 89, 40), (0, 0, 0)]

	data = {}
	current = None

	with open(alignmentFile) as f:
		for line in f:
			line = line.rstrip()

			if line.startswith('>'):
				current = line.split()[0][1:]
				data[current] = []
			else:
				for char in line:
					data[current].append(dnacolors[char])

	groupcolors = {}
	order = []
	with open(groupfn, 'r') as f:
		current = 0
		for line in f:
			line = line.split('#')[0].rstrip()
			if not line:
				continue

			# Allows definition of cluster names in file
			if line.startswith('!'):
				continue

			cluster = line.split(',')
			order = order + cluster

			for elem in cluster:
				groupcolors[elem] = clustercolors[current % len(clustercolors)]

			current += 1

	height = len(order)
	width = len(next(iter(data.values())))
	borders = 10

	im = Image.new('RGB', (width + 2 * borders, height), (255, 255, 255))

	for y, elem in enumerate(order):
		for x in range(borders):
			im.putpixel((x, y), groupcolors[elem])
			im.putpixel((x + width + borders, y), groupcolors[elem])
		try:
			for x, value in enumerate(data[elem]):
				im.putpixel((x + borders, y), value)
		except KeyError:
			print(elem, alignmentFile)
			raise

	im.save(outfn)


def _makeTree(quality_level=1):
	'''
	Create Newick tree with all species represented in refseq_genomes, based on level of genome quality
	'''

	if not os.path.isfile('refseq_taxids_l{}.txt'.format(quality_level)):
		taxids(quality_level)

	Sani = NodeSanitizer()
	TF = TaxFinder()

	outfn = 'refseq_tree_l{}.tre'.format(quality_level)

	#maybe interesting
	#accs= set()

	# While this file may contain multiple entries for one taxid, this should not impact the tree generation

	with open('refseq_taxids_l{}.txt'.format(quality_level), 'r') as infile:
		next(infile)
		lineages = []

		for line in infile:
			tax, acc, dummy = line.split('\t')
			#accs.add(acc)
			li = TF.getLineage(int(tax), display='both')
			if li:
				lineages.append(li)
			else:
				print('Could not find TaxID `{}` in TaxFinder.'.format(tax))

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


def grouptree(protein, skipclusters, genomequality):
	'''
	wrapper to call phylotree.py and plot cluster distribution on tree

	:param protein:
	:param skipclusters:
	'''

	groupfn = os.path.join(PATH, protein + '_promoters_groups.csv')
	outfn = os.path.join(PATH, protein + '_promoters_grouptree.pdf')
	treefn = 'refseq_tree_l{}.tre'.format(genomequality)


	if not os.path.isfile(groupfn):
		print('Grouping file not found. Run step 4 for this file to be created.')
		sys.exit()

	if not os.path.isfile(treefn):
		_makeTree(genomequality)

	prunefn = 'tree_to_prune.txt'


	drop = set(int(i) for i in skipclusters if i.isdecimal())
	if len(drop) != len(skipclusters):
		skipnames = [i for i in skipclusters if i.isalpha()]
		with open(groupfn, 'r') as f:
			i = 0
			for line in f:
				line = line.split('#')[0].rstrip()
				if not line:
					continue

				if line.startswith('!'):
					for n in skipnames:
						if n in line[1:]:
							drop.add(i)
				else:
					i += 1

	t = Tree(treefn, format=8)
	# tree, treefile, prune, sqrt = False, collapse = True, show_empty_species = True, startnode = None, countTotal = False
	TM = Clusters(tree=t, treefile=treefn, prune=prunefn, filename=groupfn, dropclusters=drop, show_empty_species=False)
	TM.renderTree(outfn, 1000, 90)


def motifAnalysis(protein, skipclusters):
	'''Generate motifs of each group within the MSA

	:uses: {protein}_promoters_groups.csv, {protein}_promoters_aligend.fasta
	:param protein: protein name from phylogenetics workflow
	:creates: motifs/*
	'''

	fastafn = os.path.join(PATH, protein + '_promoters_aligned.fasta')
	groupfn = os.path.join(PATH, protein + '_promoters_groups.csv')

	if not os.path.isfile(fastafn):
		print('MSA fasta file not found. Run step 1 for this file to be created.')
		sys.exit()
	if not os.path.isfile(groupfn):
		print('Grouping file not found. Run step 4 for this file to be created.')
		sys.exit()

	#dict like object fasta-headers as keys
	seqs = SeqIO.index(fastafn, 'fasta', alphabet=DNAAlphabet())

	skipnames = [i for i in skipclusters if i.isalpha()]

	groups = {}
	with open(groupfn, 'r') as f:
		i = 0
		name = ''
		for line in f:
			line = line.split('#')[0].rstrip()
			if not line:
				continue

			if str(i) in skipclusters:
				if not line.startswith('!'):
					i += 1
				continue

			# Allows definition of cluster names in file
			if line.startswith('!'):
				name = line[1:]
				continue

			cluster = line.split(',')

			if len(cluster) == 1:
				if VERBOSITY:
					print("Skipped '{}' because there's only one entry.".format(name or 'cluster{:02d}'.format(i)))
				name = ''
				continue

			#allow to drop clusters based on partial name (like 'remerged')
			if name and any((i in name for i in skipnames)):
				continue
			elif name:
				groups[name] = [seqs[n] for n in cluster]
				name = ''
			else:
				groups['cluster{:02d}'.format(i)] = [seqs[n] for n in cluster]

			i += 1

	motifdict = {}

	#function to get score for a window of n sequences (no gaps)
	def groupsimilarity(seqs):
		score = 0
		for pos in zip(*seqs):
			# # score +1 only if 100% conserved
			# score += pos.count(pos[0]) == len(pos)

			# score +1 if at least x % conserved
			score += pos.count(pos[0])/float(len(pos)) >= 0.8

			# # score +x%
			# score += pos.count(pos[0]) / float(len(pos))
		return score

	for name, gr in groups.items():
		# split group into sections with high similarity

		name = name.replace(' ', '')

		# First remove ('conserved') gaps from the group - they are irrelevant & make problems
		# remember positions of gaps to later on adapt the motif positions for the fasta-file
		seqs = [''] * len(gr)
		gappos = []
		gapcount = 0

		for i, pos in enumerate(zip(*[str(rec.seq) for rec in gr])):
			if pos.count('-') == len(pos):
				gapcount += 1
			else:
				gappos.append(gapcount)
				for j, base in enumerate(pos):
					seqs[j] += base

		motifpos = []
		winlen = 40
		scoremin = 35 #meaning this number of conserved bases inside the window
		start = None

		for i in range(len(seqs[0]) - winlen):
			window = [seq[i:i+winlen] for seq in seqs]
			# if ANY gap/non base char [N etc] in window, close current motif
			# opt: be less restrictive here, allow 1-x % gaps ?
			if set('ATGC') != set(''.join(window)):
				if start is not None:
					motifpos.append([start, i - 1 + winlen])
					start = None
				continue

			if start is None and groupsimilarity(window) >= scoremin:
				start = i

			elif start is not None and groupsimilarity(window) <= scoremin:
				motifpos.append([start, i - 1 + winlen])
				start = None


		# A lot of the position pairs overlap
		# go through them and fuse as necessary
		abs_overlap = 3
		per_overlap = 0.66

		if len(motifpos) > 1:
			clearpos = []
			current = motifpos[0][:]
			for start, end in motifpos[1:]:
				prev_start, prev_end = current
				#check for overlap, fuse if
				# - end of second/current hits is shifted by max `abs_overlap` to first
				# - overlap is the major part (>= per_overlap) of either hit
				if start < prev_end and (end - prev_end <= abs_overlap  or (prev_end - start) / float(end - start) >= per_overlap or (prev_end - start) / float(prev_end - prev_start) >= per_overlap):
					current[1] = end
				# No or not sufficient overlap, save previous, start with new current
				else:
					clearpos.append(current)
					current = [start, end]
			# save last current
			clearpos.append(current)
		else:
			clearpos = motifpos

		# make motif objects
		for i, pos in enumerate(clearpos):
			motifdict[name+'-{:02d}'.format(i)] = motifs.create([seq[pos[0]:pos[1]] for seq in seqs], alphabet=IUPAC.unambiguous_dna)
			# the attributes 'matirx_id' and 'name' are internally used by the Bio.motifs.jaspar module and are written to the header
			# save motif name (as for seq-logo files) & position in fasta file to header
			setattr(motifdict[name + '-{:02d}'.format(i)], 'matrix_id', name + '-{:02d}'.format(i))
			setattr(motifdict[name + '-{:02d}'.format(i)], 'name', 'Start/End (in {}): {}/{}'.format(fastafn, pos[0]+gappos[pos[0]], pos[1]+gappos[pos[1]]))

	motiffolder = 'motifs'
	# results with different parameters may yield less files, therefore all old files should be removed
	if os.path.isdir(os.path.join(PATH, motiffolder)):
		subprocess.call(['rm', '-r', os.path.join(PATH, motiffolder)])

	os.makedirs(os.path.join(PATH, motiffolder))


	#Opt Todo: make a single file with all consensus sequences?!
	for name, motif in motifdict.items():
		#the heigth (bit value) in logos doesnt seem quite correct currently
		motif.weblogo(os.path.join(PATH, motiffolder, '{}_{}_seqlogo.png'.format(protein, name)))

	for name, gr in groups.items():
		grmotifs = sorted([v for k, v in motifdict.items() if name.replace(' ', '') == k.split('-')[0]], key=lambda x: x.matrix_id)

		#dont make a file for groups where no conserved motif was created
		if not len(grmotifs):
			continue

		with open(os.path.join(PATH, motiffolder, '{}_{}_motifs.jas'.format(protein, name)), 'w') as out:
			#header is '>{matrix_id} {name}'
			out.write(motifs.write(grmotifs, 'jaspar'))



tasknames = ['acquire', 'clustal', 'autosort', 'megacc', 'grouping', 'msaview', 'grouptree', 'motifs']

tasks = {'acquire': ('find and download promoter sequences', getPromoterSeqs, ['evalue', 'length', 'genomequality', 'phylopath']),
		 'clustal': ('run MSA with clustalo', runClustal, []),
		 'autosort': ('similarity sort MSA for clustering', autosort, ['drop']),
		 'megacc': ('run MSA clustering with Megacc', runMegacc, []),
		 'grouping': ('find groups in clustered MSA', treeClusters, ['threshold']),
		 'msaview': ('Make plot of the MSA including the groups', msaview, []),
		 'grouptree': ('Show distributions of the cluster over a tree.', grouptree, ['skipclusters', 'genomequality']),
		 'motifs': ('analyse MSA cluster groups for motifs', motifAnalysis, ['skipclusters'])}

def runWorkflow(addargs, start, end=''):
	'''
	Starts the workflow from `start` until `end`. If `end` is empty, the workflow is run until the last element of the workflow.

	:param addargs: args objects from argparser with flags for differents steps
	:param start: The name/number of the first element to run.
	:param end: The name/number of the last element to run or an empty string if the workflow shall be run until the end.
	'''

	if isinstance(start, int) or start.isdecimal():
		if int(start) >= len(tasknames):
			raise ValueError('Only {} tasks exist, number {} is too high.'.format(len(tasknames), int(start)))
		else:
			startidx = int(start)
	elif start not in tasknames:
		raise ValueError('{} is no valid task.'.format(start))
	else:
		startidx = tasknames.index(start)

	if isinstance(end, int) or end.isdecimal():
		if int(end) >= len(tasknames):
			raise ValueError('Only {} tasks exist, number {} is too high.'.format(len(tasknames), int(end)))
		else:
			endidx = max(int(end)+1, startidx+1)
	elif end == '':
		endidx = len(tasknames)
	elif end not in tasknames:
		raise ValueError('{} is no valid task.'.format(end))
	else:
		endidx = tasknames.index(end) + 1


	for taskname in tasknames[startidx:endidx]:
		if VERBOSITY:
			print('{}: "{}"'.format(taskname, tasks[taskname][0]))
		task = tasks[taskname][1]
		kwargs = {'protein': addargs.protein}
		kwargs.update({x: getattr(addargs, x) for x in tasks[taskname][2]})
		task(**kwargs)


if __name__ == '__main__':
	import argparse
	import textwrap

	workflow = '\n'.join(textwrap.wrap("""The following is a list of the workflow. The names or numbers can be used for the -s, -e or -o arguments.""",	width=80))

	workflow += '\n\n' + '\n'.join(('{:>2}. {:<10} {}'.format(i, name, tasks[name][0]) for i, name in enumerate(tasknames)))

	parser = argparse.ArgumentParser(description='Script for analysis of Promoters. Run for single protein (from Phylogenetics workflow) or use config file')

	reqgroup = parser.add_mutually_exclusive_group(required=True)

	reqgroup.add_argument('-p', '--protein', help='Run the workflow with a proteiname from phylogenetics')
	reqgroup.add_argument('-c', '--config', action='store_true', help='Run the workflow using the config file, may overwrite all other flags')
	reqgroup.add_argument('-l', '--list', action='store_true', help='Shows the whole workflow with information and exits')
	reqgroup.add_argument('-i', '--init', action='store_true', help='Initiates the working directory with necessary files and folders')

	parser.add_argument('-v', '--verbose', action='store_true', help='Be verbose.')
	parser.add_argument('-s', '--startfrom', default=0, help='Run from and including this step') # [e.g. 7 or hist]
	parser.add_argument('-e', '--end', default='', help='Stop after this step, ignored if smaller than start') # [e.g. 7 or hist]
	parser.add_argument('-o', '--only', default='', help='Run only the given step, takes precedence over -s and -e') # [e.g. 4 or unique]
	parser.add_argument('-pa', '--path', default='', help='Path Prefix which will be added to all generated files & folders')

	# only step 0 - loading
	group0 = parser.add_argument_group('0.', 'Optional arguments applied during step 0: loading.')

	group0.add_argument('--evalue', default=50, type=int, help='Lowest exponent allowed as evalue during hit-collection. Default: 50')
	group0.add_argument('--length', default=1500, type=int, help='Number of bases in front of CDS to be acquired as promoter sequence. Default: 1500')
	group0.add_argument('-q', '--genomequality', default=1, type=int, help='Quality level for refseq genomes (also in step 6). 0: all, 1: atleast chromosome, 2: representative/reference')
	# Optional: reconfigure / add extra option to use specific file instead
	group0.add_argument('--phylopath', default='', help='Specify the filepath for phylogenetics workflow')

	# only step 2 - autosort
	group2 = parser.add_argument_group('2.', 'Optional arguments applied during step 2: autosort.')

	group2.add_argument('-d', '--drop', nargs='+', default=[], help='Remove these Accs (or index_from_0) during sorting')

	#only step 4 - grouping
	group4 = parser.add_argument_group('4.', 'Optional arguments applied during step 4: grouping.')

	group4.add_argument('-t', '--threshold', type=float, default=0.5, help='Similarity threshold for grouping of clustered MSA. Default: 0.5')

	#only step 6/7
	group6 = parser.add_argument_group('6/7.', 'Optional arguments applied during step 6 & 7: grouptree, motif.')

	group6.add_argument('-sk', '--skipclusters', nargs='+', default=[], help="Don't analyse these groups (ints with index_from_0 or parts of group names)")


	args = parser.parse_args()

	if args.list:
		parser.print_help()
		print('')
		print(workflow)
		sys.exit()

	if args.init:
		init()
		sys.exit()

	# check if init has been run

	if args.config:
		readConfig(args)
		sys.exit()

	PATH = args.path
	VERBOSITY = args.verbose

	if PATH:
		os.makedirs(PATH, exist_ok=True)

	if args.only:
		runWorkflow(args, args.only, args.only)
	else:
		runWorkflow(args, args.startfrom, args.end)