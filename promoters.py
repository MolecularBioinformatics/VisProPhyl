#!/usr/bin/env python3
import os
import sys
from ete3 import Tree, TreeStyle, faces, AttrFace
import pandas as pd
import numpy as np
from scipy.stats.mstats import gmean
from time import strftime

# For now only a few functions
# However, this script is intended to become a complete workflow (like phylogenetics)
# This will (at the least) require implementation of GetPromoterSeqs & automated calls of clustalo and megacc
# Either integrate msa-viewer functions or import them (or make the scripts work well after each other)

# ToDo: make & read config file
#	- makes all flags (extra)optional OR overwrites all flags OR only used with flag
#	- all proteins that should be used
# 	- make & use appropiated subfolders
#	- functions (& threshhold) used in treeClusters
#	- numbers & normalisation used in _similarity
#	- genome quality level (getPromoterSeq & treefile for later)

# TODO: change all scripts (really all of them) to always use acc^tax (or another consequent format) [instead of saving the tax information in the annotation files]
# This also/especially includes the used fasta-headers in GetPromoterSeq

def _myGetBasename(name):
	'''
	Strips the path and the extension from a filename. E.g. /a/path/to/file.png -> file

	:param name: The file and path to be processed
	:returns: The processed filename
	'''

	return os.path.splitext(os.path.basename(name))[0]


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


def autosort(basename, drop, verbosity=True):
	'''
	Take a fasta file with multiple aligned sequences & remove entries that have no similarity score with other entries (greedy alg.)

	:param basename: basename for .fasta file with aligned sequences
	:param drop: list (or similar) of entries that are to be dropped regardless of score
	:param verbosity: Bool
	:create: (pw-similarity) _matrix.tsv, _kept.fasta, _dropped.fasta
	'''
	seqs = []
	names = []
	fname = basename + '.fasta'

	with open(fname, 'r') as f:
		current = ''
		for line in f:
			line = line.rstrip()
			if line.startswith('>'):
				if current:
					seqs.append(current)
				names.append(line.split('|')[1])
				current = ''
				continue

			if line:
				current += line
		seqs.append(current)

	if verbosity:
		print('Calculating pairwise similarity matrix for {}.'.format(fname))

	matrix = _pw_similarity(seqs, names)

	with open(basename + '_matrix.tsv', 'w') as out:
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

	if verbosity:
		print('Following IDs were dropped to remove N/A from the similarity matrix:')
		print(', '.join(dropped))

	keep = list(matrix.index)

	i = 0

	with open(fname, 'r') as f, open(basename+'_kept.fasta', 'w') as out, open(basename+'_dropped.fasta', 'w') as out2:
		for line in f:
			if line.startswith('>'):
				name = line.rstrip().split('|')[1]
				if name in drop or str(i) in drop:
					keeping = False
				else:
					keeping = name in keep
				i += 1

			if keeping:
				out.write(line)
			else:
				out2.write(line)


# optional/possible: return Tree object (for msa viewer ?)

#(opt.?) ToDo: addition of removed files as extra cluster group (+ [recursive?] grouping of dropped files)

#opt. Todo: reclustering of merged groups by 'original' method (see/c&p msa_viewer)
def treeClusters(basename, threshold=None, verbosity=True):
	'''
	Use pw-similarity to divide a clustered MSA in different groups

	:param treefile: Newick tree with clustered MSA
	:param matrixfile: filename, pw similarity matrix (.tsv)
	:param threshold: similarity threshhold for generation of groups
	:param verbosity: bool
	:creates: csv-file with clusters
	'''

	# update treefile name once a wrapper function for megacc is there to be more uniform
	treefile = basename.split('_')[0] + '_megacc.nwk'
	matrixfile = basename + '_matirx.tsv'


	if threshold is None:
		#read from config
		threshold = 0.5

	#These should be specified by a config file or flag or whatever
	meanfunc1 = gmean
	meanfunc2 = gmean

	t = Tree(treefile)

	# matrix conatins also the dropped elements, even though probably/not yet needed

	pws = pd.read_csv(matrixfile, sep='\t', header=1, index_col=0, na_values='-')

	for node in t.traverse('postorder'):
		# Only leaves are named in a newick file!
		node.add_features(sim=None)

		if node.is_leaf():
			# names in tree are 'promoter|acc' (equal to fasta header)
			# names in matix are only acc
			# This should be changed at some point
			node.name = node.name.split('|')[1]

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
				if any(l.cl for l in x.iter_leaves()):
					continue

				# add a list with all node instances to clusters
				# Alternative (currenlty needed in msa_viewer), add: 'promoter|' + node.name to said lists
				clusters.append(x.get_leaves())

				#can be deleted since, merging overwrites this anyway, only needed for tests / direct passing of tree
				for leaf in clusters[-1]:
					leaf.cl = len(clusters)-1

	# go through clusters, >=2 neighbors have len <= 2, fuse them
	cleanclusters = []
	current = []
	merged = []
	for cluster in clusters:
		if len(cluster) <= 2:
			current += cluster
		else:
			if current:
				merged.append(len(cleanclusters)) #-> index of new (ex-)small (merged) cluster
				cleanclusters.append(current)
				current = []
			cleanclusters.append(cluster)


	with open(basename+'_clustering.csv', 'w') as clout:
		clout.write('# Clusters for {} determined with similarity threshhold {} on {}.\n'.format(os.path.basename(matrixfile), threshold, strftime('%x')))

		for i, cluster in enumerate(cleanclusters):
			if i in merged:
				note = ' (remerged)'
			else:
				note = ''

			clout.write('!Cluster {}{}\n'.format(i, note))
			clout.write(','.join(cluster))

			#For tests / returning Tree
			for leaf in cluster:
				leaf.cl = i
				# if i in merged:
				# 	leaf.add_features(remerged=True)
				# else:
				# 	leaf.add_features(remerged=False)

	#simple visualisation for tests
	def layout(node):
		if node.is_leaf():
			faces.add_face_to_node(AttrFace("cl"), node, column=0)
		else:
			faces.add_face_to_node(AttrFace("sim"), node, column=0)
	ts = TreeStyle()
	ts.layout_fn = layout
	t.show(tree_style=ts)
	# test vis end

	# for passing to next function (msa viewer) - maybe not needed
	# return t


def motifAnalysis(basename, verbosity=True):

	fastafn = basename + '.fasta'
	clusterfn = basename + '_clustering.csv'

	# use SeqIO / Bio-something parser on (complete) msa fasta
	# split into groups according to cluster file (dict with !names as keys)

	# use Bio.motifs on each group

	# generate sequence logo, consensus sequencce & save motive in Jaspar format





tasknames = ['autosort', 'grouping', 'motifs']

tasks = {
	#	  'getPromSeq': ('get promoter seqs', wrapper_func, []),
	#	  'runClustal': ('run MSA with clustalo', clustalfunc, []),
		 'autosort': ('sort MSA for clustering', autosort, ['drop']),
	#	  'runMegacc': ('run MSA clustering via megacc', megaccfunc, []),
		 'grouping': ('find groups in clustered MSA', treeClusters, ['threshold']),
		 'motifs': ('analyse MSA cluster groups for motifs', motifAnalysis, [])}


def runWorkflow(addargs, start, end=''):
	'''
	Starts the workflow from `start` until `end`. If `end` is empty, the workflow is run until the last element of the workflow.

	:param addargs: args objects from argparser with flags for differents steps
	:param start: The name of the first element to run.
	:param end: The name of the last element to run or a falsy value if the workflow shall be run until the end.
	'''

	if start not in tasknames:
		raise ValueError('{} is no valid task.'.format(start))

	if end == '':
		endidx = len(tasknames)
	elif end not in tasknames:
		raise ValueError('{} is no valid task.'.format(end))
	else:
		endidx = tasknames.index(end) + 1

	startidx = tasknames.index(start)
	for taskname in tasknames[startidx:endidx]:
		print('{}: "{}"'.format(taskname, tasks[taskname][0]))
		task = tasks[taskname][1]
		kwargs = {'basename': addargs.filename, 'verbosity': addargs.verbose}
		kwargs.update({x: getattr(addargs, x) for x in tasks[taskname][2]})
		task(**kwargs)



if __name__ == '__main__':
	import argparse

	parser = argparse.ArgumentParser(description='Script for analysis of Promoters. Run for single protein (from Phylogenetics workflow) or use config file')

	# add mutually exclusive group (required = True)
	# filename & config into group
	parser.add_argument('filename', help='base filename for all files')

	parser.add_argument('-v', '--verbose', action='store_true', help='Be verbose.')
	# parser.add_argument('-l', '--list', action='store_true', help='Shows the whole workflow with information and exits')
	# parser.add_argument('-i', '--init', action='store_true', help='Initiates the working directory with necessary files and folders')
	parser.add_argument('-s', '--startfrom', default=0, help='Run from and including this step') # [e.g. 7 or hist]
	parser.add_argument('-e', '--end', default='', help='Stop after this step') # [e.g. 7 or hist]
	parser.add_argument('-o', '--only', default='', help='Run only the given step, takes precedence over -s and -e') # [e.g. 4 or unique]

	# only step 1 - autosort
	parser.add_argument('-d', '--drop', nargs='+', default=[], help='Remove these Accs (or index_from_0) during sorting')

	#only step 2 - grouping
	parser.add_argument('-t', '--threshold', type=float, default=None, help='Similarity threshold for grouping of clustered MSA')

	args = parser.parse_args()

	# if args.list:
	# 	parser.print_help()
	# 	print('')
	# 	print(workflow)
	# 	sys.exit()

	# if args.init:
	# 	init()
	# 	sys.exit()

	# check if init has been run


	# if config: start parser


	if args.only:
		try:
			a = int(args.only)
		except ValueError:
			if args.only in tasknames:
				runWorkflow(args, args.only, args.only)
			else:
				parser.print_help()
		else:
			if a >= len(tasknames):
				parser.print_help()
	else:
		try:
			a = int(args.startfrom)
		except ValueError:
			if args.startfrom in tasknames:
				pass
			else:
				parser.print_help()
				sys.exit()
		else:
			if a >= len(tasknames):
				parser.print_help()
				sys.exit()
		if args.end:
			try:
				a = int(args.end)
			except ValueError:
				if args.end in tasknames:
					pass
				else:
					parser.print_help()
					sys.exit()
			else:
				if a >= len(tasknames):
					parser.print_help()
					sys.exit()

		runWorkflow(args, args.startfrom, args.end)

