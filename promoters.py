#!/usr/bin/env python3
import os
import sys
from ete3 import Tree, TreeStyle, faces, AttrFace
import pandas as pd
import numpy as np
from scipy.stats.mstats import gmean

# For now only a few functions
# However, this script is intended to become a complete workflow (like phylogenetics)
# This will (at the least) require impelementation of GetPromoterSeqs, sort_fastas & automated calls of mclustalo and megacc
# Either integrate msa-viewer functions or import them (or make the scripts work well after each other)

# ToDo: make & read config file
#	- all proteins that should be used
# 	- make & use appropiated subfolders

# TODO: change all scripts (really all of them) to always use acc^tax (or another consequent format) [instead of saving the tax information in the annotation files]


def _pw_distance(seqs, names):
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

# The specific numbers here will have a huge impact on downstream results
# In the context of an MSA the normalisation of the score to len(seq) is pointless, therefore removed
# Consider: Normalisation to len without gaps (not sure if thats a good idea though)
def _similarity(seq1, seq2):
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


def autosort(fname, drop, verbosity=True):
	seqs = []
	names = []
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
		print('Calculating pairwise distance matrix for {}.'.format(fname))

	matrix = _pw_distance(seqs, names)

	if '.fa' in fname:
		mname = fname.split('.fa')[0] + '_matrix.tsv'
	else:
		mname = fname.split('.')[0] + '_matrix.tsv'
	with open(mname, 'w') as out:
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
		print('Following IDs were dropped to remove N/A from the distance matrix:')
		print(', '.join(dropped))

	keep = list(matrix.index)

	i = 0

	with open(fname, 'r') as f, open(fname.split('.fasta')[0]+'_kept.fasta', 'w') as out, open(fname.split('.fasta')[0]+'_dropped.fasta', 'w') as out2:
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


def treeClusters(treefile, matrixfile, verbosity=True):

	#These should be specified by a config file or flag or whatever
	meanfunc1 = gmean
	meanfunc2 = gmean

	t = Tree(treefile)

	# matrix conatins also the dropped elements, even though probably/not yet not needed
	pws = pd.read_csv(matrixfile, sep='\t', header=1, index_col=0, na_values='-')

	## !discarded (would be needed if no data is written on nodes & instead always pulled from matrix)
	# Drop upper half of (symmetric!) matrix to get rid of duplicates
	# pws.values[np.triu_indices_from(pws)] = np.nan
	##

	for node in t.traverse('postorder'):
		# Only leaves are named!

		node.add_features(sim=None)

		if node.is_leaf():
			# names in tree are 'promoter|acc'
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

	# Cluster nodes based on similarity threshhold
	# From bottom up (meaning a cluster will contain all leafes under the point where it was created)
	# If a node has sim lower than th, put each child-branch in a separate cluster
	clusters = []
	threshhold = 0.5 #test-value, has to be a vaibale obsvly
	for node in t.traverse('postorder'):
		if node.is_leaf():
			continue

		if node.sim < threshhold or node.is_root():
			for x in node.iter_descendants():
				if any(l.cl for l in x.iter_leaves()):
					continue

				# add a list with all node instances to clusters
				# Alternative (currenlty needed in msa_viewer), add: 'promoter|' + node.name to said lists
				clusters.append(x.get_leaves())
				#can be deleted since, merging overwrites this anyway
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

	# either write feature on node to indicate 'merged' or write this to a file
	for i, cluster in enumerate(cleanclusters):
		for leaf in cluster:
			leaf.cl = i





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




# flags are very very much in test state
if __name__ == '__main__':
	import argparse

	parser = argparse.ArgumentParser(description='Script to sort out certain entries of a (MSA) fasta file.')

	parser.add_argument('fasta', help='Filename of a fasta to sort through')

	parser.add_argument('-v', '--verbose', action='store_true', help='Be verbose.')
	parser.add_argument('-d', '--drop', nargs='+', help='Never keep these (accIDs or index_from_0)')


	parser.add_argument('-a', '--auto', action='store_true', help='Automatic calculation of entries to drop by calculating pw distance matrix. Saves dropped files aswell.')

	#parser.add_argument('-c', '--tree_matrix', nargs=2, help='tree file & matrix file')

	args = parser.parse_args()
	if args.auto:
		autosort(args.fasta, args.drop, args.verbose)
		treeClusters('Promoters/clk_megacc.nwk', 'Promoters/clk_promoters_aligned_matrix.tsv')
	#else:
		#treeClusters(*args.tree_matrix)
