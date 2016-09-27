#!/usr/bin/env python3
import pandas as pd
import numpy as np

def sortfasta(fname, drop, verbose):
	with open(fname, 'r') as f:
		fastas = ['>'+fasta for fasta in f.read().split('>')][1:]
	dontkeep = set(drop)
	with open(fname.split('.fasta')[0]+'_sorted.fasta', 'w') as out:
		for i, fasta in enumerate(fastas):
			flines = fasta.split('\n')
			if i+1 in dontkeep:
				if verbose:
					print('Dropped:', flines[0])
				continue
			if any(base in ''.join(flines[1:]) for base in ('R', 'Y', 'K', 'W', 'S', 'M', 'N')):
				if verbose:
					print('Bad base char in Nr {}: {}'.format(i+1, flines[0]))
				continue
			out.write(fasta)


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
			if i==0 and j ==3:
				print(d, names[i], names[j])
	mat.index = names
	mat.columns = names
	return mat


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
		return (score/len(seq1))
	else:
		return np.NaN



def autosort(fname, verbose):
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

	if verbose:
		print('Calculating pairwise distance matrix for {}.'.format(fname))

	matrix = _pw_distance(seqs, names)

	dropped = []
	# (true/false NaN for each cell).(count per col).(count total)
	while matrix.isnull().sum().sum():
		nans = matrix.isnull().sum()
		#Dropping all columns with the highest number of NaNs + the respective rows
		worst = nans[nans == max(nans)].index
		matrix.drop(worst, axis=0, inplace=True)
		matrix.drop(worst, axis=1, inplace=True)
		dropped += list(worst)

	if verbose:
		print('Following IDs were dropped to remove N/A from the distance matrix:')
		print(', '.join(dropped))

	keep = list(matrix.index)

	with open(fname, 'r') as f, open(fname.split('.fasta')[0]+'_kept.fasta', 'w') as out, open(fname.split('.fasta')[0]+'_dropped.fasta', 'w') as out2:
		for line in f:
			if line.startswith('>'):
				keeping = line.rstrip().split('|')[1] in keep

			if keeping:
				out.write(line)
			else:
				out2.write(line)








if __name__ == '__main__':
	import argparse

	parser = argparse.ArgumentParser(description='Script to sort out certain entries of a (MSA) fasta file.')

	parser.add_argument('fasta', help='Filename of a fasta to sort through')

	parser.add_argument('-v', '--verbose', action='store_true', help='Be verbose.')

	group = parser.add_mutually_exclusive_group(required=True)
	group.add_argument('-n', '--numbers', nargs='+', type=int, help='Index numbers of entries to drop (first index = 1). Only saves non-dropped files.')
	group.add_argument('-a', '--auto', action='store_true', help='Automatic calculation of entries to drop by calculating pw distance matrix. Saves dropped files aswell.')

	args = parser.parse_args()

	if args.auto:
		autosort(args.fasta, args.verbose)
	else:
		sortfasta(args.fasta, args.numbers, args.verbose)