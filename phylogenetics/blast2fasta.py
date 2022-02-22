#!/usr/bin/env python3

'''
This creates a multi-fasta file suitable for multiple sequence alignment from Blast result. Expects a blast result as single file XML2 format (must include taxid!) and needs internet connection to download sequences.\nThe TaxFinder module can be used to unify subspecies but it is not necessary.
'''

import sys
import os
import argparse
try:
	from taxfinder import TaxFinder
except ImportError:
	pass

import phylogenetics as phylo


def main():
	'''
	Main entry point.
	'''

	parser = argparse.ArgumentParser(description=__doc__)

	parser.add_argument('xml',
		help='The Blast result as single file XML2 format (outfmt 16).')

	parser.add_argument('-m', '--mail',
		help='Please state your (real!) email address. Alternatively, '
		'you can define the environment variable BLASTMAIL.')

	parser.add_argument('-o', '--outfile', default='',
		help='Outfile name. Leave empty to write to stdout.')

	parser.add_argument('-l', '--logfile', default='',
		help='Logfile name. Leave empty to write to stderr.')

	parser.add_argument('-s', '--strip', action='store_true',
		help='If given, stop codons are stripped off.')

	parser.add_argument('-e', '--evalue', type=float, default=1e-30,
		help='Evalue cutoff for including entries [1e-30]')

	parser.add_argument('-t', '--title', type=int, default=0,
		help='Shorten the title of the entries to this length. Default '
		'is 0 (no shortening).')

	parser.add_argument('-p', '--protein', action='store_true',
		help='If given, protein accessions instead of nucleotide '
		'accessions are assumed.')

	args = parser.parse_args()

	if args.logfile:
		logfile = open(args.logfile, 'w')
	else:
		logfile = sys.stderr

	try:
		mail = os.environ['BLASTMAIL']
	except KeyError:
		mail = args.mail
	if not mail:
		print('\033[1;31mPlease set your email address using -m or the '
		'environmental variable BLASTMAIL!\033[0;0m', file=logfile)
		sys.exit()

	try:
		TF = TaxFinder()
	except NameError:
		TF = None
		print('Taxfinder module not found. Script continues, but '
		'unifying subspecies will not work.', file=logfile)

	# Get the relevant entries based on evalue cutoff
	entries = phylo.get_entries_from_blast_result(args.xml, args.evalue, TF)

	# Download the sequences
	if args.protein:
		fasta = phylo.download_protein_sequences(entries, mail=mail, title=args.title)
	else:
		fasta = phylo.download_nucleotide_sequences(entries, mail=mail, title=args.title, strip=args.strip)

	# And save the sequences
	if args.outfile:
		open(args.outfile, 'w').write(fasta)
	else:
		print(fasta)

	if args.logfile:
		logfile.close()
