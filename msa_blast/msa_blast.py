#!/usr/bin/env python3

mail = 'mathias.bockwoldt@uit.no' # Enter your email address here!

import sys
import os
import re
from Bio import Entrez, SeqIO
try:
	from taxfinder import TaxFinder
except ImportError:
	pass


logfile = sys.stderr


def get_entries(filename, cutoff, TF=None):

	entries = {}

	acc = None
	taxid = None
	sciname = None
	evalue = None
	hitfrom = None
	hitto = None

	with open(filename) as xml:
		for line in xml:
			if '<accession>' in line:
				acc = line.strip().replace('<accession>', '').replace('</accession>', '')

			elif not acc:
				continue

			elif '<taxid>' in line:
				taxid = int(line.strip().replace('<taxid>', '').replace('</taxid>', ''))

				try:
					taxid = TF.getSpeciesFromSubspecies(taxid)
				except ValueError:
					print('Taxid {} not found!'.format(taxid), file=sys.stderr)
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

	return entries


def dl_sequences(entries, strip, title):

	out = []
	total = len(entries)
	cdss = {}

	for current, taxid in enumerate(entries, start=1):
		print('\r' + ' '*75, end='\r', flush=True, file=sys.stderr)

		print('Running seq {:3d} of {:3d}: {:<15}'.format(current, total, entries[taxid][0]), end='', flush=True, file=sys.stderr)

		print('dl:', end='', flush=True, file=sys.stderr)
		try:
			handle = Entrez.efetch(db='nuccore', id=entries[taxid][0], rettype='gb', retmode='text')
		except Exception as err:
			print('\r{}: {}'.format(entries[taxid][0], err), file=sys.stderr)
			continue

		print('ok parse:', end='', flush=True, file=sys.stderr)
		seqobj = SeqIO.parse(handle, 'gb')
		record = next(seqobj)
		sequence = record._seq

		featurecds = None

		for feature in record.features:
			if feature.type == 'CDS':
				try:
					name = feature.qualifiers['product'][0]
					start = int(feature.location._start)
					end = int(feature.location._end)
				except (AttributeError, KeyError):
					continue

				if start <= entries[taxid][3] and end >= entries[taxid][4]:
					featurecds = (name, start, end)
					break
		else:
			print('\r{}: No CDS{}'.format(entries[taxid][0], ' '*40), file=logfile)
			continue

		cds = str(sequence[start:end])

		if not cds[:3] == 'ATG' or not cds[-3:] in ['TAA', 'TGA', 'TAG']:
			if cds[-3:] == 'CAT' and cds[:3] in ['CTA', 'TCA', 'TTA']:
				cds = str(sequence[start:end].reverse_complement())
			else:
				print('\r{}: No ATG or Stop codon found! Sequence will be omitted{}'.format(entries[taxid][0], ' '*30), file=logfile)
				continue

		if len(cds) % 3:
			print('\r{}: Possible frameshit! Sequence will be omitted{}'.format(entries[taxid][0], ' '*40), file=logfile)
			continue

		if re.search(r'[^ACGT]', cds):
			print('\r{}: Non canonical DNA base! Sequence will be included in output.{}'.format(entries[taxid][0], ' '*40), file=logfile)

		if strip and cds[-3:] in ['TAA', 'TGA', 'TAG']:
			cds = cds[:-3]

		fasta_head = '>{}_{}'.format(entries[taxid][0], entries[taxid][2].replace(' ', '-'))
		if title:
			fasta_head = fasta_head[:title]

		out.append('{}\n{}'.format(fasta_head, cds))

	print('', file=logfile)

	return '\n'.join(out)


def dl_protein_seqs(entries, title):

	out = []
	total = len(entries)

	for current, taxid in enumerate(entries, start=1):
		print('\r' + ' '*75, end='\r', flush=True, file=sys.stderr)

		print('Running seq {:3d} of {:3d}: {:<15}'.format(current, total, entries[taxid][0]), end='', flush=True, file=sys.stderr)

		print('dl:', end='', flush=True, file=sys.stderr)
		try:
			handle = Entrez.efetch(db='protein', id=entries[taxid][0], rettype='fasta', retmode='text')
		except Exception as err:
			print('\r{}: {}'.format(entries[taxid][0], err), file=sys.stderr)
			continue

		fasta = str(SeqIO.read(handle, 'fasta')._seq)
		handle.close()

		print('ok parse:', end='', flush=True, file=sys.stderr)

		fasta_head = '>{}_{}'.format(entries[taxid][0], entries[taxid][2].replace(' ', '-'))
		if title:
			fasta_head = fasta_head[:title]

		out.append('{}\n{}'.format(fasta_head, fasta))

	print('', file=logfile)

	return '\n'.join(out)


if __name__ == '__main__':
	import argparse

	parser = argparse.ArgumentParser(description='Makes multi-fasta file suitable for multiple sequence alignment from Blast result. Expects a blast result as single file XML2 format (must include taxid!) and needs internet connection.\nThe TaxFinder module might be very usefull although it is not necessary.')

	parser.add_argument('xml', help='The Blast result as single file XML2 format (outfmt 16).')
	parser.add_argument('-m', '--mail', help='Please state your (real!) email address. Alternatively, you can hard-code it in the script on line 3 or define the environment variable BLASTMAIL.')
	parser.add_argument('-o', '--outfile', default='', help='Outfile name. Leave empty to write to stdout.')
	parser.add_argument('-l', '--logfile', default='', help='Logfile name. Leave empty to write to stderr.')
	parser.add_argument('-s', '--strip', action='store_true', help='If given, stop codons are stripped off.')
	parser.add_argument('-e', '--evalue', type=float, default=1e-30, help='Evalue cutoff for including entries')
	parser.add_argument('-t', '--title', type=int, default=0, help='Shorten the title of the entries to this length. Default is 0 (no shortening).')
	parser.add_argument('-p', '--protein', action='store_true', help='If given, protein accessions are assumed.')

	args = parser.parse_args()

	if args.logfile:
		logfile = open(args.logfile, 'w')
	else:
		logfile = sys.stderr

	if 'BLASTMAIL' in os.environ:
		mail = os.environ['BLASTMAIL']
	if args.mail:
		mail = args.mail
	if not mail:
		print('\033[1;31mPlease change your email address in the script before running or set it via -m or the environmental variable BLASTMAIL!\033[0;0m', file=logfile)
		sys.exit()

	try:
		TF = None#TaxFinder()
	except NameError:
		TF = None
		print('Taxfinder module not found. Script continues, but unifying subspecies will not work.', file=logfile)

	Entrez.email = mail

	entries = get_entries(args.xml, args.evalue, TF)
	if args.protein:
		fasta = dl_protein_seqs(entries, args.title)
	else:
		fasta = dl_sequences(entries, args.strip, args.title)

	if args.outfile:
		open(args.outfile, 'w').write(fasta)
	else:
		print(fasta)

	if args.logfile:
		logfile.close()
