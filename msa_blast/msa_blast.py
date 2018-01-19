#!/usr/bin/env python3


### config:
mail = '' # Enter your email address here!
cutoff = 1e-30 # E-value cutoff


import sys
from Bio import Entrez, SeqIO
try:
	from taxfinder import TaxFinder
except ImportError:
	pass


if mail:
	mail_warning = ''
else:
	mail_warning = '\n\033[1;31mPlease change your email address in the script before running!\033[0;0m'

if len(sys.argv) < 2 or sys.argv[1].startswith('-'):
	print('''Makes multi-fasta file suitable for multiple sequence alignment. Expects a
blastresult file in XML format (the new xml format; must include taxid!) and
needs internet connection.{}
Dependecies: BioPython and optionally TaxFinder
Usage:
python3 {} blastresult.xml'''.format(mail_warning, __file__))
	sys.exit()

if not mail:
	print(mail_warning.strip())
	sys.exit()

try:
	TF = TaxFinder()
except NameError:
	TF = None
	print('Taxfinder module not found. Scripts continues, but unifying subspecies will not work.')

Entrez.email = mail

entries = {}

acc = None
taxid = None
sciname = None
evalue = None

with open(sys.argv[1]) as xml:
	for line in xml:
		if '<accession>' in line:
			acc = line.strip().replace('<accession>', '').replace('</accession>', '')

		elif not acc:
			continue

		elif '<taxid>' in line:
			taxid = int(line.strip().replace('<taxid>', '').replace('</taxid>', ''))

			if TF:
				taxid = TF.getSpeciesFromSubspecies(taxid)


		elif '<sciname>' in line:
			sciname = line.strip().replace('<sciname>', '').replace('</sciname>', '')

		elif '<evalue>' in line:
			evalue = float(line.strip().replace('<evalue>', '').replace('</evalue>', ''))

			if taxid not in entries or entries[taxid][1] > evalue:
				entries[taxid] = (acc, evalue, sciname)

			acc = None
			taxid = None
			sciname = None
			evalue = None

total = len(entries)
cdss = {}

with open('all_seqs.fasta', 'w') as outfile:
	for current, taxid in enumerate(entries, start=1):

		print(' '*75, end='\r', flush=True)

		print('Running seq {:3d} of {:3d}: {:<15}'.format(current, total, entries[taxid][0]), end='', flush=True)

		if entries[taxid][1] > cutoff:
			continue

		print('dl:', end='', flush=True)
		try:
			handle = Entrez.efetch(db='nuccore', id=entries[taxid][0], rettype='gb', retmode='text')
		except Exception as err:
			print('\r{}: {}'.format(entries[taxid][0], err))
			continue

		print('ok parse:', end='', flush=True)
		seqobj = SeqIO.parse(handle, 'gb')
		print('ok', end='\r', flush=True)
		record = next(seqobj)
		sequence = str(record._seq)

		featurecds = []
		sequences = {}
		taxids = {}

		for feature in record.features:
			if feature.type == 'CDS':
				name = feature.qualifiers['product'][0]
				start = int(feature.location._start)
				end = int(feature.location._end)
				featurecds.append((name, start, end))

		if len(featurecds) == 0:
			print('\r{}: No CDS{}'.format(entries[taxid][0], ' '*40))
			continue

		if len(featurecds) > 1:
			cdss[entries[taxid][0]] = tuple(featurecds)
			sequences[entries[taxid][0]] = sequence
			taxids[entries[taxid][0]] = taxid
			continue

		cds = sequence[start:end]
		if not cds.startswith('ATG') or not cds[-3:] in ['TAA', 'TGA', 'TAG']:
			print('\r{}: No ATG or Stop codon found!{}'.format(entries[taxid][0], ' '*30))
			continue

		outfile.write('>{}_{}\n{}\n'.format(entries[taxid][0], entries[taxid][2].replace(' ', '-'), cds))

print('')

if cdss:
	print('There were accessions with multiple cds:')
	for acc in cdss:
		print(acc)
		if len(cdss[acc]) > 20:
			print('More than 20 CDS found. Please investigate yourself.')
			continue

		for i, feature in cdss[acc]:
			print('{:2d}: {}'.format(i, feature[0]))

		try:
			chosen = int(input('Which one do you want (x for none)? '))
		except ValueError:
			print('Chose none.')
			continue

		print('Chose {}'.format(cdss[acc][chosen][0]))

		start = cdss[acc][chosen][1]
		end = cdss[acc][chosen][2]

		cds = sequences[acc][start:end]
		if not cds.startswith('ATG') or not cds[-3:] in ['TAA', 'TGA', 'TAG']:
			print('No ATG or Stop codon found!'.format(entries[taxid][0]))
			continue

		outfile.write('>{}_{}\n{}\n'.format(entries[taxids[acc]][0], entries[taxids[acc]][2].replace(' ', '-'), cds))
