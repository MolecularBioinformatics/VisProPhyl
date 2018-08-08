#!/usr/bin/env python3

import sys

def run_blastp(query, outfilename, db, evalue = 1, maxthreads = cpu_count()):
	'''
	Run Blastp. This may take a long time (several minutes up to half-an-hour depending on the computer power and the size of the query and the database.

	:param query: The filename of the query sequence (may include a path)
	:param outfilename: The filename (including path) to save the result to
	:param db: The database file
	:param evalue: The e-value cut-off (should usually be very high, as the results will be filtered out afterwards)
	:param maxthreads: The maximum number of threads to use by Blast
	:creates: `outfilename`
	'''

	stdout, stderr = NcbiblastpCommandline(query = query, db = db, evalue = evalue, outfmt = 5, out = outfilename, num_threads = maxthreads, max_target_seqs = 20000)

	if stdout:
		print(stdout, file=sys.stderr)
	if stderr:
		print(stderr, file=sys.stderr)


def parse_blast_result(blastXML, TF, top = 0, exclude=None, newHeader=True):
	'''
	Parses Blast result XML files and writes the best or all results with less information in a tsv file.

	:param blastXML: The filename of the Blast output (must be output type 5)
	:param TF: An instance of the TaxFinder class
	:param top: Write only the best `top` hits to file. If `top` is 0, all hits are saved.
	:param contaminants: Set with taxids of species to exclude from the results
	:param newHeader: Where the Blast results produced with new headers (database from 2016 and newer)?
	:returns: csv table as string with the results
	'''

	if contaminants is None:
		contaminants = set()

	if top < 0:
		top = 0

	result = []

	with open(blastXML, 'r') as f, open(outfile, 'w') as out:
		records = NCBIXML.parse(f)

		result.append('\t'.join(('Tax-ID', 'Acc', 'Species', 'Rank', 'e-value', 'Length', 'Lineage', 'Prot-Name', 'Query-Protein')))

		for record in records:
			for i, alignment in enumerate(record.alignments):
				if top and i > top:
					break

				infos = TF.getInfoFromHitDef(alignment.hit_id, alignment.hit_def, newHeader = newHeader)

				for info in infos:
					if info['taxid'] in contaminants:
						continue

					lineage = ', '.join(TF.getLineage(info['taxid'], display = 'name'))

					for hsp in alignment.hsps:
						try:
							line = '\t'.join((str(info['taxid']), info['acc'], info['name'], info['rank'], str(hsp.expect), str(hsp.align_length), lineage, info['protname'], record.query.split('|')[1]))
						except IndexError:
							line = '\t'.join((str(info['taxid']), info['acc'], info['name'], info['rank'], str(hsp.expect), str(hsp.align_length), lineage, info['protname'], record.query))

						result.append(line)

	return '\n'.join(result)
