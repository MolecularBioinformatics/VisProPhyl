#!/usr/bin/env python3

from taxfinder import TaxFinder
from Bio.Blast import NCBIXML
import sys
import os
import re
import math
import numpy as np
import subprocess
import matplotlib.pyplot as plt
from matplotlib import gridspec
import pickle
from Bio import Entrez, SeqIO


Entrez.email = 'mathias.bockwoldt@uit.no'


def dl_protein_seq(acc):

	try:
		return open('cache/{}.seq'.format(acc)).read()
	except IOError:
		pass

	try:
		handle = Entrez.efetch(db='protein', id=acc, rettype='fasta', retmode='text')
	except Exception as err:
		print('{}: {}'.format(acc, err), file=sys.stderr)
		raise

	seq = str(SeqIO.read(handle, 'fasta')._seq)
	handle.close()

	open('cache/{}.seq'.format(acc), 'w').write(seq)

	return seq


nodedict = {}


def get_node_from_taxid(loi, taxid, TF):
	if taxid in nodedict:
		return nodedict[taxid]

	lineage = TF.get_lineage_fast(taxid)

	if lineage == loi:
		node = loi[-1]
	else:
		node = None
		for i, taxon in enumerate(lineage):
			if taxon not in loi:
				if i > 0:
					node = lineage[i-1]
				break

		if node == loi[0]:
			node = None

	nodedict[taxid] = node

	return node


def get_evalues(loi, blastXMLlist, TF):
	loidict = {x: None for x in loi}

	for blastXML in blastXMLlist:
		with open(blastXML, 'r') as f:
			records = NCBIXML.parse(f)

			for record in records:
				for i, descr in enumerate(record.descriptions):
					for hit in descr.items:
						taxid = hit.taxid

						node = get_node_from_taxid(loi, taxid, TF)
						if node is None:
							continue

						for hsp in record.alignments[i].hsps:
							try:
								elog = -1 * math.log10(hsp.expect)
							except ValueError:
								elog = 200 # if e == 0.0

							if loidict[node] is None or elog > loidict[node]:
								loidict[node] = elog

	return loidict


def get_best_hit_ids(loi, blastXMLlist, TF):
	loidict = {x: {} for x in loi}

	for blastXML in blastXMLlist:
		with open(blastXML, 'r') as f:
			records = NCBIXML.parse(f)

			for record in records:
				for i, descr in enumerate(record.descriptions):
					for hit in descr.items:
						taxid = hit.taxid
						acc = hit.accession
						node = get_node_from_taxid(loi, taxid, TF)
						if node is None:
							continue

						for hsp in record.alignments[i].hsps:
							try:
								elog = -1 * math.log10(hsp.expect)
							except ValueError:
								elog = 200 # if e == 0.0

							if taxid not in loidict[node]:
								loidict[node][taxid] = [(elog, acc)]
							elif elog > loidict[node][taxid][0][0]:
								loidict[node][0] = (elog, acc)

	return loidict


def run_needle(acc1, acc2):
	subprocess.run(['needle', '-asequence', 'cache/{}.seq'.format(acc1), '-bsequence', 'cache/{}.seq'.format(acc2), '-gapopen', '10.0', '-gapextend', '0.5', '-outfile', 'tmp/aln'], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

	res = open('tmp/aln').read()
	sim = float(re.search(r'Similarity:[ ]+[0-9]+/[0-9]+[ ]+\(([0-9\.]+)%\)', res).group(1))
	#sim = float(re.search(r'Identity:[ ]+[0-9]+/[0-9]+[ ]+\(([0-9\.]+)%\)', res).group(1))

	return sim


similarity_cache = {}

def get_similarity(acc1, acc2):
	if acc1 == acc2:
		return 100

	if acc1 > acc2:
		acc1, acc2 = acc2, acc1

	if (acc1, acc2) in similarity_cache:
		return similarity_cache[(acc1, acc2)]

	try:
		seq1 = dl_protein_seq(acc1)
		seq2 = dl_protein_seq(acc2)
	except Exception:
		raise
		return 0

	similarity = run_needle(acc1, acc2)

	similarity_cache[(acc1, acc2)] = similarity

	return similarity



def get_simvalues(loi, blastXMLlist, ref_hits, TF):
	values = {x: None for x in loi}
	best_hits = get_best_hit_ids(loi, blastXMLlist, TF)

	for node in best_hits:
		if len(best_hits[node]) < 2:
			continue
		for taxid in best_hits[node]:
			if taxid in ref_hits[node]:
				try:
					value = get_similarity(best_hits[node][taxid][0][1], ref_hits[node][taxid][0][1])
					if not values[node] or value > values[node]:
						values[node] = value
					if value == 100:
						break
				except IndexError:
					print(best_hits[node], ref_hits[node])
					raise

	return values


def get_loi(target1, target2):
	lineage1 = TF.get_lineage_fast(int(target1))
	lineage2 = TF.get_lineage_fast(int(target2))

	raise NotImplementedError


if __name__ == '__main__':
	from glob import glob

	TF = TaxFinder()

	panels = []
	interesting = {}

	with open(sys.argv[1]) as f:
		loi = next(f).split()
		if len(loi) == 2:
			lineage_of_interest = get_loi(*loi)
		else:
			lineage_of_interest = list(TF.get_lineage_fast(int(loi[0])))

		for line in f:
			line = line.rstrip()
			if not line:
				panels.append([])
				continue

			lline = line.split()
			fns = []
			for fn in lline[1].split(','):
				fns.extend(glob(fn))
			panels[-1].append(lline[0])
			interesting[lline[0]] = tuple(fns)

	picklename = sys.argv[1] + '.p'

	#try:
	#	values = pickle.load(open(picklename, 'rb'))
	#except IOError:

	values = {}

	for panel in panels:

		reference_hits = get_best_hit_ids(lineage_of_interest, interesting[panel[0]], TF)

		for key in panel[1:]:
			values[key] = get_simvalues(lineage_of_interest, interesting[key], reference_hits, TF)

		values_buffer = {}
		for taxid in lineage_of_interest:
			value = None
			for name in values:
				if values[name][taxid] is not None:
					value = 100
					break
			values_buffer[taxid] = value

		values[panel[0]] = values_buffer # {x: 100 for x in lineage_of_interest}


	#	pickle.dump(values, open(picklename, 'wb'))


	#for i in range(len(lineage_of_interest)):
	#	print('{: <20}: {}'.format(lineage_names[i], values[lineage_of_interest[i]]))

	to_delete = {x: True for x in lineage_of_interest}

	for name in values:
		for taxid in values[name]:
			if values[name][taxid] is not None:
				to_delete[taxid] = False

	for taxid in to_delete:
		if to_delete[taxid]:
			lineage_of_interest.remove(taxid)

	# Can't use the simple get_lineage, because some parts were probably deleted a few lines up
	lineage_names = [TF.get_name_from_id(taxid) for taxid in lineage_of_interest]

	img_width = 16
	img_height = 8*len(panels)
	fig = plt.figure(figsize = (img_width, img_height))
	gs = gridspec.GridSpec(len(panels), 2, width_ratios=[3, 1])
	axlist = []
	for num, panel in enumerate(panels):
		axlist.append(plt.subplot(gs[(num-1)*2]))
		ax = axlist[-1]

		colors = plt.cm.rainbow(np.linspace(0,1,len(panel)))
		markers = ['.', 'o', 'v', '^', '<', '>', '1', '2', '3', '4', '8', 's', 'p', '*', 'h', 'H', '+', 'x', 'D', 'd']

		x_values = np.arange(0, len(lineage_of_interest))
		legendlines = []
		legendnames = []

		for i, name in enumerate(sorted(panel)):
			value_list = np.array([values[name][x] for x in lineage_of_interest], dtype=np.float)
			pl = ax.plot(x_values, value_list, linestyle='-', color=colors[i], marker=markers[i%len(panel)])[0]

			legendlines.append(pl)
			legendnames.append(name)

			mask = np.isfinite(value_list)
			ax.plot(x_values[mask], value_list[mask], linestyle='--', color=colors[i])

		#ax.set_ylabel('-log(e-value)')
		#ax.set_ylim((0, 210))

		plt.sca(ax)
		plt.xticks(range(len(lineage_names)), lineage_names, rotation=90)

		axlist.append(plt.subplot(gs[(num-1)*2+1]))
		ax = axlist[-1]
		ax.axis('off')
		ax.legend(legendlines, legendnames, loc='center left')

	plt.tight_layout()

	if len(sys.argv) > 2 and sys.argv[2].startswith('sh'):
		plt.show()
	else:
		fn = sys.argv[1] + '.pdf'
		plt.savefig(fn)
		print('File saved to "{}"'.format(fn))


