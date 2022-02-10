#!/usr/bin/env python3

from taxfinder import TaxFinder
from Bio.Blast import NCBIXML
import sys
import os
import math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec
import pickle


nodedict = {}


def get_node_from_taxid(loi, taxid, TF):
	if taxid in nodedict:
		return nodedict[taxid]

	lineage = TF.getLineageFast(taxid)

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


def get_values(loi, blastXMLlist, TF):
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
			lineage_of_interest = list(TF.getLineageFast(int(loi[0])))

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

	try:
		values = pickle.load(open(picklename, 'rb'))
	except IOError:
		values = {}
		for key in interesting:
			values[key] = get_values(lineage_of_interest, interesting[key], TF)

		pickle.dump(values, open(picklename, 'wb'))

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

	# Can't use the simple getLineage, because some parts were probably deleted a few lines up
	lineage_names = [TF.getNameFromID(taxid) for taxid in lineage_of_interest]

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

		ax.set_ylabel('-log(e-value)')
		ax.set_ylim((0, 210))

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


