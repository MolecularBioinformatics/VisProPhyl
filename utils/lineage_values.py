#!/usr/bin/env python3


'''
This script is meant to be changed to fit your needs.
It will write a plot to "input_filename.pdf".
This plot will contain the best e-value for a seed that is not in the
remaining lineage.
Run like this:
$ python3 lineage_values.py input_filename
`input_filename` should point to a file that looks like this:

```
taxid

prot1	path/to/blastresult_prot1a.xml path/to/blastresult_prot1b.xml
prot2	path/to/blastresult_prot2.xml
...
```
'''


import math
import pickle
import sys

from glob import glob

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec
from Bio.Blast import NCBIXML

from taxfinder import TaxFinder

nodedict = {}


def get_node_from_taxid(loi, taxid, TF):
	'''
	Returns the best node for a given taxid in the lineage of interest (loi).
	'''

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


def get_values(loi, blast_xml_list, TF):
	'''
	Get evalues for the lineage of interest (loi) in each blast XML file.
	'''

	loidict = {x: None for x in loi}

	for blast_xml in blast_xml_list:
		with open(blast_xml, 'r') as infile:
			records = NCBIXML.parse(infile)

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


def main():
	'''
	Main entry point.
	'''

	TF = TaxFinder()

	panels = []
	interesting = {}

	with open(sys.argv[1]) as infile:
		loi = int(next(infile).strip())
		lineage_of_interest = list(TF.get_lineage_fast(loi))

		for line in infile:
			line = line.strip()
			if not line:
				panels.append([])
				continue

			lline = line.split()
			filenames = []
			for filename in lline[1:]:
				filenames.extend(glob(filename))
			panels[-1].append(lline[0])
			interesting[lline[0]] = tuple(filenames)

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

	# Can't use the simple get_lineage, because some parts were probably deleted a few lines up
	lineage_names = [TF.get_name_from_id(taxid) for taxid in lineage_of_interest]

	img_width = 16
	img_height = 8*len(panels)
	plt.figure(figsize = (img_width, img_height))
	grid = gridspec.GridSpec(len(panels), 2, width_ratios=[3, 1])
	axlist = []
	for num, panel in enumerate(panels):
		axlist.append(plt.subplot(grid[(num-1)*2]))
		axis = axlist[-1]

		colors = plt.cm.get_cmap('rainbow', np.linspace(0,1,len(panel)))
		markers = ['.', 'o', 'v', '^', '<', '>', '1', '2', '3', '4',
		'8', 's', 'p', '*', 'h', 'H', '+', 'x', 'D', 'd']

		x_values = np.arange(0, len(lineage_of_interest))
		legendlines = []
		legendnames = []

		for i, name in enumerate(sorted(panel)):
			value_list = np.array([values[name][x] for x in lineage_of_interest], dtype=float)
			plot = axis.plot(x_values, value_list, linestyle='-', color=colors[i], marker=markers[i%len(panel)])[0]

			legendlines.append(plot)
			legendnames.append(name)

			mask = np.isfinite(value_list)
			axis.plot(x_values[mask], value_list[mask], linestyle='--', color=colors[i])

		axis.set_ylabel('-log(e-value)')
		axis.set_ylim((0, 210))

		plt.sca(axis)
		plt.xticks(range(len(lineage_names)), lineage_names, rotation=90)

		axlist.append(plt.subplot(grid[(num-1)*2+1]))
		axis = axlist[-1]
		axis.axis('off')
		axis.legend(legendlines, legendnames, loc='center left')

	plt.tight_layout()

	if len(sys.argv) > 2 and sys.argv[2].startswith('sh'):
		plt.show()
	else:
		outfile = sys.argv[1] + '.pdf'
		plt.savefig(outfile)
		print('File saved to "{}"'.format(outfile))



if __name__ == '__main__':
	main()
