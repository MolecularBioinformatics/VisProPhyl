import sys
#import taxfinder as tf
from PIL import Image
from argparse import ArgumentParser
from ete3 import Tree
import pandas as pd
import numpy as np
from scipy.stats.mstats import gmean

def getCategory(line):
	if line.startswith('GPT'):
		return (255,0,0)
	elif line.startswith('PPT'):
		return (0,255,0)
	elif line.startswith('TPT'):
		return (0,0,255)
	elif line.startswith('XPT'):
		return (255,255,0)
	elif line.startswith('algae'):
		return (0,255,255)

	return (0,0,0)


def toList(s):
	return s.split(',')


parser = ArgumentParser(description='Show multiple sequence alignments pixel-wise')
parser.add_argument('msa', help='The alignment file in fasta format')
parser.add_argument('-o', '--order', help='A file with the order in which the alignments should be shown. One element per line.')
parser.add_argument('-n', '--newick', help='A newick file with all elements from which the order can be derived and clusters shown.')
parser.add_argument('-e', '--exon', help='A file with exons. Only proteins with known exon borders will be shown.')
parser.add_argument('-t', '--thres', type=float, default=0, help='Distance threshold for when a new cluster shall be opened.')
parser.add_argument('-d', '--seeds', type=toList, default=[], help='Seeds for cluster determination, separated by comma')
parser.add_argument('-b', '--border', type=int, default=10, help='Border with in pixels.')
parser.add_argument('-a', '--nucacid', action='store_true', help='Is it Nucleic acid? (Otherwise, it is a protein sequence.)')
#colors should be default, replace by nocolor ? (thereby freeing up -c for clusters or catgeories)
parser.add_argument('-c', '--color', action='store_true', help='Plot colors?')
parser.add_argument('-q', '--quiet', action='store_true', help='Do not be verbose')
parser.add_argument('-s', '--save', default='msa.png', help='Set the output file name. File ending determines file type. Defaults to msa.png.')
parser.add_argument('--clusterout', nargs=2, help='Filename for annotation tsv-file with TaxID information & filename for output file.')
args = parser.parse_args()


alignmentFile = args.msa
orderFile = args.order
newickFile = args.newick
exonFile = args.exon
threshold = args.thres
seeds = args.seeds		# for all e.g. ['XPT_monocots_Phoenix-dactylifera_gi672166121', 'TPT_dicots_Eucalyptus-grandis_gi702348122', 'GPT_dicots_Camelina-sativa_gi727606022', 'PPT_dicots_Citrus-clementina_gi567908459']
printColor = args.color
dna = args.nucacid
quiet = args.quiet
saveFN = args.save
saveclusters = args.clusterout


borders = 10

if dna:
	aacolors = {
	'A': (203, 0, 0),
	'T': (0, 0, 203),
	'U': (0, 0, 203),
	'G': (0, 203, 0),
	'C': (203, 0, 191),
	'N': (0, 0, 0),
	'-': (255, 255, 255),
	#these are sometimes added by clustalo aswell
	'R': (0, 0, 0), #puRine, A/G
	'Y': (0, 0, 0), #pYrimidine, C/T
	'W': (0, 0, 0),
	'K': (0, 0, 0),
	'M': (0, 0, 0),
	'S': (0, 0, 0)

	}
else:
	aacolors = {
	'A': (203, 0, 0),
	'V': (203, 0, 0),
	'F': (203, 0, 0),
	'P': (203, 0, 0),
	'M': (203, 0, 0),
	'I': (203, 0, 0),
	'L': (203, 0, 0),
	'W': (203, 0, 0),
	'D': (0, 0, 203),
	'E': (0, 0, 203),
	'R': (203, 0, 191),
	'K': (203, 0, 191),
	'S': (0, 203, 0),
	'T': (0, 203, 0),
	'Y': (0, 203, 0),
	'H': (0, 203, 0),
	'C': (0, 203, 0),
	'N': (0, 203, 0),
	'G': (0, 203, 0),
	'Q': (0, 203, 0),
	'X': (0, 0, 0),			# Unknown
	'B': (170, 170, 170),	# Asx (N or D)
	'Z': (170, 170, 170),	# Glx (Q or E)
	'J': (203, 0, 0),		# Xle (I or L)
	'-': (255, 255, 255)	# Gap
	}

# Easily distinguishable colors from colorbrewer2.org, some way to figure out which one to use ?
# 7
clustercolors = [(228,26,28),(55,126,184),(77,175,74),(152,78,163),(255,127,0),(255,255,51),(166,86,40)]
# 14 (12 + black + grey/white)
clustercolors = [(166,206,227),(31,120,180),(178,223,138),(51,160,44),(251,154,153),(227,26,28),(253,191,111),(255,127,0),(202,178,214),(106,61,154),(255,255,153),(177,89,40),(0,0,0),(224,224,224)]

#print('initializing TaxFinder... ')
#
#TF = tf.TaxFinder()

if not quiet and (orderFile or newickFile or exonFile):
	print('reading order and categorizing... ')

order = []
category = {}
clusters = {}
## added for clusterout(), should probably be a class.variable
clusterlist = []
##
exons = {}
orderDone = False

if exonFile:
	ecolors = [(200,147,45), (99,93,220)]
	with open(exonFile, 'r') as f:
		for line in f:
			if line.startswith('#'):
				continue
			acc, e = line.rstrip().split('\t')
			# Turns the string e.g. ((1, 4), (4, 8)) to an actual tuple of tuples (basically the same as eval(e) would do)
			e = tuple(tuple(int(y) for y in x.split(', ')) for x in e[1:-1].split('), ('))
			ex = [(0,0,0) for _ in range(e[-1][1])]
			num = 0
			for exon in e:
				num = not num
				for pos in range(exon[0] - 1, exon[1]):
					if ex[pos] == (0,0,0):
						ex[pos] = ecolors[num]
					else:
						ex[pos] = (28,208,28)
			exons[acc] = tuple(ex)

if orderFile:
	with open(orderFile) as f:
		for line in f:
			line = line.rstrip()
			acc = '_'.join(line.split('_')[3:])
			if not exons or acc in exons:
				order.append(line)
			category[line] = getCategory(line)
	#?? orderDone = True ??

# if '.fa' in alignmentFile:
# 	mname = alignmentFile.split('.fa')[0] + '_matrix.tsv'
# else:
# 	mname = alignmentFile.split('.')[0] + '_matrix.tsv'
mname = 'clk_promoters_aligned_matrix.tsv'

def treeClusters(t, threshold, clustercolors, matrixfile=mname):

	#These should be specified by a config file or flag or whatever
	meanfunc1 = gmean
	meanfunc2 = gmean

	#t = Tree(newickFile)

	# matrix conatins also the dropped elements, even though probably/not yet not needed
	pws = pd.read_csv(matrixfile, sep='\t', header=1, index_col=0, na_values='-')

	## !discarded (would be needed if no data is written on nodes & instead always pulled from matrix)
	# Drop upper half of (symmetric!) matrix to get rid of duplicates
	# pws.values[np.triu_indices_from(pws)] = np.nan
	##

	for node in t.traverse('postorder'):
		# Only leaves are named!

		# maybe needs to be 1 (to allow single sequences to become a cluster
		node.add_features(sim=None)

		if node.is_leaf():
			# names in tree are 'promoter|acc'
			# names in matix are only acc
			# This should be changed at some poin
			#node.name = node.name.split('|')[1]

			node.add_features(cl=None)

			node.add_features(distances = pws[node.name.split('|')[1]])

		else:
			# differentiate between children that are leaves (get mean of values for all other current leaves from pw similarity matrix)
			# and children that are not (already have values)
			leaves = node.get_leaf_names()

			values = []

			for n in node.children:
				if n.name in leaves:
					meansim = meanfunc1(n.distances[[n.split('|')[1] for n in leaves]].dropna())
					#print(n.distances[leaves].dropna())
					if not np.isnan(meansim):
						values.append(meansim)
				else:
					values.append(n.sim)

			#print(values)
			node.sim = meanfunc2(values)

	# Cluster nodes based on min/max similarity threshhold
	cl = []
	for node in t.traverse('postorder'):
		if node.is_leaf():
			continue

		if node.sim < threshold or node.is_root():
			for x in node.iter_descendants():
				if any(l.cl for l in x.iter_leaves()):
					continue

				# 'promoter|' + node.name
				cl.append(x.get_leaf_names())
				for leaf in x.get_leaves():
					leaf.cl = len(cl)-1

	cleanclusters = []
	current = []
	for cluster in cl:
		if len(cluster) <= 2:
			current += cluster
		else:
			if current:
				cleanclusters.append(current)
				current = []
			cleanclusters.append(cluster)

	global clusterlist
	clusterlist = cleanclusters


	clusters = {}
	for i, c in enumerate(cleanclusters):
		clusterNo = i % len(clustercolors)
		for name in c:
			clusters[name] = clustercolors[clusterNo]

	return clusters



def getClusters(t, threshold, clustercolors):
	clusters = {}
	currentCluster = None
	clusterNo = 0

	for node in t.traverse('postorder'):
		if node.is_leaf():
			if currentCluster == None:
				currentCluster = node
				clusterlist.append(list())

			if t.get_distance(currentCluster, node) > threshold:
				currentCluster = node
				clusterlist.append(list())

				clusterNo = int((clusterNo + 1) % len(clustercolors))

			clusterlist[-1].append(node.name)
			clusters[node.name] = clustercolors[clusterNo]

	return clusters


def getFixedClusters(t, seedNames, clustercolors):
	seeds = []
	for node in t.traverse('postorder'):
		if node.is_leaf():
			for seedName in seedNames:
				if seedName in node.name:
					seeds.append(node)
					break

	#seedNames = ['XPT_monocots_Phoenix-dactylifera_gi672166121',
	#'TPT_dicots_Eucalyptus-grandis_gi702348122',
	#'GPT_dicots_Camelina-sativa_gi727606022',
	#'PPT_dicots_Citrus-clementina_gi567908459']
	#'GPT_lower_Klebsormidium-flaccidum_gi971511185']

	#seeds = [t.search_nodes(name = name)[0] for name in seedNames]

	clusters = {}

	#added for clusterout()
	#This function is 'default' so at least one cluster needs to be opened, even without seeds
	global clusterlist
	if seeds:
		clusterlist += [list() for _ in seeds]
	else:
		clusterlist = [list()]

	for node in t.traverse('postorder'):
		if node.is_leaf():

			dist = 1e6
			idx = 0

			for i, seed in enumerate(seeds):
				d = t.get_distance(seed, node)
				if d < dist:
					dist = d
					idx = i

			clusterlist[idx].append(node.name)
			clusters[node.name] = clustercolors[idx]

	return clusters


if newickFile:
	mytree = Tree(newickFile)

	if threshold:
		#clusters = getClusters(mytree, threshold, clustercolors)
		clusters = treeClusters(mytree, threshold, clustercolors)
	else:
		clusters = getFixedClusters(mytree, seeds, clustercolors)

	if not orderDone:
		for node in mytree.traverse('postorder'):
			if node.is_leaf():
				nname = node.name
				#acc format hard-coded
				acc = '_'.join(nname.split('_')[3:])
				if not exons or acc in exons:
					order.append(nname)
				category[nname] = getCategory(nname)
		orderDone = True


if not quiet:
	print('reading alignment... ')

data = {}
current = None
pos = 0

with open(alignmentFile) as f:
	for line in f:
		line = line.rstrip()

		if line.startswith('>'):
			pos = 0
			current = line.split()[0][1:]
			data[current] = []
			#Hard coded for exon stuff ?!
			acc = '_'.join(current.split('_')[3:])
			if not orderDone and (not exons or acc in exons):
				order.append(current)
				category[current] = getCategory(current)
		else:
			if exons and acc not in exons:
				del(data[current])
				continue
			for char in line:
				if exons:
					if char == '-':
						data[current].append((255,255,255))
					else:
						try:
							data[current].append(exons[acc][pos])
						except IndexError:
							data[current].append((0,0,0))
						pos += 1
				elif printColor:
					data[current].append(aacolors[char])
				else:
					if char == '-':
						data[current].append((255,255,255))
					else:
						data[current].append((0,0,0))


def clusterout(anFN, outFN):

	mapping = {}

	with open(anFN, 'r') as f:
		#The elements of order *should* be the fasta-headers/Newick names (without '>')
		elements = list(map(lambda x: x.split('|')[1], order))
		next(f)
		for l in f:
			lline = l.rstrip().split('\t')
			if lline[0] in elements:
				mapping[lline[0]] = lline[1]


	with open(outFN, 'w') as out:
		out.write('# Names/Acc-IDs/... of MSA sequences sorted by clusters. Each lines containes comma separated entries for one cluster.\n')
		out.write("# Order of clusters is defined by Newick file '{}', unless the -o (--order) option was given to msa_viewer this also corresponds to the top-down order of the result picture.\n".format(newickFile))
		for cluster in clusterlist:
			out.write(','.join(list(map(lambda x: x.split('|')[1]+'^'+mapping[x.split('|')[1]], cluster)))+'\n')


if saveclusters:
	clusterout(*saveclusters)


height = len(order)
width = len(next(iter(data.values())))

if not quiet:
	print('creating image... ')

im = Image.new('RGB', (width + 4*borders, height), (255,255,255))

for y, elem in enumerate(order):
	for x in range(borders):
		#should be optional
		im.putpixel((x + borders, y), category[elem])
		im.putpixel((x + width + 2*borders, y), category[elem])

		if clusters:
			im.putpixel((x, y), clusters[elem])
			im.putpixel((x + width + 3*borders, y), clusters[elem])
	try:
		for x, value in enumerate(data[elem]):
			im.putpixel((x + 2*borders, y), value)
	except KeyError:
		print(elem, alignmentFile)
		raise

im.save(saveFN)
