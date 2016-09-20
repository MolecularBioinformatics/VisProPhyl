#!/usr/bin/env python3

from ete3 import Tree, TreeStyle, NodeStyle, faces, PieChartFace, TextFace, CircleFace
import sys
import math
from operator import add
from matplotlib import colors as mpl_colors
from matplotlib import cm as mpl_colormap
from datetime import datetime
from taxfinder import TaxFinder

class TreeMaker(object):

	def __init__(self, tree, treefile, combo, prune, tf, bins=(30, 150, 5), sqrt = False, collapse = True, show_empty_species = True, startnode = None):

		if startnode is None:
			self.t = tree
		else:
			node = tree.search_nodes(name=startnode)
			if not node:
				raise KeyError('Node {} not found in tree.'.format(startnode))
			self.t = node[0]

		if combo.startswith('crosshits/'):
			combo = combo
		else:
			combo = 'crosshits/' + combo
		if combo.endswith('.tsv'):
			self.combo = combo
		else:
			self.combo = combo + '.tsv'


		self.sqrt = sqrt
		self.collapse = collapse
		self.empty = show_empty_species
		self.treefile = treefile
		self.prune = prune
		self.startnode = startnode

		# Read tree attributes into dict
		self.elements = {}
		with open(self.combo, 'r') as f:
			next(f)
			for line in f:
				ev, num, accs = line.rstrip('\n').split('\t')
				for acc in accs.split(','):
					tax = tf.getTaxID(acc)
					self.elements[tax] = int(ev)

		# self.featurelist: strs with feature names to be used (f_0, f_1, ...)
		self.featurelist = ['f_{}'.format(i) for i in range(bins[2])]

		# self.f_0, ... equals range of corresponding evalues (exp only)
		steps = int((bins[1] - bins[0]) / bins[2])
		for i, f in enumerate(self.featurelist):
			setattr(self, f, set(range(bins[0] + steps*i, bins[0] + steps*(i+1))))
		# extend last range up to bins[1]
		lastf = getattr(self, self.featurelist[-1])
		lastf.update(range(bins[1]-steps, bins[1]+1))


		#min(getattr(self, f)), max(getattr(self, f))
		self.featurenames = ['e-values between 1e-{} and 1e-{}'.format(min(getattr(self, f)), max(getattr(self, f))) for f in self.featurelist]


		# Colors for the features
		#                          0          1          2          3          4          5          6          7
		self.standardcolors = ['#AA0000', '#00AA00', '#55AAFF', '#DDDD00', '#AA00AA', '#3333FF', '#000000', '#c8c8c8']

		self.colors = self.standardcolors

		self.addFeatures()

		self.readPruningFile()

		self.pruneTree()

		self.treeLayout()

	# end __init__()


	def layout(self, node):
		try:
			percents = [round(100.0*getattr(node, f)/node.total) for f in self.featurelist] + [round(100.0*node.f_none/node.total)]
		except ZeroDivisionError:
			txt = TextFace(' ' + node.plainName, ftype='Arial', fgcolor='#000000')
		else:
			i=0
			if sum(percents) != 100:
				x = sorted(percents)[-1]
				i = percents.index(x)
				percents[i] += 100 - sum(percents)

			if self.sqrt:
				size = int(math.sqrt(node.total)*5+1)
			else:
				size = int(math.log10(node.total)*20+10)

			P = PieChartFace(percents, size, size, self.colors)
			faces.add_face_to_node(P, node, 0, position="branch-top")

			#fgcol = mpl_colors.rgb2hex(mpl_colormap.brg(col/10.0)) # foreground/text color depending on NAD consumer density
			fgcol = '#000000'

			# if self.count -> print node.total

			# if node.counter:
			# 	cnt = float(sum(node.counter))/float(node.total)
			# 	txt = TextFace(' {}\n {:.2f}'.format(node.plainName, cnt), ftype='Arial', fgcolor=fgcol)
			# else:
			txt = TextFace(' ' + node.plainName, ftype='Arial', fgcolor=fgcol)
			# #txt.margin_right = -30

		if node.is_leaf():
			faces.add_face_to_node(txt, node, 0, position="branch-right")
		else:
			faces.add_face_to_node(txt, node, 0, position="branch-bottom")

	# end layout


	def treeLayout(self):
		# Add a tree style and a legend

		self.ts = TreeStyle()
		self.ts.show_scale = False
		self.ts.show_leaf_name = False
		self.ts.extra_branch_line_type = 0
		self.ts.extra_branch_line_color = 'black'
		self.ts.layout_fn = self.layout

		self.ts.legend_position = 4

		for name, color in zip(self.featurenames, self.colors):
			self.ts.legend.add_face(CircleFace(10, color), column=0)
			self.ts.legend.add_face(TextFace(name), column=1)

		if self.empty:
			self.ts.legend.add_face(CircleFace(10, self.colors[7]), column=0)
			self.ts.legend.add_face(TextFace('none'), column=1)

		self.ts.legend.add_face(TextFace(datetime.now().strftime('%a, %d.%m.%Y; %H:%M:%S'), fsize=8), column=1)
		self.ts.legend.add_face(TextFace('Tree: ' + self.treefile, fsize=8), column=1)
		self.ts.legend.add_face(TextFace('Crosshits: ' + self.combo[10:], fsize=8), column=1)
		self.ts.legend.add_face(TextFace('Pruning: ' + self.prune, fsize=8), column=1)

		# Order the graphical tree output from highly branched nodes on top and less branched nodes to the bottom
		self.t.ladderize(1)

	# end treeLayout


	def readPruningFile(self):
		self.pruneAfter = set()
		self.pruneAfterChildren = set()
		self.dontPrune = set()
		self.detach = set()
		self.delete = set()

		try:
			with open(self.prune, 'r') as f:
				active = self.pruneAfter
				for line in f:
					line = line.strip()

					if not line:
						continue

					if line.startswith('#'):
						continue

					if line.startswith('!'):
						keyword = line[1:].lower()
						if keyword == 'pruneafter':
							active = self.pruneAfter
						elif keyword == 'pruneafterchildren':
							active = self.pruneAfterChildren
						elif keyword == 'dontprune':
							active = self.dontPrune
						elif keyword == 'detach':
							active = self.detach
						elif keyword == 'delete':
							active = self.delete
						else:
							print('Keyword {} not found!'.format(keyword))

						continue

					line = line.split('#')[0].rstrip()

					if '^' in line:
						active.add(int(line.split('^')[1]))
					else:
						try:
							active.add(int(line))
						except ValueError:
							active.add(line)

		except IOError:
			print('File {} not found!'.format(self.prune))

	# end readPruningFile()


	def addFeatures(self):
		# Traverse the tree and add features to each node
		for n in self.t.traverse('postorder'):
			if n.name == 'NoName' or not n.name:
				n.name = 'Biota^1'

			n.add_features(plainName = n.name.split('^')[0])
			n.add_features(taxid = int(n.name.split('^')[1]))

			kwargs = {f: 0 for f in self.featurelist}
			n.add_features(**kwargs)

			n.add_features(f_none = 0)
			n.add_features(total = 0)

			if n.is_leaf():
				countThis = True
				try:
					#Tree is derived from one protein, has more entries than crosshits
					ev = self.elements[n.taxid]
					for f in self.featurelist:
						if ev in getattr(self, f):
							# n.'f' += 1
							setattr(n, f, getattr(n, f) + 1)
							break

				except KeyError:
					#no crosshit at this leaf
					if self.empty:
						n.f_none += 1
					else:
						countThis = False

				if countThis:
					n.total += 1


			else:
				for x in n.children:
					for f in self.featurelist:
						# n.'f' += x.'f'
						setattr(n, f, getattr(n, f) + getattr(x, f))

					n.f_none += x.f_none
					n.total += x.total

	# end addFeatures()


	def pruneTree(self):
		# Start pruning the tree

		# First, we need to collect the nodes to prune and afterwards prune them.
		# Otherwise, tree traversing doesn't work.
		toPrune = []
		toDetach = []
		toFinalDelete = []
		for n in self.t.traverse('postorder'):
			# Prune after children
			if n.taxid in self.pruneAfterChildren or n.plainName in self.pruneAfterChildren:
				for m in n.children:
					toPrune.append(m)

			# Prune directly
			if n.taxid in self.pruneAfter or n.plainName in self.pruneAfter:
				toPrune.append(n)

			# Detach directly
			if n.taxid in self.detach or n.plainName in self.detach:
				toDetach.append(n)

			# Delete in the end
			if n.taxid in self.delete or n.plainName in self.delete:
				toFinalDelete.append(n)

		# Now, we have the nodes to keep. We need a list of nodes to detach.
		for n in toPrune:
			for m in n.children:
				toDetach.append(m)

		# Finally, we can detach the nodes.
		for n in toDetach:
			if n.taxid not in self.dontPrune and n.plainName not in self.dontPrune:
				n.detach()

		# Now that the tree is pruned, we may go through it ones again and see, whether
		# we can collapse some nodes (if wanted by user)
		if self.collapse:
			toDelete = []
			for n in self.t.traverse('postorder'):
				if not n.is_leaf() and (len(n.get_sisters()) == 0 or len(n.children) == 1):
					toDelete.append(n)
				elif n.is_leaf() and n.total == 0 and n not in toPrune:
					toDelete.append(n)

			for n in toDelete:
				n.delete()

		for n in toFinalDelete:
			n.delete()

	# end pruneTree()


	def showTree(self):
		self.t.show(tree_style=self.ts)


	def renderTree(self, outfile, width, dpi):
		self.t.render(outfile, w = width, dpi = dpi, tree_style = self.ts)


if __name__ == '__main__':
	import argparse

	parser = argparse.ArgumentParser(description='Create a tree with piecharts. You need: Newick tree file, config file, attributes file, and a pruning file')

	parser.add_argument('combo', help='Filename of .tsv file with the crosshits of two proteins.')

	parser.add_argument('-t', '--tree', default='trees/general.tre', help='Filename of .tre file for the tree. The file must be in Newick format and the names in the format: name^taxid')
	parser.add_argument('-p', '--prune', default='tree_to_prune.txt', help='Filename of the file with information, which nodes to prune and delete.')
	parser.add_argument('--startnode', default=None, help='If the tree should not start from root, give the exact name of the node that should be the root of the tree.')
	parser.add_argument('--sqrt', action='store_true', help='Give this flag to use square root instead of log10 for the circle size')
	parser.add_argument('-n', '--nocollapse', action='store_true', help='If this flag is given, dont collapse the tree.')
	parser.add_argument('-e', '--empty', action='store_true', help='If this flag is given, show also species that do not have any of the given features.')
	parser.add_argument('-s', '--show', action='store_true', help='Show tree in ETE Tree Viewer (good for inspection)')
	parser.add_argument('-o', '--outfile', default='', help='If a filename is given, save the tree to that file. Valid extensions are: svg, pdf, png')
	parser.add_argument('-w', '--width', type=int, default=1000, help='Width of the resulting picture in pixels.')
	parser.add_argument('-b', '--bins', type=int, default=[5], nargs=1, help='Number of bins/features to create for e-value range')
	parser.add_argument('-r', '--range', type=int, nargs=2, default=[30,150], help='Range of evalues to use for bins')

	args = parser.parse_args()

	t = Tree(args.tree, format=8)

	TM = TreeMaker(tree = t, treefile = args.tree, combo=args.combo, prune = args.prune, tf=TaxFinder(), bins=args.range+args.bins, sqrt = args.sqrt, collapse = not args.nocollapse, show_empty_species = args.empty, startnode = args.startnode)

	if args.show:
		TM.showTree()

	if args.outfile:
		TM.renderTree(args.outfile, width=args.width, dpi=90)

	if not args.show and not args.outfile:
		print('The script finished successfully but no output was created.\nUse -o or -s to create some output.')