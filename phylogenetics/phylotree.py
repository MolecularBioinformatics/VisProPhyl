#!/usr/bin/env python3

'''
This is the command line tool phylotree. This is not supposed to be imported.
'''

import argparse
import logging
import math
import os

from datetime import datetime
from operator import add

from ete3 import Tree, TreeStyle, faces, PieChartFace, TextFace, CircleFace


class TreeMaker():
	'''
	This class collects all information for showing or saving a tree.
	'''

	def __init__(self, tree, treefile, prune, config, attr, collapse=True, show_empty_species=True, start_node=None, count_total=False):
		'''
		Initialize the TreeMaker.

		:param tree: ete3.Tree object. The basis for the tree.
		:param treefile: The filename of the tree (only for the legend)
		:param prune: Filename of the file with information on how to
			prune the tree (automatically generated by `phylogenetics -i`)
		:param config: Filename of the config file (automatically
			generated by `phylogenetics -i`)
		:param attr: Filename of the attributes file. Generated by
			`phylogenetics`
		:param collapse: If True, collapse nodes with only a single child
		:param show_empty_species: If True, also show taxa that do not
			contain any of the given attributes
		:param start_node: If given, use this node name as root for the
			tree. If None (default), use the whole tree.
		:param count_total: TODO: What does this do???
		'''

		if start_node is None:
			self.tree = tree
		else:
			node = tree.search_nodes(name=start_node)
			if not node:
				raise KeyError(f'start_node "{start_node}" not found in tree.')
			self.tree = node[0]

		self.config = config
		self.attr = attr
		self.scaling_methods = {'log10', 'sqrt', 'linear'}
		self.scaling = 'log10'
		self.collapse = collapse
		self.show_empty_species = show_empty_species
		self.treefile = treefile
		self.prune = prune
		self.start_node = start_node
		self.count = count_total

		# Default colors for the features
		self.default_colors = [
			'#AA0000', # A
			'#00AA00', # B
			'#55AAFF', # C
			'#DDDD00', # A+B
			'#AA00AA', # A+C
			'#3333FF', # B+C
			'#000000', # A+B+C
			'#c8c8c8', # none
		]
		self.colors = self.default_colors[:]

		self.f_a = None  # Code for feature 1
		self.f_b = None  # Code for feature 2
		self.f_c = None  # Code for feature 3

		self.name_a = None  # Name for feature 1
		self.name_b = None  # Name for feature 2
		self.name_c = None  # Name for feature 3

		self.counter_elements = []

		self.read_attributes()
		self.read_config()
		self.add_features()
		self.read_pruning_file()
		self.prune_tree()
		self.tree_layout()


	def read_attributes(self):
		'''
		Read tree attributes into dict.
		'''

		# {1234: {'aA', 'aF', 'bE'}, ...}
		self.taxid2attr = {}
		self.name2attr = {}
		with open(self.attr) as attr_file:
			mode = 'name'
			for line in attr_file:
				lline = line.rstrip().split(maxsplit=1)

				if mode == 'name':
					if lline[0] == '---':
						mode = 'attrs'
					else:
						self.name2attr[lline[0]] = lline[1]

				else:
					taxid = int(lline[0])
					attributes = {lline[1][i:i+2] for i in range(0, len(lline[1]), 2)}
					self.taxid2attr[taxid] = attributes


	def read_config(self):
		'''
		Reads the config file (saved in self.config) and sets various
		configuration variables accordingly.
		'''

		self.counter_elements = []
		self.colors = []

		mode = ''
		modes = {
			'FEATURES': 'features',
			'COUNT': 'count',
			'COLORS': 'colors',
			'SCALE': 'scale',
		}

		with open(self.config) as config_file:
			for linenum, line in enumerate(config_file, start=1):
				line = line.split('#', maxsplit=1)[0].strip()
				if not line:
					continue

				if line in modes:
					mode = modes[line]

				elif mode == 'features':
					if line in self.name2attr:
						attr = self.name2attr[line]
					else:
						raise ValueError(f'tree_config l.{linenum}: Unknown attribute "{line}"')

					if self.f_a is None:
						self.f_a = attr
						self.name_a = line
					elif self.f_b is None:
						self.f_b = attr
						self.name_b = line
					elif self.f_c is None:
						self.f_c = attr
						self.name_c = line
					else:
						raise ValueError('tree_config l.{linenum}: More than three attributes given')

				elif mode == 'count':
					if line in self.name2attr:
						attr = self.name2attr[line]
					else:
						raise ValueError(f'tree_config l.{linenum}: Unknown attribute "{line}"')
					self.counter_elements.append(attr)

				elif mode == 'colors':
					self.colors.append(f'#{line}')

				elif mode == 'scale':
					line = line.lower()
					if line in self.scaling_methods:
						self.scaling = line
					else:
						raise ValueError(f'tree_config l.{linenum}: Unknown scaling method: {line}!')

				else:
					raise ValueError(f'tree_config l.{linenum}: No mode selected and "{line}" found!')

		num_colors = len(self.colors)
		if num_colors == 0:
			self.colors = self.default_colors[:]
		if num_colors > 8:
			logging.warning(f'Too many colors in tree_config. Found {num_colors}, expected 8.')
			self.colors = self.colors[:8]
		elif num_colors < 8:
			logging.warning(f'Too few colors in tree_config. Found {num_colors},'
			'expected 8. Filling with default colors.')
			self.colors.extend(self.default_colors[num_colors:])

		if not self.scaling:
			self.scaling = 'log10'


	def layout(self, node):
		'''
		Defines the layout of the given node.

		:param node: An ete3 Node to configure.
		'''

		font = 'Arial'
		fgcol = '#000000'

		try:
			percents = [
				round(100.0*node.f_a/node.total),
				round(100.0*node.f_b/node.total),
				round(100.0*node.f_c/node.total),
				round(100.0*node.f_ab/node.total),
				round(100.0*node.f_ac/node.total),
				round(100.0*node.f_bc/node.total),
				round(100.0*node.f_all/node.total),
				round(100.0*node.f_none/node.total),
			]
		except ZeroDivisionError:
			txt = TextFace(f' {node.plainName}', ftype=font, fgcolor=fgcol)
		else:
			idx = 0
			if sum(percents) != 100:
				sorted_percents = sorted(percents)[-1]
				idx = percents.index(sorted_percents)
				percents[idx] += 100 - sum(percents)

			if self.scaling == 'sqrt':
				size = int(math.sqrt(node.total)*5+1)
			elif self.scaling == 'linear':
				size = int(node.total + 5)
			else:
				size = int(math.log10(node.total)*20+10)


			piechart = PieChartFace(percents, size, size, self.colors)
			faces.add_face_to_node(piechart, node, 0, position="branch-top")

			if node.counter:
				cnt = float(sum(node.counter))/float(node.total)
				txt = TextFace(f' {node.plainName}\n {cnt:.2f}', ftype=font, fgcolor=fgcol)
			elif self.count:
				cnt = node.total
				txt = TextFace(f' {node.plainName}\n {cnt:d}', ftype=font, fgcolor=fgcol)
			else:
				txt = TextFace(f' {node.plainName}', ftype=font, fgcolor=fgcol)
			#txt.margin_right = -30

		if node.is_leaf():
			faces.add_face_to_node(txt, node, 0, position="branch-right")
		else:
			faces.add_face_to_node(txt, node, 0, position="branch-bottom")


	def tree_layout(self):
		'''
		Defines the layout of the whole tree and adds a legend.
		'''

		self.tstyle = TreeStyle()
		self.tstyle.show_scale = False
		self.tstyle.show_leaf_name = False
		self.tstyle.extra_branch_line_type = 0
		self.tstyle.extra_branch_line_color = '#000000'
		self.tstyle.layout_fn = self.layout
		self.tstyle.branch_vertical_margin = -10

		self.tstyle.legend_position = 4

		if self.name_a:
			self.tstyle.legend.add_face(CircleFace(10, self.colors[0]), column=0)
			self.tstyle.legend.add_face(TextFace(self.name_a), column=1)

		if self.name_b:
			self.tstyle.legend.add_face(CircleFace(10, self.colors[1]), column=0)
			self.tstyle.legend.add_face(TextFace(self.name_b), column=1)

		if self.name_c:
			self.tstyle.legend.add_face(CircleFace(10, self.colors[2]), column=0)
			self.tstyle.legend.add_face(TextFace(self.name_c), column=1)

		if self.name_a and self.name_b:
			self.tstyle.legend.add_face(CircleFace(10, self.colors[3]), column=0)
			self.tstyle.legend.add_face(TextFace(f'{self.name_a} + {self.name_b}'), column=1)

		if self.name_a and self.name_c:
			self.tstyle.legend.add_face(CircleFace(10, self.colors[4]), column=0)
			self.tstyle.legend.add_face(TextFace(f'{self.name_a} + {self.name_c}'), column=1)

		if self.name_b and self.name_c:
			self.tstyle.legend.add_face(CircleFace(10, self.colors[5]), column=0)
			self.tstyle.legend.add_face(TextFace(f'{self.name_b} + {self.name_c}'), column=1)

		if self.name_a and self.name_b and self.name_c:
			self.tstyle.legend.add_face(CircleFace(10, self.colors[6]), column=0)
			self.tstyle.legend.add_face(TextFace(f'{self.name_a} + {self.name_b} + {self.name_c}'), column=1)

		if self.show_empty_species:
			self.tstyle.legend.add_face(CircleFace(10, self.colors[7]), column=0)
			self.tstyle.legend.add_face(TextFace('none'), column=1)

		self.tstyle.legend.add_face(
			TextFace(datetime.now().strftime('%a, %Y-%m-%d %H:%M:%S'), fsize=8), column=1)
		self.tstyle.legend.add_face(TextFace(f'Path: {os.getcwd()}', fsize=8), column=1)
		self.tstyle.legend.add_face(TextFace(f'Tree: {self.treefile}', fsize=8), column=1)
		self.tstyle.legend.add_face(TextFace(f'Attributes: {self.attr}', fsize=8), column=1)
		self.tstyle.legend.add_face(TextFace(f'Config: {self.config}', fsize=8), column=1)
		self.tstyle.legend.add_face(TextFace(f'Pruning: {self.prune}', fsize=8), column=1)

		# Order the graphical tree output from highly branched nodes on
		# top and less branched nodes to the bottom
		self.tree.ladderize(1)


	def add_features(self):
		'''
		Traverses the tree and add features to each node.
		'''

		for node in self.tree.traverse('postorder'):
			if node.name == 'NoName' or not node.name:
				node.name = 'Biota^1'

			node_name, node_taxid = node.name.split('^', maxsplit=1)

			node.add_features(plainName = node_name)
			node.add_features(taxid = int(node_taxid))

			node.add_features(f_a = 0)
			node.add_features(f_b = 0)
			node.add_features(f_c = 0)
			node.add_features(f_ab = 0)
			node.add_features(f_ac = 0)
			node.add_features(f_bc = 0)
			node.add_features(f_all = 0)
			node.add_features(f_none = 0)
			node.add_features(total = 0)
			node.add_features(counter = [0]*len(self.counter_elements))
			if node.is_leaf():
				count_this = True
				elems = self.taxid2attr[node.taxid]

				if self.f_a in elems and self.f_b in elems and self.f_c in elems:
					node.f_all += 1
				elif self.f_b in elems and self.f_c in elems:
					node.f_bc += 1
				elif self.f_a in elems and self.f_c in elems:
					node.f_ac += 1
				elif self.f_a in elems and self.f_b in elems:
					node.f_ab += 1
				elif self.f_c in elems:
					node.f_c += 1
				elif self.f_b in elems:
					node.f_b += 1
				elif self.f_a in elems:
					node.f_a += 1
				elif self.show_empty_species:
					node.f_none += 1
				else:
					count_this = False

				if count_this:
					node.total += 1

					for i, count in enumerate(self.counter_elements):
						if count in elems:
							node.counter[i] += 1

			else:
				for child in node.children:
					node.f_a += child.f_a
					node.f_b += child.f_b
					node.f_c += child.f_c
					node.f_ab += child.f_ab
					node.f_ac += child.f_ac
					node.f_bc += child.f_bc
					node.f_all += child.f_all
					node.f_none += child.f_none
					node.total += child.total

					node.counter = list(map(add, node.counter, child.counter))


	def read_pruning_file(self):
		'''
		Reads the pruning file (if any) and saves the pruning configuration.
		'''

		self.prune_after = set()
		self.prune_after_children = set()
		self.dont_prune = set()
		self.detach = set()
		self.delete = set()

		try:
			with open(self.prune) as prune_file:
				active = self.prune_after
				for linenum, line in enumerate(prune_file, start=1):
					line = line.split('#', maxsplit=1)[0].rstrip()

					if not line:
						continue

					if line.startswith('!'):
						keyword = line[1:].lower()
						if keyword == 'pruneafter':
							active = self.prune_after
						elif keyword == 'pruneafterchildren':
							active = self.prune_after_children
						elif keyword == 'dontprune':
							active = self.dont_prune
						elif keyword == 'detach':
							active = self.detach
						elif keyword == 'delete':
							active = self.delete
						else:
							raise ValueError(f'tree_to_prune l. {linenum}: Keyword {keyword} not found!')

						continue

					if '^' in line:
						active.add(int(line.split('^', maxsplit=1)[1]))
					else:
						try:
							active.add(int(line))
						except ValueError:
							active.add(line)

		except IOError:
			logging.warning(f'Pruning file "{self.prune}" not found!')


	def prune_tree(self):
		'''
		Prunes the tree according to saved options (likely loaded by
		`self.read_pruning_file()`)
		'''

		# First, we need to collect the nodes to prune and afterwards prune them.
		# Otherwise, tree traversing doesn't work.
		to_prune = []
		to_detach = []
		to_final_delete = []
		for node in self.tree.traverse('postorder'):
			# Prune after children
			if node.taxid in self.prune_after_children or node.plainName in self.prune_after_children:
				for child in node.children:
					to_prune.append(child)

			# Prune directly
			if node.taxid in self.prune_after or node.plainName in self.prune_after:
				to_prune.append(node)

			# Detach directly
			if node.taxid in self.detach or node.plainName in self.detach:
				to_detach.append(node)

			# Delete in the end
			if node.taxid in self.delete or node.plainName in self.delete:
				to_final_delete.append(node)

		# Now, we have the nodes to keep. We need a list of nodes to detach.
		for node in to_prune:
			for child in node.children:
				to_detach.append(child)

		# Finally, we can detach the nodes.
		for node in to_detach:
			if node.taxid not in self.dont_prune and node.plainName not in self.dont_prune:
				node.detach()

		# Now that the tree is pruned, we may go through it ones again
		# and see, whether we can collapse some nodes (if wanted by user).

		# Node attributes that should lead to collapsing
		# (= no information in node) [marked *]
		#
		# * len(n.children) == 1 -> reduce long singular lineages to one
		#   step, leaves dont have children
		# 	- if not self.show_empty_species only children with n.total > 0 have to
		#     be considered: len([x for x in n.children if x.total])
		#	- the root node can be 'deleted', therefore at that point
		#     complete collapsing fails right now
		#
		# * n.total == 0 [and n not in to_prune] -> empty nodes,
		#   (Pruning makes new leaves, keeping these) [^1]
		#	- this obviously doesn't apply showing empty nodes
		#	- I think in principle pruning commands should not affect
		#     collapsing (collapse should only make tree more readable,
		#     but dont change information). If (default) empty nodes
		#     should disappear, what about nodes that are empty & pruned?
		#     Keep them just because they were mentioned ? (currently for
		#     leaves made by pruning ^1). It would probably be more
		#     consequent to collapse emtpy nodes regardless of pruning.
		#     Maybe a '!dontcollapse' option for specific nodes could be
		#     introduced for the pruningfile

		if self.collapse:
			to_delete = []
			for node in self.tree.traverse('postorder'):

				if node.is_leaf():
					if not self.show_empty_species and node.total == 0 and node not in to_prune:
						to_delete.append(node)
				elif self.show_empty_species:
					if len(node.children) == 1:
						to_delete.append(node)
				elif node.total == 0 or len([x for x in node.children if x.total]) == 1:
					to_delete.append(node)

			for node in to_delete:
				node.delete()

		for node in to_final_delete:
			node.delete()


	def show_tree(self):
		'''
		Shows the tree in an interactive window (provided by ete3)
		'''

		self.tree.show(tree_style=self.tstyle)


	def render_tree(self, outfile, width, dpi=90):
		'''
		Redners the tree to a file.

		:param outfile: The file to save the tree to. The filename extension
			defines the file type. Valid extensions are defined by ete3,
			but svg, pdf, and png should be fine.
		:param width: The width of the file in pixels (height is calculated
			from the tree size)
		:param dpi: Image resolution in dots per inch.
		'''

		self.tree.render(outfile, w=width, dpi=dpi, tree_style=self.tstyle)


def main():
	'''
	This is the start point for the command line tool.
	'''

	parser = argparse.ArgumentParser(
		description='Create a tree with protein distribution as piecharts.')

	# General options
	parser.add_argument('-t', '--tree', default='trees/general.tre',
		help='Filename of .tre file for the tree. The file must be in '
		'Newick format and the names in the format: name^taxid '
		'[default: %(default)s]')

	parser.add_argument('-p', '--prune', default='tree_to_prune.txt',
		help='Filename of the file with information, which nodes to prune '
		'and delete. [default: %(default)s]')

	parser.add_argument('-c', '--config', default='tree_config.txt',
		help='Filename of the file that contains the config with the '
		'tree features and the colors. [default: %(default)s]')

	parser.add_argument('-a', '--attr', default='attributes.txt',
		help='Filename of the file with the tree attributes. '
		'[default: %(default)s]')

	parser.add_argument('--start_node', default=None,
		help='If the tree should not start from root, give the exact name '
		'of the node that should be the root of the tree.')

	parser.add_argument('-n', '--nocollapse', action='store_true',
		help='If this flag is given, dont collapse the tree.')

	parser.add_argument('-e', '--empty', action='store_true',
		help='If this flag is given, show also species that do not have '
		'any of the given features.')

	parser.add_argument('-s', '--show', action='store_true',
		help='Show tree in ETE Tree Viewer (good for inspection)')

	parser.add_argument('--count_total', action='store_true',
		help='Print the total number of hits in addition to Piechart size.')

	parser.add_argument('-o', '--outfile', default='',
		help='If a filename is given, save the tree to that file. '
		'Valid extensions are: svg, pdf, png')

	parser.add_argument('-w', '--width', type=int, default=1000,
		help='Width of the resulting picture in pixels. [default: %(default)s]')

	args = parser.parse_args()

	tree = Tree(args.tree, format=8)

	tree_maker = TreeMaker(
		tree=tree,
		treefile=args.tree,
		prune=args.prune,
		collapse=not args.nocollapse,
		show_empty_species=args.empty,
		start_node=args.start_node,
		count_total=args.count_total,
		config=args.config,
		attr=args.attr,
	)

	if args.show or not args.outfile:
		tree_maker.show_tree()

	if args.outfile:
		tree_maker.render_tree(args.outfile, width=args.width)
