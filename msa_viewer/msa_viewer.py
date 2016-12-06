#!/usr/bin/env python3

import sys
import os
from PIL import Image, ImageDraw, ImageFont # install pillow, not pil!
from ete3 import Tree # For easy Newick tree traversal

#TODO: add_categories & add_clusters is missing

class MSA():

	def __init__(self, alignment_file, print_color=False, is_dna=False, save_fn='msa.png', border_width=10):
		self.alignment_file = alignment_file
		self.print_color = print_color
		self.is_dna = is_dna
		self.save_fn = save_fn
		self.border_width = border_width

		self.alignment = {}

		self.entity_colors = self._determine_colors()
		# Easily distinguishable colors from colorbrewer2.org
		self.cluster_colors = [(228,26,28),(55,126,184),(77,175,74),(152,78,163),(255,127,0),(255,255,51),(166,86,40)]
		if self.print_color:
			self.exon_colors = {False: (200,169,109), True: (154,151,220), 'notfound': (208,28,28), 'double': (28,208,28), 'gap': (255,255,255)}
		else:
			self.exon_colors = {False: (60,60,60), True: (120,120,120), 'notfound': (180,180,180), 'double': (0,0,0), 'gap': (255,255,255)}

		self.order = []
		self.category = {}
		self.clusters = {}
		self.exons = {}
		self.acc = None
		self.longest_name = 0
		self.overlay_acc = ''
		self.overlay_spans = []
		self.overlay_coords = []


	def make_alignment(self, overwrite_order = False, fontsize = 0):
		self.read_alignment_file(overwrite_order)
		self.print_alignment(fontsize)


	def read_alignment_file(self, overwrite_order = False):
		data = []
		alignment = []
		alignments = []
		current = None
		nongap_pos = 0
		overlay_pos = 0
		overlay_real_pos = 0
		overlay_coords = []

		overlay_spans = self.overlay_spans

		self.alignment = {}

		def bundle_alignment(data, alignment, name):
			return {'name': name, 'colors': tuple(data), 'letters': ''.join(alignment)}

		with open(self.alignment_file, 'r') as f:
			for line in f:
				line = line.rstrip()

				if line.startswith('>'):
					if data:
						self.alignment[acc] = bundle_alignment(data, alignment, current)
					data = []
					current = line[1:]
					pos = 0
					alignment = []
					# A bit too much hardcoding here I guess
					acc = '_'.join(current.split('_')[3:])
					self.longest_name = max(len(current), self.longest_name)
					if overwrite_order:
						self.order.append(acc)
				else:
					alignment.extend(line)
					for char in line:
						if self.exons:
							if char == '-':
								data.append(self.exon_colors['gap'])
							else:
								try:
									data.append(self.exons[acc][pos])
								except IndexError:
									data.append(self.exon_colors['notfound'])
								nongap_pos += 1
						else:
							data.append(self.entity_colors[char])


					if acc == self.overlay_acc:
						for char in line:
							try:
								if overlay_pos == overlay_spans[0][0]:
									overlay_coords.append(overlay_real_pos)
								elif overlay_pos == overlay_spans[0][1]:
									overlay_coords.append(overlay_real_pos)
									del(overlay_spans[0])
							except IndexError:
								break
							if char != '-':
								overlay_pos += 1
							overlay_real_pos += 1

		self.alignment[acc] = bundle_alignment(data, alignment, current)
		self.overlay_coords = [overlay_coords[n:n+2] for n in range(0, len(overlay_coords), 2)]


	def print_alignment(self, fontsize=0):
		if not self.alignment:
			raise ValueError('No alignment was read in. Read it by executing read_alignment_file()')

		if self.save_fn.endswith('svg'):
			self.print_svg(fontsize)
		elif self.save_fn.endswith('html') or self.save_fn.endswith('htm'):
			self.print_html(fontsize)
		else:
			self.print_pil(fontsize)


	def print_svg(self, fontsize=0):
		if not fontsize:
			raise NotImplementedError('SVG output is not defined for fontsize 0')

		width_to_height = 0.65
		y_correction = fontsize * 0.2
		x_correction = fontsize * 0.05

		names = []
		alignments = []
		total_colors = {}
		cur_color_code = 0
		delta_x = fontsize * width_to_height
		xoffset = self.longest_name * delta_x

		total_height = fontsize * len(self.order)
		total_width = xoffset + delta_x * len(self.alignment[self.order[0]]['colors'])

		for y, elem in enumerate(self.order):
			name = self.alignment[elem]['name']
			letters = self.alignment[elem]['letters']
			colors = self.alignment[elem]['colors']
			for c in set(colors):
				if c not in total_colors:
					total_colors[c] = chr(cur_color_code + ord('a'))
					cur_color_code += 1
			text = ['<text x="0" y="{}">{}</text>\n'.format(fontsize * (y+1) - y_correction, name)]
			oldcol = None
			stretch = []
			stretch_start = 0
			for x in range(len(colors)):
				if colors[x] == oldcol:
					stretch.append(letters[x])
				else:
					if oldcol:
						svg_stretch = []
						for i, char in enumerate(stretch):
							svg_stretch.append('\t\t<text x="{:.1f}" y="{:.1f}">{}</text>\n'.format(xoffset + (stretch_start+i) * delta_x + x_correction, (y+1) * fontsize - y_correction, char))
						text.append('\t<rect x="{:.1f}" y="{:.1f}" width="{:.1f}" height="{}" class="c{}" />\n{}'.format(xoffset + stretch_start * delta_x, y * fontsize, len(stretch) * delta_x, fontsize, total_colors[oldcol], ''.join(svg_stretch)))
					oldcol = colors[x]
					stretch = [letters[x]]
					stretch_start = x
			svg_stretch = []
			for i, char in enumerate(stretch):
				svg_stretch.append('\t\t<text x="{:.1f}" y="{:.1f}">{}</text>\n'.format(xoffset + (stretch_start+i) * delta_x + x_correction, (y+1) * fontsize - y_correction, char))
			text.append('\t<rect x="{:.1f}" y="{:.1f}" width="{:.1f}" height="{}" class="c{}" />\n{}'.format(xoffset + stretch_start * delta_x, y * fontsize, len(stretch) * delta_x, fontsize, total_colors[oldcol], ''.join(svg_stretch)))

			alignments.append(''.join(text))

		# The template file must be in the same folder as this script.
		# The complicated access is due to importability of this script.
		template_fn = font_fn = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'svg_template')
		template = open(template_fn, 'r').read()

		css_code = ['.c{} {{fill: rgb{};}}'.format(name, repr(value)) for value, name in total_colors.items()]

		template = template.format(WIDTH=total_width+1, HEIGHT=total_height, FONTSIZE=fontsize, CSS='\n'.join(css_code), ALIGNMENTS='\n'.join(alignments))

		with open(self.save_fn, 'w') as f:
			f.write(template)


	def print_html(self, fontsize=0):
		if not fontsize:
			raise NotImplementedError('HTML output is not defined for fontsize 0')

		names = []
		alignments = []
		total_colors = {}
		cur_color_code = 0
		for y, elem in enumerate(self.order):
			name = self.alignment[elem]['name']
			letters = self.alignment[elem]['letters']
			colors = self.alignment[elem]['colors']
			for c in set(colors):
				if c not in total_colors:
					total_colors[c] = chr(cur_color_code + ord('a'))
					cur_color_code += 1
			names.append('<nobr>{}&nbsp;</nobr>'.format(name))
			text = []
			oldcol = None
			stretch = []
			for x in range(len(colors)):
				if colors[x] == oldcol:
					stretch.append(letters[x])
				else:
					if oldcol:
						text.append('<span class="c{}">{}</span>'.format(total_colors[oldcol], ''.join(stretch)))
					oldcol = colors[x]
					stretch = [letters[x]]
			text.append('<span class="c{}">{}</span>'.format(total_colors[oldcol], ''.join(stretch)))

			alignments.append('<nobr>{}</nobr>'.format(''.join(text)))

		# The template file must be in the same folder as this script.
		# The complicated access is due to importability of this script.
		template_fn = font_fn = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'html_template')
		template = open(template_fn, 'r').read()

		css_code = ['.c{} {{background-color: rgb{};}}'.format(name, repr(value)) for value, name in total_colors.items()]

		template = template.format(FONTSIZE=fontsize, TITLE=os.path.basename(self.save_fn), CSS='\n'.join(css_code), NAMES='<br>\n'.join(names), ALIGNMENTS='<br>\n'.join(alignments))

		#template = template.replace('{FONTSIZE}', str(fontsize))
		#template = template.replace('{TITLE}', os.path.basename(self.save_fn))
		#template = template.replace('{CSS}', '\n'.join(css_code))
		#template = template.replace('{NAMES}', '<br>\n'.join(names))
		#template = template.replace('{ALIGNMENTS}', '<br>\n'.join(alignments))

		with open(self.save_fn, 'w') as f:
			f.write(template)



	def print_pil(self, fontsize=0):
		if fontsize:
			# The font file must be in the same folder as this script.
			# The complicated access is due to importability of this script.
			font_fn = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'FreeMono.ttf')

			# Read in font
			monofont = ImageFont.truetype(font_fn, fontsize)

			# Determine real font size (in pixels)
			font_width, font_height = ImageDraw.Draw(Image.new('RGB', (100, 100), (0,0,0))).textsize('X', font=monofont)
			xoffset = (self.longest_name + 1) * font_width
			font_height = int(font_height * 1.25)
			font_width = int(font_width * 1.4)
			yoffset = font_height if self.overlay_acc else 0
			height = len(self.order) * font_height + 2 * yoffset
			width = len(self.alignment[self.order[0]]['letters']) * font_width + xoffset #TODO: width depends on presence of categories & clusters

			im = Image.new('RGB', (width + 4*self.border_width, height), (255,255,255))
			draw = ImageDraw.Draw(im)

			for y, elem in enumerate(self.order):
				if self.border_width:
					#TODO: Categories should also be optional
					draw.rectangle(((self.border_width, y*font_height + yoffset), (2*self.border_width, (y+1)*font_height + yoffset)), self.category[elem])
					draw.rectangle(((width + 2*self.border_width, y*font_height + yoffset), (width + 3*self.border_width, (y+1)*font_height + yoffset)), self.category[elem])
					if self.clusters:
						draw.rectangle(((0, y*font_height + yoffset), (self.border_width, (y+1)*font_height + yoffset)), self.clusters[elem])
						draw.rectangle(((width + 3*self.border_width, y*font_height + yoffset), (width + 4*self.border_width, (y+1)*font_height + yoffset)), self.clusters[elem])
				draw.text((1 + 2*self.border_width, y*font_height + yoffset), self.alignment[elem]['name'], fill = (0,0,0), font = monofont)
				try:
					colors = self.alignment[elem]['colors']
					letters = self.alignment[elem]['letters']
					for x in range(len(colors)):
						draw.rectangle(((xoffset + x*font_width + 2*self.border_width, y*font_height + yoffset), (xoffset + (x+1)*font_width + 2*self.border_width, (y+1)*font_height + yoffset)), colors[x])
						draw.text((xoffset + x*font_width + 2*self.border_width, y*font_height + yoffset), letters[x], fill = (0,0,0), font = monofont)
				except KeyError:
					pass
		else:
			yoffset = 1 if self.overlay_acc else 0
			height = len(self.order) + 2 * yoffset
			width = len(self.alignment[self.order[0]]['letters'])

			im = Image.new('RGB', (width + 4*self.border_width, height), (255,255,255))

			for y, elem in enumerate(self.order):
				for x in range(self.border_width):
					# TODO: Categories should also be optional
					im.putpixel((x + self.border_width, y + yoffset), self.category[elem])
					im.putpixel((x + width + 2*self.border_width, y + yoffset), self.category[elem])

					if self.clusters:
						im.putpixel((x, y + yoffset), self.clusters[elem])
						im.putpixel((x + width + 3*self.border_width, y + yoffset), self.clusters[elem])
				try:
					for x, value in enumerate(self.alignment[elem]['colors']):
						im.putpixel((x + 2*self.border_width, y + yoffset), value)
				except KeyError:
					pass

		#if membranecoords:
		#	draw = ImageDraw.Draw(im, mode = 'RGBA')
		#	if fontsize:
		#		for x0, x1 in membranecoords:
		#			draw.rectangle(((xoffset + (x0-1)*font_width + 2*self.border_width, 0), (xoffset + x1*font_width + 2*self.border_width, height + yoffset)), (0,0,0,120))
		#	else:
		#		for x0, x1 in membranecoords:
		#			draw.rectangle(((x0 + 2*self.border_width, 0), (x1 + 2*self.border_width + 1, height + yoffset)), (0,0,0,120))

		im.save(self.save_fn)


	def _determine_colors(self):
		if not self.print_color:
			colors = {chr(x): (0, 0, 0) for x in range(ord('A'), ord('Z')+1)}
			colors['-'] = (255, 255, 255)
		elif self.is_dna:
			colors = {
				'A': (203, 0, 0),
				'T': (0, 0, 203),
				'U': (0, 0, 203),
				'G': (0, 203, 0),
				'C': (203, 0, 191),
				'N': (0, 0, 0),
				'-': (255, 255, 255)}
		else:
			colors = {
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
				'-': (255, 255, 255)}	# Gap

		return colors


	def _read_file(self, fn):
		with open(fn, 'r') as f:
			for line in f:
				line = line.split('#')[0].rstrip()
				if not line:
					continue
				yield line


	def add_order(self, fn):
		self.order = list(self._read_file(fn))


	def add_newick_info(self, fn, overwrite_order = False):
		mytree = Tree(fn)

		#if threshold:
		#	clusters = getClusters(mytree, threshold, clustercolors)
		#else:
		#	clusters = getFixedClusters(mytree, seeds, clustercolors)

		if overwrite_order:
			self.order = []
			for node in mytree.traverse('postorder'):
				if node.is_leaf():
					nname = node.name
					acc = '_'.join(nname.split('_')[3:])
					self.order.append(nname)


	def add_exons(self, fn):
		for line in self._read_file(fn):
			acc, e = line.split('\t')
			# Remove the outer brackets and spaces
			e = e[1:-1].replace(' ', '')
			# Turns the string e.g. (1,4),(4,8) to an actual tuple of tuples (basically the same as eval(e) would do)
			e = tuple(tuple(int(y) for y in x.split(',')) for x in e.split('),('))
			ex = [self.exon_colors['notfound'] for _ in range(e[-1][1])]
			current = False
			for exon in e:
				current = not current
				for pos in range(exon[0] - 1, exon[1]):
					if ex[pos] == self.exon_colors['notfound']:
						ex[pos] = self.exon_colors[current]
					else:
						ex[pos] = self.exon_colors['double']
			self.exons[acc] = tuple(ex)


	def add_overlays(self, fn):
		handle = self._read_file(fn)
		self.overlay_acc = next(handle)
		for line in handle:
			self.overlay_spans.append(tuple(map(int, line.split())))


	def add_cluster_info(self, threshold=0, seeds=None):
		raise NotImplementedError



if __name__ == '__main__':
	from argparse import ArgumentParser

	def toList(s):
		return s.split(',')

	parser = ArgumentParser(description='Show multiple sequence alignments pixel-wise')
	parser.add_argument('msa', help='The alignment file in fasta format')
	parser.add_argument('-o', '--order', help='A file with the order in which the alignments should be shown. One element per line. [None]')
	parser.add_argument('-n', '--newick', help='A newick file with all elements from which the order can be derived and clusters shown. [None]')
	parser.add_argument('-e', '--exon', help='A file with exons. Only proteins with known exon borders will be shown. [None]')
	parser.add_argument('-v', '--overlay', help='A file with overlay elements (like membrane-spanning). [None]')
	parser.add_argument('-t', '--thres', type=float, default=0, help='Distance threshold for when a new cluster shall be opened. [0]')
	parser.add_argument('-d', '--seeds', type=toList, default=None, help='Seeds for cluster determination, separated by comma. [None]')
	parser.add_argument('-b', '--border', type=int, default=10, help='Border width in pixels. [10]')
	parser.add_argument('-f', '--fontsize', type=int, default=0, help='Font size of letters. If font size is 0, only single pixels are shown [0]')
	parser.add_argument('-a', '--nucacid', action='store_true', help='Is it Nucleic acid? (Otherwise, it is a protein sequence.) [no/protein]')
	parser.add_argument('-c', '--color', action='store_true', help='Plot colors? [no]')
	parser.add_argument('-s', '--save', default='msa.png', help='Set the output file name. File ending determines file type. [msa.png]')
	args = parser.parse_args()

	alignment_file = args.msa
	order_file = args.order
	newick_file = args.newick
	exon_file = args.exon
	overlay_file = args.overlay
	threshold = args.thres
	seeds = args.seeds or []
	print_color = args.color
	fontsize = args.fontsize
	is_dna = args.nucacid
	save_fn = args.save
	border_width = args.border

	order_done = False

	msa = MSA(alignment_file, print_color, is_dna, save_fn, border_width)
	if order_file:
		msa.add_order(order_file)
		order_done = True
	if newick_file:
		msa.add_newick_info(newick_file, overwrite_order = not order_done)
		order_done = True
	if exon_file:
		msa.add_exons(exon_file)
	if overlay_file:
		msa.add_overlays(overlay_file)

	msa.make_alignment(overwrite_order = not order_done, fontsize = fontsize)
	#msa.print_alignment(fontsize = fontsize)
