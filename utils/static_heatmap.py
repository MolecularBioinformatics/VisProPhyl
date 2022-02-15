#!/usr/bin/env python3

'''
This script is meant to be changed to fit your needs.
It will write a static heatmap as svg to STDOUT.
Run in a folder, where a whole Phylogenetics run is finished like this:
$ python3 static_heatmap.py > heatmap.svg
'''


import sys

from PIL import ImageFont



def get_string_pixel_length(string, font='LiberationSans-Regular', size=20):
	container = ImageFont.truetype(font, size)
	dimensions = container.getsize(string)
	return dimensions[0]


def staticHeatmap(taxa, proteins, cutoff = None):
	'''
	Create a static heatmap as SVG. Taxa and proteins are lists of strings. Taxa must be in the form "printable_name^taxid" where underscores will be replaced by spaces. Both taxa and proteins will be shown in order of the list.
	If cutoff is None (default), the heatmap will be shown as greyscale (the better the hit, the darker). Else, the heatmap will be shown black/white where values above the cutoff will be filled. The cutoff must be an integer, the negative log of the evalue (e.g. 1e-30 → cutoff = 30).
	'''


	# width, height, proteins, taxa, squares
	template = '''<?xml version="1.0" encoding="UTF-8" standalone="no"?>
<!DOCTYPE svg PUBLIC "-//W3C//DTD SVG 1.1//EN" "http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd">
<svg xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink" viewBox="0 0 {width} {height}" width="100%" height="100%">
{proteins}
{taxa}
{squares}
</svg>
'''

	# y, x, text
	temp_protein = '''<text
transform="rotate(-45 {rotx} {roty})" x="{x}" y="{y}"
style="font-style:normal;font-weight:normal;font-size:35px;line-height:125%;font-family:sans-serif;letter-spacing:0px;word-spacing:0px;fill:#000000;fill-opacity:1;stroke:none;stroke-width:1px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1"
xml:space="preserve">
<tspan style="font-size:20px" x="{x}" y="{y}">{text}</tspan></text>'''

	# y, x, text
	temp_taxon = '''<text x="{x}" y="{y}"
style="font-style:normal;font-weight:normal;font-size:35px;line-height:125%;font-family:sans-serif;letter-spacing:0px;word-spacing:0px;fill:#000000;fill-opacity:1;stroke:none;stroke-width:1px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1"
xml:space="preserve">
<tspan style="font-size:20px" x="{x}" y="{y}">{text}</tspan></text>'''

	# y, x, size, color
	temp_square = '''<rect x="{x}" y="{y}" height="{size}" width="{size}"
style="fill:rgb{color};fill-opacity:1;stroke:none;" />'''


	taxids = []
	taxnames = []
	for i, tax in enumerate(taxa):
		taxname, taxid = tax.split('^')
		taxnames.append(taxname.replace('_', ' '))
		taxids.append(taxid)

	matrix = {}
	for i, protein in enumerate(proteins):
		matrix[protein] = {x: 0 for x in taxids}
		with open('interactivetables/{}.tsv'.format(protein), 'r') as f:
			for line in f:
				lline = line.split()
				if lline[0] in taxids:
					matrix[protein][lline[0]] = int(lline[1])

	width = 18
	dist = 10

	upper_limit = 150

	print_proteins = []
	print_taxa = []
	print_squares = []

	offset_x_squares = dist
	offset_x_proteins = offset_x_squares + int(dist*1.4)
	offset_x_taxa = len(matrix) * (width + dist) + offset_x_squares

	offset_y_proteins = int(max(get_string_pixel_length(x) for x in proteins) * 0.9)
	offset_y_squares = offset_y_proteins + int(dist*0.7)
	offset_y_taxa = offset_y_squares + int(width*0.9)

	for i, taxon in enumerate(taxnames):
		t = temp_taxon.format(x=offset_x_taxa, y=offset_y_taxa + i*(width+dist), text=taxon)
		print_taxa.append(t)

	for a, protein in enumerate(proteins):
		xpos = offset_x_proteins + a*(width+dist)
		ypos = offset_y_proteins
		rotx = xpos
		roty = ypos
		t = temp_protein.format(x=xpos, y=ypos, rotx=rotx, roty=roty, text=protein)
		print_proteins.append(t)

		for b, taxid in enumerate(taxids):
			if cutoff is None:
				e = 255 - int((min(matrix[protein][taxid], upper_limit) / 150)*255)
				c = (e, e, e)
			else:
				if matrix[protein][taxid] < cutoff:
					c = (255,255,255)
				else:
					c = (0,0,0)
			t = temp_square.format(y=offset_y_squares + b*(width + dist), x=offset_x_squares + a*(width + dist), size=width, color=c)
			print_squares.append(t)

	total_width = offset_x_taxa + dist + max(get_string_pixel_length(x) for x in taxa)
	total_height = offset_y_squares + (width + dist)*len(taxa)

	final = template.format(width=total_width, height=total_height, proteins='\n'.join(print_proteins), taxa='\n'.join(print_taxa), squares='\n'.join(print_squares))

	print(final)


if __name__ == '__main__':

	taxa_to_show = []
	with open('heatmap_config.txt') as heatmap_file:

		active = False
		keywords = {'TAXA_TO_SHOW', 'COLORS', 'ALGO'}

		for line in heatmap_file:
			line = line.split('#', maxsplit=1)[0].strip()
			if not line:
				continue

			if not active and line == 'TAXA_TO_SHOW':
				active = True
				continue

			if line in keywords:
				break

			self.taxa_to_show.append(line)

	proteins_to_show = []
	with open('proteinlist.txt') as protein_file:

		for line in protein_file:
			line = line.split('#', maxsplit=1)[0].strip()
			if not line:
				continue

			name = line.split()[0]

			self.proteins_to_show.append(name)


	staticHeatmap(taxa_to_show, proteins_to_show, None)
