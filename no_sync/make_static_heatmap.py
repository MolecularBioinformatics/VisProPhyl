#!/usr/bin/env python3

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
<svg viewBox="0 0 {width} {height}" width="{width}" height="{height}">
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
	taxa_to_show = ['Homo_sapiens_(Primate)^9606',
		'Mus_musculus_(Rodent)^10116',
		'Gallus_gallus_(Bird)^9031',
		'Chelonia_mydas_(Turtle)^8469',
		'Python_bivittatus_(Snake)^176946',
		'Xenopus_laevis_(Frog)^8355',
		'Danio_rerio_(Fish)^7955',
		'Branchiostoma_floridae_(Lancelet)^7739',
		'Strongylocentrotus_purpuratus_(Sea_urchin)^7668',
		'Drosophila_melanogaster_(Fly)^7227',
		'Caenorhabditis_elegans_(Nematode)^6239',
		'Helobdella_robusta_(Annelid)^6412',
		'Schistosoma_haematobium_(Platyhelminthes)^6185',
		'Saccharomyces_cerevisiae_(Fungi)^4932',
		'Arabidopsis_thaliana_(Plant)^3702',
		'Acanthamoeba_castellanii_(Amoeba)^5755']
	proteins_to_show = ['PRMT1',
		'PRMT2',
		'PRMT3',
		'PRMT4',
		'PRMT5',
		'PRMT6',
		'PRMT7',
		'PRMT8',
		'PRMT9',
		'PRMT11',
		'mTOR']

	staticHeatmap(taxa_to_show, proteins_to_show, None)
