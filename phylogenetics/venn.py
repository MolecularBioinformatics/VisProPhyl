#!/usr/bin/env python3

'''
Module to create Venn diagrams. Import `venn` to use it.

This code was derived from https://github.com/tctianchi/pyvenn
'''

import sys

from itertools import chain

import matplotlib.pyplot as plt
from matplotlib import patches

LOGIC = 1
NUMBER = 2
PERCENT = 4

default_colors = [
	# r, g, b, a
	[92, 192, 98, 0.5],
	[90, 155, 212, 0.5],
	[246, 236, 86, 0.6],
	[241, 90, 96, 0.4],
	[255, 117, 0, 0.3],
	[82, 82, 190, 0.2]]

default_colors = [[i[0] / 255.0, i[1] / 255.0, i[2] / 255.0, i[3]] for i in default_colors]

text_colors = [[i[0], i[1], i[2], 1] for i in default_colors]


def _draw_ellipse(axis, x, y, width, height, angle, fillcolor):
	'''
	Draw an ellipse on axis ax, at position x, y, with given width, height, angle, and fillcolor
	'''

	ellipse = patches.Ellipse(
		xy=(x, y),
		width=width,
		height=height,
		angle=angle,
		color=fillcolor)
	axis.add_patch(ellipse)


def _draw_triangle(axis, x1, y1, x2, y2, x3, y3, fillcolor):
	'''
	Draw a triangle on axis ax with corners (x1, y1), (x2, y2), (x3, y3) and fillcolor
	'''

	xy = [
		(x1, y1),
		(x2, y2),
		(x3, y3)
	]
	polygon = patches.Polygon(
		xy=xy,
		closed=True,
		color=fillcolor)
	axis.add_patch(polygon)


def _draw_text(axis, x, y, text, color=(0, 0, 0, 1)):
	'''
	Draw a text on axis ax at position x, y and a given color
	'''

	axis.text(
		x, y, text,
		horizontalalignment='center',
		verticalalignment='center',
		fontsize=14,
		color=color)


def get_labels(data, fill=NUMBER):
	'''
	Get a dict of labels for groups in data

	@type data: list[Iterable]
	@rtype: dict[str, str]

	:param data: List of lists with the elements to analyse
	:param fill: NUMBER, LOGIC, or PERCENT to define what label to show. Can be chained with |

	:returns: a dict of labels for all possible sets

	:example:
	>>> get_labels([range(10), range(5,15), range(3,8)], fill=NUMBER)

		{'001': '0',
		 '010': '5',
		 '011': '0',
		 '100': '3',
		 '101': '2',
		 '110': '2',
		 '111': '3'}
	'''

	num_data = len(data)

	sets_data = [set(data[i]) for i in range(num_data)]   # sets for separate groups
	s_all = set(chain(*data))						 # union of all sets

	set_collections = {}
	formatter = f'0{num_data}b'
	for num in range(1, 2**num_data):
		key = format(num, formatter)
		value = s_all
		sets_for_intersection = [sets_data[i] for i in range(num_data) if key[i] == '1']
		sets_for_difference = [sets_data[i] for i in range(num_data) if key[i] == '0']
		for current_set in sets_for_intersection:
			value = value & current_set
		for current_set in sets_for_difference:
			value = value - current_set
		set_collections[key] = value

	labels = {key: '' for key in set_collections}
	if fill & LOGIC:
		for key in set_collections:
			labels[key] = f'{key}: '
	if fill & NUMBER:
		for key, value in set_collections.items():
			labels[key] += str(len(value))
	if fill & PERCENT:
		data_size = len(s_all)
		for key, value in set_collections.items():
			perc = 100 * len(value) / data_size
			labels[key] += f'({perc:.1f}%)'

	return labels


def _venn2(labels, names=('A', 'B'), colors=None, figsize=(9, 7), dpi=96):
	'''
	Plot a 2-set Venn diagram

	:param labels: a label dict where keys are identified via binary
		codes ('01', '10', '11'), hence a valid set could look like:
		{'01': 'text 1', '10': 'text 2', '11': 'text 3'}. Unmentioned
		codes are considered as ''.
	:param names: List with the names of the groups
	:param colors: List of two lists with colors for the Venn elements.
		Colors must be defined like [r, g, b, alpha], where all four
		elements must be floats between 0 and 1.
	:param figsize: Tuple with two numbers giving the size of the figure
	:param dpi: int with dpi

	:returns: Tuple with Figure and AxesSubplot object
	'''

	if colors is None:
		colors = default_colors

	fig = plt.figure(figsize=figsize, dpi=dpi)
	axis = fig.add_subplot(111, aspect='equal')
	axis.set_axis_off()
	axis.set_ylim(bottom=0.0, top=0.7)
	axis.set_xlim(left=0.0, right=1.0)

	# body
	_draw_ellipse(axis, 0.375, 0.3, 0.5, 0.5, 0.0, colors[0])
	_draw_ellipse(axis, 0.625, 0.3, 0.5, 0.5, 0.0, colors[1])
	_draw_text(axis, 0.74, 0.30, labels.get('01', ''))
	_draw_text(axis, 0.26, 0.30, labels.get('10', ''))
	_draw_text(axis, 0.50, 0.30, labels.get('11', ''))

	# legend
	_draw_text(axis, 0.20, 0.56, names[0], text_colors[0])
	_draw_text(axis, 0.80, 0.56, names[1], text_colors[1])
	leg = axis.legend(names, loc='best', fancybox=True)
	leg.get_frame().set_alpha(0.5)

	return fig, axis


def _venn3(labels, names=('A', 'B', 'C'), colors=None, figsize=(9, 9), dpi=96):
	'''
	Plot a 3-set Venn diagram

	:param labels: a label dict where keys are identified via binary
		codes ('001', '010', '011', ...), hence a valid set could look
		like: {'001': 'text 1', '010': 'text 2', '011': 'text 3', ...}.
		Unmentioned codes are considered as ''.
	:param names: List with the names of the groups
	:param colors: List of three lists with colors for the Venn elements.
		Colors must be defined like [r, g, b, alpha], where all four
		elements must be floats between 0 and 1.
	:param figsize: Tuple with two numbers giving the size of the figure
	:param dpi: int with dpi

	:returns: Tuple with Figure and AxesSubplot object
	'''

	if colors is None:
		colors = default_colors

	fig = plt.figure(figsize=figsize, dpi=dpi)
	axis = fig.add_subplot(111, aspect='equal')
	axis.set_axis_off()
	axis.set_ylim(bottom=0.0, top=1.0)
	axis.set_xlim(left=0.0, right=1.0)

	# body
	_draw_ellipse(axis, 0.333, 0.633, 0.5, 0.5, 0.0, colors[0])
	_draw_ellipse(axis, 0.666, 0.633, 0.5, 0.5, 0.0, colors[1])
	_draw_ellipse(axis, 0.500, 0.310, 0.5, 0.5, 0.0, colors[2])
	_draw_text(axis, 0.50, 0.27, labels.get('001', ''))
	_draw_text(axis, 0.73, 0.65, labels.get('010', ''))
	_draw_text(axis, 0.61, 0.46, labels.get('011', ''))
	_draw_text(axis, 0.27, 0.65, labels.get('100', ''))
	_draw_text(axis, 0.39, 0.46, labels.get('101', ''))
	_draw_text(axis, 0.50, 0.65, labels.get('110', ''))
	_draw_text(axis, 0.50, 0.51, labels.get('111', ''))

	# legend
	_draw_text(axis, 0.15, 0.87, names[0], text_colors[0])
	_draw_text(axis, 0.85, 0.87, names[1], text_colors[1])
	_draw_text(axis, 0.50, 0.02, names[2], text_colors[2])
	leg = axis.legend(names, loc='best', fancybox=True)
	leg.get_frame().set_alpha(0.5)

	return fig, axis


def _venn4(labels, names=('A', 'B', 'C', 'D'), colors=None, figsize=(12, 12), dpi=96):
	'''
	Plot a 4-set Venn diagram

	:param labels: a label dict where keys are identified via binary
		codes ('0001', '0010', '0011', ...), hence a valid set could look
		like: {'0001': 'text 1', '0010': 'text 2', '0011': 'text 3', ...}.
		Unmentioned codes are considered as ''.
	:param names: List with the names of the groups
	:param colors: List of four lists with colors for the Venn elements.
		Colors must be defined like [r, g, b, alpha], where all four
		elements must be floats between 0 and 1.
	:param figsize: Tuple with two numbers giving the size of the figure
	:param dpi: int with dpi

	:returns: Tuple with Figure and AxesSubplot object
	'''

	if colors is None:
		colors = default_colors

	fig = plt.figure(figsize=figsize, dpi=dpi)
	axis = fig.add_subplot(111, aspect='equal')
	axis.set_axis_off()
	axis.set_ylim(bottom=0.0, top=1.0)
	axis.set_xlim(left=0.0, right=1.0)

	# body
	_draw_ellipse(axis, 0.350, 0.400, 0.72, 0.45, 140.0, colors[0])
	_draw_ellipse(axis, 0.450, 0.500, 0.72, 0.45, 140.0, colors[1])
	_draw_ellipse(axis, 0.544, 0.500, 0.72, 0.45, 40.0, colors[2])
	_draw_ellipse(axis, 0.644, 0.400, 0.72, 0.45, 40.0, colors[3])
	_draw_text(axis, 0.85, 0.42, labels.get('0001', ''))
	_draw_text(axis, 0.68, 0.72, labels.get('0010', ''))
	_draw_text(axis, 0.77, 0.59, labels.get('0011', ''))
	_draw_text(axis, 0.32, 0.72, labels.get('0100', ''))
	_draw_text(axis, 0.71, 0.30, labels.get('0101', ''))
	_draw_text(axis, 0.50, 0.66, labels.get('0110', ''))
	_draw_text(axis, 0.65, 0.50, labels.get('0111', ''))
	_draw_text(axis, 0.14, 0.42, labels.get('1000', ''))
	_draw_text(axis, 0.50, 0.17, labels.get('1001', ''))
	_draw_text(axis, 0.29, 0.30, labels.get('1010', ''))
	_draw_text(axis, 0.39, 0.24, labels.get('1011', ''))
	_draw_text(axis, 0.23, 0.59, labels.get('1100', ''))
	_draw_text(axis, 0.61, 0.24, labels.get('1101', ''))
	_draw_text(axis, 0.35, 0.50, labels.get('1110', ''))
	_draw_text(axis, 0.50, 0.38, labels.get('1111', ''))

	# legend
	_draw_text(axis, 0.13, 0.18, names[0], text_colors[0])
	_draw_text(axis, 0.18, 0.83, names[1], text_colors[1])
	_draw_text(axis, 0.82, 0.83, names[2], text_colors[2])
	_draw_text(axis, 0.87, 0.18, names[3], text_colors[3])
	leg = axis.legend(names, loc='best', fancybox=True)
	leg.get_frame().set_alpha(0.5)

	return fig, axis


def _venn5(labels, names=('A', 'B', 'C', 'D', 'E'), colors=None, figsize=(13, 13), dpi=96):
	'''
	Plot a 5-set Venn diagram

	:param labels: a label dict where keys are identified via binary
		codes ('00001', '00010', '00011', ...), hence a valid set could
		look like: {'00001': 'text 1', '00010': 'text 2', '00011': 'text 3', ...}.
		Unmentioned codes are considered as ''.
	:param names: List with the names of the groups
	:param colors: List of five lists with colors for the Venn elements.
		Colors must be defined like [r, g, b, alpha], where all four
		elements must be floats between 0 and 1.
	:param figsize: Tuple with two numbers giving the size of the figure
	:param dpi: int with dpi

	:returns: Tuple with Figure and AxesSubplot object
	'''

	if colors is None:
		colors = default_colors

	fig = plt.figure(figsize=figsize, dpi=dpi)
	axis = fig.add_subplot(111, aspect='equal')
	axis.set_axis_off()
	axis.set_ylim(bottom=0.0, top=1.0)
	axis.set_xlim(left=0.0, right=1.0)

	# body
	_draw_ellipse(axis, 0.428, 0.449, 0.87, 0.50, 155.0, colors[0])
	_draw_ellipse(axis, 0.469, 0.543, 0.87, 0.50, 82.0, colors[1])
	_draw_ellipse(axis, 0.558, 0.523, 0.87, 0.50, 10.0, colors[2])
	_draw_ellipse(axis, 0.578, 0.432, 0.87, 0.50, 118.0, colors[3])
	_draw_ellipse(axis, 0.489, 0.383, 0.87, 0.50, 46.0, colors[4])
	_draw_text(axis, 0.27, 0.11, labels.get('00001', ''))
	_draw_text(axis, 0.72, 0.11, labels.get('00010', ''))
	_draw_text(axis, 0.55, 0.13, labels.get('00011', ''))
	_draw_text(axis, 0.91, 0.58, labels.get('00100', ''))
	_draw_text(axis, 0.78, 0.64, labels.get('00101', ''))
	_draw_text(axis, 0.84, 0.41, labels.get('00110', ''))
	_draw_text(axis, 0.76, 0.55, labels.get('00111', ''))
	_draw_text(axis, 0.51, 0.90, labels.get('01000', ''))
	_draw_text(axis, 0.39, 0.15, labels.get('01001', ''))
	_draw_text(axis, 0.42, 0.78, labels.get('01010', ''))
	_draw_text(axis, 0.50, 0.15, labels.get('01011', ''))
	_draw_text(axis, 0.67, 0.76, labels.get('01100', ''))
	_draw_text(axis, 0.70, 0.71, labels.get('01101', ''))
	_draw_text(axis, 0.51, 0.74, labels.get('01110', ''))
	_draw_text(axis, 0.64, 0.67, labels.get('01111', ''))
	_draw_text(axis, 0.10, 0.61, labels.get('10000', ''))
	_draw_text(axis, 0.20, 0.31, labels.get('10001', ''))
	_draw_text(axis, 0.76, 0.25, labels.get('10010', ''))
	_draw_text(axis, 0.65, 0.23, labels.get('10011', ''))
	_draw_text(axis, 0.18, 0.50, labels.get('10100', ''))
	_draw_text(axis, 0.21, 0.37, labels.get('10101', ''))
	_draw_text(axis, 0.81, 0.37, labels.get('10110', ''))
	_draw_text(axis, 0.74, 0.40, labels.get('10111', ''))
	_draw_text(axis, 0.27, 0.70, labels.get('11000', ''))
	_draw_text(axis, 0.34, 0.25, labels.get('11001', ''))
	_draw_text(axis, 0.33, 0.72, labels.get('11010', ''))
	_draw_text(axis, 0.51, 0.22, labels.get('11011', ''))
	_draw_text(axis, 0.25, 0.58, labels.get('11100', ''))
	_draw_text(axis, 0.28, 0.39, labels.get('11101', ''))
	_draw_text(axis, 0.36, 0.66, labels.get('11110', ''))
	_draw_text(axis, 0.51, 0.47, labels.get('11111', ''))

	# legend
	_draw_text(axis, 0.02, 0.72, names[0], text_colors[0])
	_draw_text(axis, 0.72, 0.94, names[1], text_colors[1])
	_draw_text(axis, 0.97, 0.74, names[2], text_colors[2])
	_draw_text(axis, 0.88, 0.05, names[3], text_colors[3])
	_draw_text(axis, 0.12, 0.05, names[4], text_colors[4])
	leg = axis.legend(names, loc='best', fancybox=True)
	leg.get_frame().set_alpha(0.5)

	return fig, axis


def _venn6(labels, names=('A', 'B', 'C', 'D', 'E'), colors=None, figsize=(20, 20), dpi=96):
	'''
	Plot a 6-set Venn diagram

	:param labels: a label dict where keys are identified via binary
		codes ('00001', '00010', '00011', ...), hence a valid set could
		look like: {'00001': 'text 1', '00010': 'text 2', '00011': 'text 3', ...}.
		Unmentioned codes are considered as ''.
	:param names: List with the names of the groups
	:param colors: List of six lists with colors for the Venn elements.
		Colors must be defined like [r, g, b, alpha], where all four
		elements must be floats between 0 and 1.
	:param figsize: Tuple with two numbers giving the size of the figure
	:param dpi: int with dpi

	:returns: Tuple with Figure and AxesSubplot object
	'''

	if colors is None:
		colors = default_colors

	fig = plt.figure(figsize=figsize, dpi=dpi)
	axis = fig.add_subplot(111, aspect='equal')
	axis.set_axis_off()
	axis.set_ylim(bottom=0.230, top=0.845)
	axis.set_xlim(left=0.173, right=0.788)

	# body
	_draw_triangle(axis, 0.637, 0.921, 0.649, 0.274, 0.188, 0.667, colors[0])
	_draw_triangle(axis, 0.981, 0.769, 0.335, 0.191, 0.393, 0.671, colors[1])
	_draw_triangle(axis, 0.941, 0.397, 0.292, 0.475, 0.456, 0.747, colors[2])
	_draw_triangle(axis, 0.662, 0.119, 0.316, 0.548, 0.662, 0.700, colors[3])
	_draw_triangle(axis, 0.309, 0.081, 0.374, 0.718, 0.681, 0.488, colors[4])
	_draw_triangle(axis, 0.016, 0.626, 0.726, 0.687, 0.522, 0.327, colors[5])
	_draw_text(axis, 0.212, 0.562, labels.get('000001', ''))
	_draw_text(axis, 0.430, 0.249, labels.get('000010', ''))
	_draw_text(axis, 0.356, 0.444, labels.get('000011', ''))
	_draw_text(axis, 0.609, 0.255, labels.get('000100', ''))
	_draw_text(axis, 0.323, 0.546, labels.get('000101', ''))
	_draw_text(axis, 0.513, 0.316, labels.get('000110', ''))
	_draw_text(axis, 0.523, 0.348, labels.get('000111', ''))
	_draw_text(axis, 0.747, 0.458, labels.get('001000', ''))
	_draw_text(axis, 0.325, 0.492, labels.get('001001', ''))
	_draw_text(axis, 0.670, 0.481, labels.get('001010', ''))
	_draw_text(axis, 0.359, 0.478, labels.get('001011', ''))
	_draw_text(axis, 0.653, 0.444, labels.get('001100', ''))
	_draw_text(axis, 0.344, 0.526, labels.get('001101', ''))
	_draw_text(axis, 0.653, 0.466, labels.get('001110', ''))
	_draw_text(axis, 0.363, 0.503, labels.get('001111', ''))
	_draw_text(axis, 0.750, 0.616, labels.get('010000', ''))
	_draw_text(axis, 0.682, 0.654, labels.get('010001', ''))
	_draw_text(axis, 0.402, 0.310, labels.get('010010', ''))
	_draw_text(axis, 0.392, 0.421, labels.get('010011', ''))
	_draw_text(axis, 0.653, 0.691, labels.get('010100', ''))
	_draw_text(axis, 0.651, 0.644, labels.get('010101', ''))
	_draw_text(axis, 0.490, 0.340, labels.get('010110', ''))
	_draw_text(axis, 0.468, 0.399, labels.get('010111', ''))
	_draw_text(axis, 0.692, 0.545, labels.get('011000', ''))
	_draw_text(axis, 0.666, 0.592, labels.get('011001', ''))
	_draw_text(axis, 0.665, 0.496, labels.get('011010', ''))
	_draw_text(axis, 0.374, 0.470, labels.get('011011', ''))
	_draw_text(axis, 0.653, 0.537, labels.get('011100', ''))
	_draw_text(axis, 0.652, 0.579, labels.get('011101', ''))
	_draw_text(axis, 0.653, 0.488, labels.get('011110', ''))
	_draw_text(axis, 0.389, 0.486, labels.get('011111', ''))
	_draw_text(axis, 0.553, 0.806, labels.get('100000', ''))
	_draw_text(axis, 0.313, 0.604, labels.get('100001', ''))
	_draw_text(axis, 0.388, 0.694, labels.get('100010', ''))
	_draw_text(axis, 0.375, 0.633, labels.get('100011', ''))
	_draw_text(axis, 0.605, 0.359, labels.get('100100', ''))
	_draw_text(axis, 0.334, 0.555, labels.get('100101', ''))
	_draw_text(axis, 0.582, 0.397, labels.get('100110', ''))
	_draw_text(axis, 0.542, 0.372, labels.get('100111', ''))
	_draw_text(axis, 0.468, 0.708, labels.get('101000', ''))
	_draw_text(axis, 0.355, 0.572, labels.get('101001', ''))
	_draw_text(axis, 0.420, 0.679, labels.get('101010', ''))
	_draw_text(axis, 0.375, 0.597, labels.get('101011', ''))
	_draw_text(axis, 0.641, 0.436, labels.get('101100', ''))
	_draw_text(axis, 0.348, 0.538, labels.get('101101', ''))
	_draw_text(axis, 0.635, 0.453, labels.get('101110', ''))
	_draw_text(axis, 0.370, 0.548, labels.get('101111', ''))
	_draw_text(axis, 0.594, 0.689, labels.get('110000', ''))
	_draw_text(axis, 0.579, 0.670, labels.get('110001', ''))
	_draw_text(axis, 0.398, 0.670, labels.get('110010', ''))
	_draw_text(axis, 0.395, 0.653, labels.get('110011', ''))
	_draw_text(axis, 0.633, 0.682, labels.get('110100', ''))
	_draw_text(axis, 0.616, 0.656, labels.get('110101', ''))
	_draw_text(axis, 0.587, 0.427, labels.get('110110', ''))
	_draw_text(axis, 0.526, 0.415, labels.get('110111', ''))
	_draw_text(axis, 0.495, 0.677, labels.get('111000', ''))
	_draw_text(axis, 0.505, 0.648, labels.get('111001', ''))
	_draw_text(axis, 0.428, 0.663, labels.get('111010', ''))
	_draw_text(axis, 0.430, 0.631, labels.get('111011', ''))
	_draw_text(axis, 0.639, 0.524, labels.get('111100', ''))
	_draw_text(axis, 0.591, 0.604, labels.get('111101', ''))
	_draw_text(axis, 0.622, 0.477, labels.get('111110', ''))
	_draw_text(axis, 0.501, 0.523, labels.get('111111', ''))

	# legend
	_draw_text(axis, 0.674, 0.824, names[0], text_colors[0])
	_draw_text(axis, 0.747, 0.751, names[1], text_colors[1])
	_draw_text(axis, 0.739, 0.396, names[2], text_colors[2])
	_draw_text(axis, 0.700, 0.247, names[3], text_colors[3])
	_draw_text(axis, 0.291, 0.255, names[4], text_colors[4])
	_draw_text(axis, 0.203, 0.484, names[5], text_colors[5])
	leg = axis.legend(names, loc='best', fancybox=True)
	leg.get_frame().set_alpha(0.5)

	return fig, axis


def venn(groups, names=None, fill=NUMBER, savename=None, show=False, **options):
	'''
	Create a Venn diagram.

	:param groups: If `names` is None, Dict {str: list}. Otherwise,
		`groups` must be a list of lists. Each sublist should contain
		elements that are compared among all groups to create the Venn
		diagram. The string is the name of that group (if `names` is not
		given)
	:param names: List of strings with the names of the groups. Can be
		omitted, if the groups are a dict that contains the names.
	:param fill: What to display. One or more of LOGIC, NUMBER, PERCENT.
		You may concatenate them with | to use two or all three options.
	:param savename: If given (string or file-like object), save the
		Venn diagram to this filename or file-like object.
	:param show: If True, show the Venn diagram in a dynamic matplotlib window.

	:param colors: List of lists with colors for the Venn elements.
		Colors must be defined like [r, g, b, alpha], where all four
		elements must be floats between 0 and 1.
	:param figsize: Tuple with two numbers giving the size of the figure
	:param dpi: int with dpi

	:returns: A tuple with the matplotlib.Figure object and the matplotlib.Axis object
	'''

	if not 1 < len(groups) < 7:
		raise NotImplementedError(f'Please use two to six groups (2-6). You used {len(groups)} groups.')

	venns = {2: _venn2, 3: _venn3, 4: _venn4, 5: _venn5, 6: _venn6}

	if names is None:
		groups = {name: set(groups[name]) for name in groups}
	else:
		groups = {names[i]: set(groups[i]) for i in range(len(groups))}

	names = sorted(groups.keys())
	#name_ids = {2**i: n for i, n in enumerate(names)}

	to_plot = get_labels([groups[name] for name in names], fill=fill)

	plt.close('all') # Reset any possible earlier figures
	fig, axis = venns[len(groups)](to_plot, names, **options)

	if savename is not None:
		try:
			plt.savefig(savename)
		except TypeError:
			print(f'savename "{savename}" is neither a valid string for '
			'a path nor a file-like object. The Venn diagram is not saved.', file=sys.stderr)

	if show:
		plt.show()

	return fig, axis


if __name__ == '__main__':
	A = 'abcdef'
	B = 'bcdefg'
	C = 'cdefgh'
	D = 'fghi'
	E = 'acf'
	F = 'bfh'

	# Just to check if these raise some exception
	get_labels([range(10), range(5,15), range(3,8)], fill=NUMBER)
	venn([A, B], 'AB')
	venn({'A': A, 'B': B, 'C': C}, fill = LOGIC | NUMBER, show=True)
	venn({'A': A, 'B': B, 'C': C, 'D': D}, fill = NUMBER | PERCENT, dpi=300)
	venn([A, B, C, D, E], names = ['A', 'B', 'C', 'D', 'E'], figsize = (25, 25))
	venn({'A': A, 'B': B, 'C': C, 'D': D, 'E': E, 'F': F})
