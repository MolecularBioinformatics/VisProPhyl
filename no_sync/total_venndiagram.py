#!/usr/bin/env python3

import venn
import os

with open('tree_config.txt', 'r') as f:
	codes = []
	names = []
	props = []
	for line in f:
		line = line.split('#')[0].strip()
		if not line or line.startswith('>') or line.startswith('!'):
			continue

		lline = line.split(':')
		codes.append(lline[0])
		names.append(lline[1])
		props.append(set())

if not (2 <= len(names) <= 4):
	raise ValueError('The number of proteins to include is wrong ({})!'.format(len(names)))

with open('attributes.txt', 'r') as f:
	for line in f:
		xname, xcodes = line.rstrip().split('\t')
		for i, c in enumerate(codes):
			if c in xcodes:
				props[i].add(xname)

venn.venn(props, names, fill = venn.NUMBER, save_name = 'test.png', show_plot = False)
