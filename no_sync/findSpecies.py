#!/usr/bin/env python3

from taxfinder import TaxFinder

TF = TaxFinder()

#tofind = ['aL', 'aK', 'aJ']		# NamPRT, NNMT, NMNAT3
#nottofind = []					# -
tofind = ['aL', 'aK']			# NamPRT, NNMT
nottofind = ['aJ']				# NMNAT3



with open('attributes.txt', 'r') as f:
	for line in f:
		taxid, attributes = line.split()
		isok = True

		for att in tofind:
			if att not in attributes:
				isok = False
				break

		for att in nottofind:
			if att in attributes:
				isok = False
				break

		if isok:
			info = TF.getTaxInfo(taxid)
			if info['rank'] == 'species':
				print(taxid)

