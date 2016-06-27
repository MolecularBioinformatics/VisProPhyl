#!/usr/bin/env python3

import struct
import re
import os
import sys

class TaxFinder():

	def __init__(self):

		#self.binfile = open('/home/mathias/bin/KronaTools-2.6/taxonomy/gi_taxid.dat', 'rb')
		# The database is saved in binary. Each gi number (starting with 0) has a 4 byte taxonomy id associated with it.
		# Empty gi numbers are filled with 0 (zero) taking also 4 bytes. By this, the db access is quite easy (see self.getTaxID()).
		########################################


		self.acc2taxid = open(self._getFN('acc2taxid'), 'rb')
		with open(self._getFN('numLines') ,'r') as f:
			self.numLines = int(f.read().rstrip())

		self.taxdb = {}

		with open(self._getFN('taxinfo'), 'r') as namefile:
			for line in namefile:
				# TaxID, Level, Parent, Rank, Name
				l = line.split('\t')
				self.taxdb[int(l[0])] = {'level': int(l[1]), 'parent': int(l[2]), 'rank': l[3], 'name': l[4].rstrip()}

		self.lineageCache = {}
		self.taxidCache = {}


	def __enter__(self):
		return self


	def __exit__(self, exc_type, exc_value, traceback):
		self.acc2taxid.close()
		#self.binfile.close()	#############


	def _getFN(self, fn):
		''' Gets absolute path for a given file that is in the same directory as this script '''

		return os.path.join(os.path.dirname(os.path.realpath(__file__)), fn)


	def getTaxID(self, acc):
		''' Finds the NCBI taxonomy id given an accession id '''
		# Accessions are always uppercase
		acc = acc.upper()

		# Only consider the accession without the version part
		if '.' in acc:
			acc = acc.split('.')[0]

		# If we already looked for the accesion, get it from the cache
		if acc in self.taxidCache:
			return self.taxidCache[acc]

		lo = 0
		hi = self.numLines
		x = acc.encode('utf-8')	# Turns the accession id into a bytestring

		# Simple binary search in the sorted file with constant line length
		while lo < hi:
			mid = (lo + hi) >> 1
			self.acc2taxid.seek(mid*20)
			a = self.acc2taxid.read(12)
			if x <= a:
				hi = mid
			else:
				lo = mid + 1

		self.acc2taxid.seek(lo*20)
		rawread = self.acc2taxid.read(19).decode('utf-8')
		testacc = rawread[:12].rstrip()

		if testacc != acc:
			taxid = 1
		else:
			taxid = int(rawread[12:].rstrip())

		self.taxidCache[acc] = taxid

		return taxid


	def OLD_getTaxID(self, gi):
		''' Finds the NCBI taxonomy id given a gi number '''
		self.binfile.seek(int(gi)*4, 0)		# Go to position gi*4
		taxid = self.binfile.read(4)		# Read 4 bytes from that position
		try:
			return struct.unpack('I', taxid)[0]	# Turn these 4 bytes into an unsigned 32 bit integer (I)
		except struct.error:
			print('GI: {}, Bytes: {}'.format(gi, taxid))
			return 1


	def getNameFromID(self, taxid):
		''' Returns the name fitting to the taxid '''

		return self.taxdb[taxid]['name']


	def getTaxInfo(self, taxid):
		'''
		Get taxonomy information given a taxid number.
		:returns: dict('taxid': int, 'level': int, 'parent': int, 'rank': str, 'name': str)
		'''

		try:
			taxinfo = self.taxdb[taxid]
			taxinfo['taxid'] = taxid
		except KeyError:
			print('Taxid not found:', taxid)
			taxinfo = {'taxid': 1, 'level': 0, 'parent': 0, 'rank': 'no rank', 'name': 'unclassified'}

		return taxinfo


	def getInfoFromHitDef(self, hitid, hitdef, newHeader = True):
		''' Get all taxonomy information from a hit id and hit definition (may include several species) '''

		#if newHeader:
		#	hit = hitid + '|' + hitdef
		#	reResults = re.finditer('^[^\|]+\|[^\|]+\|([0-9]+)\|[^ ]+   [^ ]+   [^ ]+   ([^\[]+).*$', hit)	# Group 1 is gi number, group 2 is protein name
		#else:
		#	hit = hitid + hitdef
		#	reResults = re.finditer('gi\|([0-9]+)\|[^\|]+\|[^\|]+\|([^\[]+)', hit)	# Group 1 is gi number, group 2 is protein name

		hit = hitid + hitdef
		reResults = re.finditer('gi\|[0-9]+\|[^\|]+\|([^\|]+)\|([^>]+)', hit)	# group 1 is accession, group 2 is protein name

		results = []

		for r in reResults:
			acc = r.group(1).strip()
			protname = r.group(2).strip()

			if '[' in protname:
				protname = protname.split('[')[0].rstrip()

			res = self.getTaxInfo(self.getTaxID(acc))
			res['acc'] = acc
			res['protname'] = protname

			results.append(res)

		return results


	def getSpeciesFromSubspecies(self, taxid):
		lineage = self.getLineageFast(int(taxid))
		for tid in lineage[::-1]:
			if self.getTaxInfo(tid)['rank'] == 'species':
				return tid

		return None


	def getLineage(self, taxid, display = 'name'):
		if taxid in self.lineageCache:
			return self.lineageCache[taxid]

		lineage = []

		while True:
			try:
				t = self.taxdb[taxid]
			except KeyError:
				self.lineageCache[taxid] = tuple()
				return []
			if display == 'taxid':
				s = taxid
			elif display == 'name':
				s = t['name']
			else:
				s = t['name'] + '^' + str(taxid)
			lineage.append(str(s))
			if taxid == 1:
				break
			taxid = t['parent']

		self.lineageCache[taxid] = tuple(lineage[::-1])

		return self.lineageCache[taxid]


	def getLineageFast(self, taxid):
		if taxid in self.lineageCache:
			return self.lineageCache[taxid]

		lineage = []

		while taxid != 1:
			try:
				t = self.taxdb[taxid]
			except KeyError:
				self.lineageCache[taxid] = tuple()
				return tuple()
			lineage.append(taxid)
			taxid = t['parent']

		lineage.append(1)

		self.lineageCache[taxid] = tuple(lineage[::-1])

		return self.lineageCache[taxid]
