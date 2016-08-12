#!/usr/bin/env python3

import re
import os

class TaxFinder():

	def __init__(self):

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

	def getNameFromID(self, taxid):
		''' Returns the taxonomic name of the given taxid '''

		return self.taxdb[int(taxid)]['name']


	def getTaxInfo(self, taxid):
		'''
		Get taxonomic information for the given taxid.
		:returns: {'taxid': int, 'level': int, 'parent': int, 'rank': str, 'name': str}
		'''

		taxid = int(taxid)

		try:
			taxinfo = self.taxdb[taxid]
			taxinfo['taxid'] = taxid
		except KeyError:
			print('Taxid not found:', taxid)
			taxinfo = {'taxid': 1, 'level': 0, 'parent': 0, 'rank': 'no rank', 'name': 'unclassified'}

		return taxinfo


	def getInfoFromHitDef(self, hitid, hitdef, newHeader = True):
		'''
		Get all taxonomy information from a hit id and hit definition (may include several species)
		:returns: [{'taxid': int, 'level': int, 'parent': int, 'rank': str, 'name': str, 'acc': str, 'protname': str}, ...]
		'''

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
		'''
		Given the taxid of a subspecies, returns the species or raises a ValueError if no species could be found.
		'''

		lineage = self.getLineageFast(int(taxid))
		for tid in lineage[::-1]:
			if self.getTaxInfo(tid)['rank'] == 'species':
				return tid

		raise ValueError('No species found for {}'.format(taxid))


	def getLowestReasonableTaxon(self, taxid):
		'''
		Given a taxid, returns the taxid the closest species or higher level (that is not `no rank`) that does not contain 'sp.' in its name. Raises a ValueError if no "reasonable taxon" could be found.
		'''

		notok = {'no rank', 'subspecies', 'forma', 'varietas'}
		lineage = self.getLineageFast(int(taxid))
		for tid in lineage[::-1]:
			info = self.getTaxInfo(tid)
			rank = info['rank']
			if rank not in notok and 'sp.' not in info['name']:
				return tid

		raise ValueError('No reasonable taxon found for {}'.format(taxid))


	def getLineage(self, taxid, display = 'name'):
		'''
		Given a taxid, returns the lineage up to `root` as tuple. `display` configures how the lineage should be shown. If `display` is 'name', the taxonomic name will be used. If it is 'taxid', the taxid will be used. If it is anything else, name^taxid will be used.
		This method uses caching. If the lineage for a taxid was already found before, it will return that lineage in the `display` mode used in the first search, ignoring the current `display` value.
		If the taxid could not be found, an empty tuple will be returned.
		'''

		taxid = int(taxid)

		if taxid in self.lineageCache:
			return self.lineageCache[taxid]

		lineage = []

		while True:
			try:
				t = self.taxdb[taxid]
			except KeyError:
				self.lineageCache[taxid] = tuple()
				return tuple
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
		'''
		Given a taxid, returns the lineage up to `root` as list. All elements will be taxids.
		This method is faster than `getLineage`, so use this when you need many lineages.
		If the taxid could not be found, an empty tuple will be returned.
		'''

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

		lin = tuple(lineage[::-1])

		self.lineageCache[taxid] = lin

		return lin
