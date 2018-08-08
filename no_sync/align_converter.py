#!/usr/bin/env python3

from Bio import AlignIO, Alphabet
from io import StringIO
import sys


def convert(from_format, to_format, in_string_or_handle, seq_count = None, alphabet = None):
	'''
	Takes a string or file handle of a sequence alignment in `from_format` and returns a string with the sequence alignment converted to `to_format`.

	*param from_format* Format of the input alignment (str)
	*param to_format* Desired format of the output (str)
	*param in_string_or_handle* Either the alignment as string, or the file handle of the file with the alignment.
	*seq_count* Optional number of sequences expected in each alignment. Recommended for fasta format files. (int)
	*alphabet* Optional alphabet of the alignment (BioPython Alphabet object)
	'''

	# If the input is a string, it must be turned into a StringIO. Otherwise, Biopython mistakes it for a file name
	if isinstance(in_string_or_handle, str):
		in_string_or_handle = StringIO(in_string_or_handle)

	out_string = StringIO()
	alignments = AlignIO.parse(in_string_or_handle, from_format, seq_count, alphabet)
	AlignIO.write(alignments, out_string, to_format)

	return out_string.getvalue()


READ_FORMATS = {'clustal', 'emboss', 'fasta', 'fasta-m10', 'ig', 'maf', 'mauve', 'nexus', 'phylip', 'phylip-sequential', 'phylip-relaxed', 'stockholm'}
WRITE_FORMATS = {'clustal', 'fasta', 'maf', 'mauve', 'nexus', 'phylip', 'phylip-sequential', 'phylip-relaxed', 'stockholm'}

ALPHABETS = {'dna': Alphabet.DNAAlphabet,
		'rna': Alphabet.RNAAlphabet,
		'nuc': Alphabet.NucleotideAlphabet,
		'prot': Alphabet.ProteinAlphabet,
		'prot3': Alphabet.ThreeLetterProtein}


formats_string = '''Name      read/write Comment
----      ---------- -------
clustal        rw    The format of the Clustal family.
emboss         r_    The EMBOSS simple/pairs alignment format.
fasta          rw    The widely used FASTA format.
fasta-m10      r_    FASTA format with -m10 activated.
ig             r_    The IntelliGenetics file format.
maf            rw    Multiple Alignment Format (MAF); e.g. UCSC genome browser.
mauve          rw    Mauve’s eXtended Multi-FastA (XMFA) file format.
nexus          rw    Also known as PAUP format.
phylip         rw    Strict interlaced PHYLIP: truncates names at 10 chars.
phylip-sequential rw Strict sequential PHYLIP: truncates names at 10 chars.
phylip-relaxed rw    Relaxed PHYLIP: allows long names.
stockholm      rw    Also known as PFAM format; supports rich annotation.

For more information, see http://biopython.org/wiki/AlignIO'''

alphabets_string = '''dna   - DNA
rna   - RNA
nuc   - Nucleotides (DNA or RNA)
prot  - Protein one letter code
prot3 - Protein three letter code'''


if __name__ == '__main__':
	import argparse
	parser = argparse.ArgumentParser(description='This script converts alignment files from one format into another. For a list of supported formats, invoke this script with -l or --list.')

	parser.add_argument('from_format', help='The format of the input file')
	parser.add_argument('to_format', help='The desired format of the output file')
	parser.add_argument('in_file', help='Filename of the input file. For stdin, use -')
	parser.add_argument('-l', '--list', action='store_true', help='Shows a list of possible formats and alphabets and exits the script')
	parser.add_argument('-a', '--alphabet', help='Alphabet of the alignment. Possible values: (nuc, dna, rna, prot, prot3)')

	args = parser.parse_args()

	if args.list:
		print(formats_string)
		print(alphabets_string)

	if args.from_format not in READ_FORMATS:
		print('The format "{}" cannot be read. This is the list of possible formats:'.format(args.from_format))
		print(formats_string)
		sys.exit()

	if args.to_format not in WRITE_FORMATS:
		print('The format "{}" cannot be written. This is the list of possible formats:'.format(args.from_format))
		print(formats_string)
		sys.exit()

	alphabet = None
	if args.alphabet:
		if args.alphabet in ALPHABETS:
			alphabet = ALPHABETS[args.alphabet]
		else:
			print('The alphabet "{}" is unknown. This is the list of possible alphabets:'.format(args.alphabet))
			print(alphabets_string)
			sys.exit()

	if args.in_file == '-':
		alignment = sys.stdin.read()
	else:
		try:
			alignment = open(args.in_file).read()
		except IOError:
			print('The file "{}" was not found.'.format(args.in_file))
			sys.exit()

	seq_count = None
	if args.from_format == 'fasta':
		seq_cound = alignment.count('>')

	print(convert(args.from_format, args.to_format, alignment, seq_count = seq_count, alphabet = alphabet))



