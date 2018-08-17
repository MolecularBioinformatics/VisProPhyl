#!/usr/bin/env bash

if [ ! -d "$TFPATH" ]
then
	echo 'TFPATH is not set or does not point to an existing path.' >&2
	echo 'Please see the readme to see what this means.' >&2
	exit 1
fi

cd ${TFPATH}

echo "Downloading accession to taxid file..."
curl -R --retry 1 -o prot.accession2taxid.gz ftp://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/prot.accession2taxid.gz

echo "Downloading PDB to taxid file..."
curl -R --retry 1 -o pdb_chain_taxonomy.tsv.gz ftp://ftp.ebi.ac.uk/pub/databases/msd/sifts/flatfiles/tsv/pdb_chain_taxonomy.tsv.gz

echo "Unzipping accession to taxid file..."
gunzip prot.accession2taxid.gz

echo "Unzipping PDB to taxid file..."
gunzip pdb_chain_taxonomy.tsv.gz

echo "Running Python to create the acc2taxid database..."
python3 <<EOF
with open('unsorted_acc2taxid', 'w') as out:
	print('py: reading pdb_chain_taxonomy.tsv and writing unsorted_acc2taxid...')
	with open('pdb_chain_taxonomy.tsv', 'r') as f:
		next(f)
		next(f)
		pdb = None
		taxid = None
		for line in f:
			lline = line.split()
			if pdb == lline[0]:
				continue
			if pdb is not None:
				out.write('{:<12}{:<7}\n'.format(pdb.upper(), taxid))
			pdb = lline[0]
			taxid = lline[2]
	print('py: reading prot.accession2taxid and writing unsorted_acc2taxid...')
	with open('prot.accession2taxid', 'r') as f:
		next(f)
		for line in f:
			lline = line.split()
			# Writes only the accession and the taxid to the new file.
			# Accession and taxid are right-padded with spaces:
			# acc12345____1234___{NEWLINE}
			# ^--- 12 ---^^- 7 -^
			# This allows later for an easier binary search in the file.
			# A field separator (like tab) is not necessary as the column width is constant
			out.write('{:<12}{:<7}\n'.format(lline[0], lline[2]))
EOF

echo "Sorting accession to taxid file..."
sort -s -k 1,1 --parallel=6 -o acc2taxid unsorted_acc2taxid

echo "Counting lines of the acc2taxid file..."
wc -l < acc2taxid > numLines

echo "Cleaning up the acc2taxid part..."
rm prot.accession2taxid
rm pdb_chain_taxonomy.tsv
rm unsorted_acc2taxid

echo "Downloading lineage file..."
curl -R --retry 1 -o taxdump.tar.gz ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz

echo "Unzipping lineage file..."
tar xzf taxdump.tar.gz
rm taxdump.tar.gz

echo "Running Python to create the taxonomy database..."
python3 <<EOF
data = {}	# data[taxid] = [depth, parent, rank, name]

print('py: reading names.dmp...')

with open('names.dmp', 'r') as f:
	for line in f:
		lline = line.split('\t|\t')
		# taxid, name, uniqueName, nameClass
		if lline[3] == 'scientific name\t|\n':
			data[lline[0]] = ['@', '@', '@', lline[1]]

print('py: reading nodes.dmp...')

with open('nodes.dmp', 'r') as f:
	for line in f:
		lline = line.split('\t|\t')
		# taxid, parentTaxid, rank, *others
		data[lline[0]][1] = lline[1]
		data[lline[0]][2] = lline[2]

print('py: writing taxinfo...')

with open('taxinfo', 'w') as out:
	for taxid in sorted(data.keys()):
		# TaxID, Level, Parent, Rank, Name
		level = 0
		tid = taxid
		while tid != '1':
			tid = data[tid][1]
			level += 1
		data[taxid][0] = str(level)
		out.write('{}\t{}\n'.format(taxid, '\t'.join(data[taxid])))
EOF

echo "Cleaning up..."
rm *.dmp
rm gc.prt
rm readme.txt

echo "Databases are created. Now it's time for some holidays..."
