Phylogenetics
=============

This repo bundles all phylogenetic scripts used during my PhD. Details for the scripts are given below. I propose linking the files to the Python path, your PATH, or alias the scripts.

All scripts can be run with `-h` or `--help` to get instructions on how to use them.

Installation dependencies are given separatly with every script.

The scripts were originally written by Mathias Bockwoldt.


phylogenetics.py
----------------

To begin, run `phylogenetics.py --init` to initialise the folder structure. Create one fasta file in `fastas/` for each sequence you want to have analysed. I suggest to run blast on some server (e.g. Stallo). Update `nr` to the latest version and run blast with the following script (adjust accordingly for multiple files):

```bash
db=/path/to/database
fn=/path/to/query_filename.fasta
threads=16	# probably adjust this value depending on the server

blastp -db $db -query $fn -out ${fn%.fasta}.xml -outfmt 5 \
    -max_target_seqs 20000 -evalue 1 -num_threads $threads
```

Move the resulting xml files to the folder `blastresults/`.

Run `phylogenetics.py --all` and examine the results. Especially `proteinlist.txt` must be changed according to your proteins.
