Phylogenetics
=============

So, the basic functions *should* work by now. But there are still some functions that I did not implement, yet. The current status is in `status.csv`.

The script `phylogenetics.py` can be imported as module, `phylotree.py` not (yet) reasonably. All `txt` files and the `html` file are needed for the initialization of new projects.

Both, `phylogenetics.py` and `phylotree.py` can be called directly. Use the `-h` or `--help` flag to get information about the usage.


phylogenetics.py
----------------

To begin, run `phylogenetics.py --init` to initialise the folder structure. Create one fasta file in `fastas/` for each sequence you want to have analysed. I suggest to run blast on some server (e.g. Stallo). Update `nr` to the latest version and run blast with the following script (adjust accordingly for multiple files):

```bash
$db=/path/to/database
$fn=/path/to/query_filename.fasta
$threads=16	# probably adjust this value depending on the server

blastp -db $db -query $fn -out ${fn%.fasta}.xml -outfmt 5 -max_target_seqs 20000 -evalue 1 -num_threads $threads
```

Move the resulting xml files to the folder `blastresults/` (you may need to create it).

Run `phylogenetics.py --all` and examine the results. Especially `proteinlist.txt` must be changed according to your proteins.
