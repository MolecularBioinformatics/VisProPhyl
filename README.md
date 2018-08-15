Phylogenetics
=============

This repo bundles all phylogenetic scripts used during my PhD. Details for the scripts are given below. I propose linking the files to the Python path, your PATH, or alias the scripts.

All scripts can be run with `-h` or `--help` to get instructions on how to use them.

Installation dependencies are given separately with every script.

The scripts were originally written by Mathias Bockwoldt.

```sh
pip install -r requirements.txt
pip install .
```

Todo
----

* phylogenetics/lineage_values* → have to be incorporated
* taxfinder → Here I have to think about how to include the update scripts
* bin/msa_viewer → How much work should I use here?
* no_sync/* → What to do with these random scripts and also other random script that fly around?
* Thorough explanation of every script and module. Why is it there? What does it do? How to use?
    * phylogenetics
    * phylotree
    * msa_blast
    * lineage_values*
    * msa_viewer?
    * taxfinder/...
    * phylogenetics/phylogenetics.py
    * phylogenetics/venn.py


phylogenetics.py
----------------

This script does many phylogenetic analyses based on Blast results. It also creates the basics for `phylotree.py`.

To begin, run `phylogenetics.py --init` to initialise the folder structure and templates. Create one fasta file in `fastas/` for each sequence you want to have analysed. I suggest to run blast on some server (e.g. Stallo). Update `nr` to the latest version and run blast with the following script (adjust accordingly for multiple files):

```bash
db=/path/to/database
fn=/path/to/query_filename.fasta
threads=16	# probably adjust this value depending on the server

blastp -db $db -query $fn -out ${fn%.fasta}.xml -outfmt 5 \
    -max_target_seqs 20000 -evalue 1 -num_threads $threads
```

Move the resulting xml files to the folder `blastresults/`. Then have a look at the various `.txt` files. They have instructions written inside them. Here is a short overview of what they are doing:

- *crosshits.txt*: Define which proteins to check for a more detailed analysis of crosshits
- *heatmap_config.txt*: Define the species to show in the heatmap
- *limits.txt*: Define e-value and length limits for proteins to be considered by the analysis
- *proteinlist.txt*: Define the proteins used. Only proteins defined here will be considered. Also, if you used multiple sequences to pose for the same protein (e.g. splice variants or from different species), this has to be defined here
- *tree_config.txt*: Necessary for `phylotree.py`
- *tree_to_prune.txt*: Necessary for `phylotree.py`

After configuring the `.txt` files, run `phylogenetics.py --all` and examine the results.

The outputs are the following.

**TODO**
