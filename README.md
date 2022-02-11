Phylogenetics
=============

TODO: Description of what this is :)


```sh
$ pip install phylogenetics ### TODO: This is not working for now!
### TODO: Use for now:
$ pip install wheel  # The package ete3 needs wheel to be installed first...
$ cd path/to/Phylogenetics
$ pip install -e .
```


Todo
----

* phylogenetics/lineage_values* → have to be incorporated
* no_sync/* → What to do with these random scripts and also other random script that fly around?
* Thorough explanation of every script and module. Why is it there? What does it do? How to use?
    * phylogenetics
    * phylotree
    * msa_blast
    * lineage_values*
    * phylogenetics/phylogenetics.py
    * phylogenetics/venn.py


phylogenetics.py
----------------

This script does many phylogenetic analyses based on Blast results. It also creates the basics for `phylotree.py`.

To begin, run `phylogenetics.py --init` to initialise the folder structure and templates. Create one fasta file in `fastas/` for each sequence you want to have analysed. I suggest to run blast online on the NCBI homepage.

Move the resulting xml files to the folder `blastresults/`. Then have a look at the various `.txt` files. They have instructions written inside them. Here is a short overview of what they are doing:

- *heatmap_config.txt*: Define the species to show in the heatmap
- *limits.txt*: Define e-value and length limits for proteins to be considered by the analysis
- *proteinlist.txt*: Define the proteins used. Only proteins defined here will be considered. Also, if you used multiple sequences to pose for the same protein (e.g. splice variants or from different species), this has to be defined here
- *tree_config.txt*: Necessary for `phylotree.py`
- *tree_to_prune.txt*: Necessary for `phylotree.py`

After configuring the `.txt` files, run `phylogenetics.py --all` and examine the results.
