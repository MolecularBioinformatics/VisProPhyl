VisProPhyl
=============

This package contains tools to map BLAST results on the NCBI taxonomy phylogenetic tree. It helps to analyse the presence/absence of proteins or genes in various taxa. This information can be useful, for example, for the analysis of two competing pathways (see e.g. [Bockwoldt et al. 2019](https://doi.org/10.1073/pnas.1902346116)).

If you used this tool in your work, please cite:

> TODO: Citation of our paper.


Installation
------------

You will need [Python 3.6+](https://www.python.org/) and, optionally, [BLAST+](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download). Everything else can be install using `pip`. We suggest to use a [virtual environment](https://docs.python.org/3/library/venv.html) to install everything into.

```sh
# Create a virtual environment. Not strictly necessary,
# but helps fiddling with different versions.
python3 -m venv phyloenv

# Activate the environment. You have to do this every time, you want
# to work with VisProPhyl and open a new shell.
source /path/to/phylotest/phyloenv/bin/activate

# Install VisProPhyl
pip install visprophyl

# You need to tell taxfinder once to download the latest NCBI
# taxonomy information
taxfinder_update

# Now, you can run all the good stuff. Start, for example, with:
phylogenetics --init
```

After the installation, when opening a new terminal, you only have to reactivate the virtual environment

```sh
source /path/to/phylotest/phyloenv/bin/activate
```

Components
----------

This package consists of four command line tools and one importable Python module. All command line tools can be run with `--help` to get help.


Workflow
--------

A typical workflow is described in (TODO: Our paper) and can be summed up to these steps:

1. Create a folder for your project. Open a terminal/shell in this folder and run `phylogenetics --init`.
2. Collect fasta files of your proteins or genes of interest. These are called seeds. Fasta files of multiple species for the same protein or gene are possible. Place the fasta files in the `fasta` folder of your project. Modify the files `proteinlist.txt` and `limits.txt`.
3. Run BLAST against the database of your choice (e.g. `nr` for proteins, `nt` for genes). This can be done using the command line tool or the NCBI BLAST website. It is important that taxonomy ids are included in the results. This is given in all NCBI databases. For your own databases, please consult the documentation for [`makeblastdb`](https://www.ncbi.nlm.nih.gov/books/NBK569841/) on how to include taxonomy ids. If you run BLAST on the NCBI Website, make sure to download the result as "Single-file XML2" and save them to the `blastresults` folder with the same filename as the corresponding fasta file but with `.xml` instead of `.fasta`.
4. Run the command `phylogenetics --all`.
5. Modify the files `tree_config.txt` and `tree_to_prune.txt`.
6. Run the command `phylotree --show`. You may want to go forth and back between modifying `tree_to_prune.txt` and visualizing the tree until the focus of the tree is right.
7. Save the tree using `phylotree --outfile tree_name.pdf`.
8. If you are interested in the interactive heatmap, modify `heatmap_config.txt` and run `phylogenetics --only intheat`.


Output files
------------

`histograms/` shows the length distribution of BLAST results for each seed. `blastmappings/` shows where the seeds map on the BLAST results dependent on the e-value cutoff. Both can be used to tweak `limits.txt`.

`trees/` contains the trees for each seed. These trees contain all species in which the seed was found by BLAST. `general.tre` is the combined tree for all seeds.

`heatmap.html` is an interactive heatmap of selected taxa. You can select taxa and other parameters in `heatmap_config.txt` and re-run `phylogenetics --only intheat`.

`matrix.csv` is a table that shows how similar two seeds are. The number of the cell where seed A and seed B cross is the e-value at which the BLAST search of seed A found seed B.


Other command line tools
------------------------

`blast2fasta` can be used to download sequences based on BLAST results.


Importable modules
------------------

`phylogenetics` can not only used as command line tool, but also imported for your own workflows without the hardcoded filenames etc.

`venn` can create Venn diagrams. There is no direct integration with the other tools in this package, but it may serve useful for your own custom workflows.

`sample_taxids` can be used to randomly sample taxids from a file of taxids, and write these out to a heatmap configuration file. Rerunning the heatmap step of phylogenetics will then redraw the heatmap with sampled taxids.


Utilities
---------

There are some utility scripts that are only on Github and not downloaded by `pip`. You can find them in the `utils` folder in the repository. These scripts are not polished and are meant to be changed before using them. They might come in handy, though, so we did not delete them outright. Each of the scripts has a little information about them in the top of the file.

`lineage_value.py` gives an overview about how present given seeds are along a phylogenetic lineage. If you, for example, given human (taxonomy id 9606) as target, it will show, how well the seeds are found in Hominidae, Simiiformes, Primates, Mammalia, etc.

`static_heatmap.py` is similar to the `intheat` step in `phylogenetics`. As the name says, it is not interactive, but can be used as a possible starting image for publication.
