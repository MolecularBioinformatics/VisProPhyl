TaxFinder
=========

Before the TaxFinder module can be used, a taxonomy database must be created and the script must be told, where the database lives. First, set the environmetal variable `TFPATH` to the path, where the database lives (or shall live). The directory must exist (to avoid creating directories due to a typo in the path) and must be readable by the TaxFinder and writable by `updateTaxonomy.sh`. For example:

```bash
# in ~/.bashrc or similar
TFPATH=~/work/taxfinder_database
```

Then run the script `updateTaxonomy.sh`. But be aware that the script will take up to two hours and will create about 10 GB data. During creation of the database, up to 30 GB disk space is needed.

After these preparations, TaxFinder can be used. Initialization may take several seconds, as the database has to be read.

```python
from taxfinder import TaxFinder

TF = TaxFinder()

print(TF.getLineage(9606))
```
