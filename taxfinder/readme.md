TaxFinder
=========

Before the TaxFinder module can be used, a taxonomy database must be created. This is done by running `updateTaxonomy.sh`. The script will take up to an hour and will create about 6 GB data. During creation of the database, up to 20 GB disk space is needed.

Initializing TaxFinder in a script may take several seconds.

```python
from taxfinder import TaxFinder

TF = TaxFinder()

TF.getLineage(9606)
```
