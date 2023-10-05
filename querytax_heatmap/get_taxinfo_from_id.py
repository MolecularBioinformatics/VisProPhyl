#!/usr/bin/env python3
"""
get_taxinfo_from_id.py

Queries NCBI via Entrez to get taxonomic names
and full lineage information according to accession IDs.

To be used to generate heatmap config files
for phyo analysis with randomised organism selection.
"""

from Bio import Entrez

## set your Entrez email
Entrez.email = "hsieh.y.chen@uit.no"


def query_ncbi(taxid, database='taxonomy'):
    """
    Submits query via Entrez to NCBI
    Returns record storing results
    
    :param taxid: int
    :param database: str

    :returns: result record
    """
    handle = Entrez.efetch(db=database, id=str(taxid), retmode='xml')
    record = Entrez.read(handle)
    handle.close()

    if not record:
        raise ValueError(f'Could not find info for taxid {taxid}.')
    return record


def get_str_details(entrezrec_obj, tag='ScientificName'):
    """
    Returns specified description of organism
    from Entrez record object

    tag = 'ScientificName', 'Lineage', 'TaxId'

    :param entrezrec_obj: Bio.Entrez.Parser.ListElement
    :returns: Bio.Entrez.Parser.StringElement
    """
    attr = entrezrec_obj[0][tag]
    return attr
