#!/usr/bin/env python3
"""
sample_taxids.py

Reads in list of NCBI taxonomic ids and a sample size.
Returns randomly sampled taxids.

Used to randomly sample orgs for phylo heatmap.
"""

from pathlib import Path
import argparse
from random import sample

from get_taxinfo_from_id import *

COLORSDICT = {'C0':'00cc00',
             'C1':'cc0000',
             'C2':'00cccc',
             'C3':'0000cc',
             'C4':'cc00cc',
             'C5':'00cc00',
             'C6':'cc0000',
             'C7':'00cccc',
             'C8':'0000cc',
             'C9':'cc00cc'}

def read_in_taxids(taxidfile):
    """
    Reads in file with taxids
    Returns list of taxids

    :param taxidfile: pathlib.PosixPath
    :returns: list
    """
    idlist = []
    with open(taxidfile, 'r') as tfile:
        for line in tfile.readlines():
            idlist.append(line.strip())
    return list(set(idlist))


def rand_sample_taxids(listoftaxids, samplesize=50):
    """
    Randomly samples list of taxids 
    Returns sampled list of a certain size

    :param listoftaxids: list
    :returns: list
    """
    return sample(listoftaxids, samplesize)


def get_taxa_name_dict(sampled_taxids):
    """
    Takes in list of sampled taxids, queries NCBI
    via Entrez, stores name results in a dict of form:
    {taxid: [sciname, commonname, lineage]}

    :param sampled_taxids: list
    :returns: dict
    """
    taxnamedict = {}
    for taxid in sampled_taxids:
        try:
            rec = query_ncbi(taxid)
        except ValueError as valerr:
            print(valerr)
            continue
        try:
            sciname = get_str_details(rec, 'ScientificName')
        except KeyError as keyerr:
            sciname = ''
        try:
            commname = get_str_details(rec, 'OtherNames')['GenbankCommonName']
        except KeyError as keyerr:
            commname = ''
        try:
            lineage = get_str_details(rec, 'Lineage')
        except KeyError as keyerr:
            lineage = ''
        taxnamedict[taxid] = [sciname, commname, lineage]
    return taxnamedict


def format_dict_to_lines(sometaxname_dict):
    """
    Takes in dict of taxonomic info.
    Returns a list of string of names of format

    'sci_name(common_name)^taxid'

    :param sometaxname_dict: dict
    :returns: list
    """
    allnamelist = []
    for key, val in sometaxname_dict.items():
        sciname = '_'.join(val[0].split())
        commname = '_'.join(val[1].split())
        if not commname:
            kingdom = val[2].split('; ')[1]
            commname = kingdom
        allnamelist.append(f'{sciname}({commname})^{key}')
    return '\n'.join(allnamelist)


def write_names_heatmap_config(pathtoconfig, taxnamelist_str):
    """
    Writes out sampled taxonomic ids and names to 
    heatmap_config.txt

    :param pathtoconfig: pathlib.PosixPath
    :param taxnamelist_str: list
    """
    with open(pathtoconfig, 'w+') as configf:
        configf.write('TAXA_TO_SHOW\n')
        configf.write(taxnamelist_str)
        configf.write('\n')
        configf.write('COLORS\n')
        for key,val in COLORSDICT.items():
            configf.write(f'{key}\t{val}\n')
        configf.write('\n')
        configf.write('ALGO\n')
        configf.write('centroid')


def sample_taxids(taxidfile, pathtoconfigfile, samplesize=50):
    """
    Reads in list of taxids. Queries them on NCBI
    via Entrez. 
    Takes a random sample of ids based on samplesize.
    Writes out new heatmap config file for step 9 of phylo analysis.

    :param taxidfile: pathlib.PosixPath
    :param pathtoconfigfile: pathlib.PosixPath
    :param samplesize: int
    """ 
    full_id_list = read_in_taxids(taxidfile)
    sampled_id_list = rand_sample_taxids(full_id_list, samplesize)

    dict_taxnames = get_taxa_name_dict(sampled_id_list)
    write_names_heatmap_config(pathtoconfigfile, format_dict_to_lines(dict_taxnames))
