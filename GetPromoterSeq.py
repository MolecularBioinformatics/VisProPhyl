#!/usr/bin/env python3

#import re
import os
import sys
from subprocess import call
from Bio import Entrez, SeqIO

Entrez.email = 'nicolai.von-kuegelgen@student.uni-tuebingen.de'
Entrez.tool = 'PullPromoterSeq-NvK'


class taxids:
    def __init__(self):
        """Provide an instance with all taxids of 'good' genomes, their accession numbers and ftp-download locations"""
        self.ids = set()
        self.map = dict()

        with open('Promoters/refseq_taxids.txt', 'r') as f:
            next(f)
            for line in f:
                tax, acc, ftp = line.rstrip().split('\t')
                self.ids.add(tax)
                self.map[tax] = (acc, ftp)


def initFiles():
    """Creates necessary files for the script to run.

    :creates: Promoters/refseq_taxids.txt"""

    print("Creating initial files")

    os.makedirs('Promoters')

    if not os.path.isfile('Promoters/assembly_summary_refseq.txt'):
        call('rsync --copy-links --times rsync://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_refseq.txt .')

    with open('Promoters/assembly_summary_refseq.txt', 'r') as f, open('Promoters/refseq_taxids.txt', 'w') as out:
        next(f)
        next(f)
        out.write('#Taxid\tGenome-Accseion\tftp-link\n')
        for line in f:
            lline = line.rstrip().split('\t')
            # lline[4] == 'na' - refseq categeory unspecified (= not a representative/reference genome)
            if lline[11] in ('Scaffold', 'Contig'): #completeness of genome at least on one chromosome
                continue

            out.write("\t".join((lline[5],lline[0],lline[19]))+'\n') # taxid, acc, ftp-link


#  If this script is to be used without the phylogenetics.py workflow (up to step 2) then this function will have to be
#  substituted with/adapted to a parser for a list of protein Accession-IDs and probably also the TaxFinder module/class
def getUniqueAccs(protein, taxid, cutoff=50):
    """Read combinedtables from Phylogenetics workflow.
     Extract Accesion numbers and taxids for hits, apply refseq as filter, save to file

     :param protein: a protein name from the phylogenetics workflow
     :param taxid: an instance of the taxids class
     :param cutoff: worst evalue (e^-{cutoff} accepted for hits
     :creates: 'Promoters/hits_{prortein}.txt'"""

    print('Getting usable blast hits for {:<10}'.format(protein), end='\r')

    with open("combinedtables/{}.tsv".format(protein), "r") as f:
        next(f)
        accs = set()
        mapping = dict()
        for line in f:
            lline = line.split('\t')
            tax = lline[0]
            acc = lline[1]
            evalue = lline[4]

            #Only work with hits, that have a (good) genome in the refseq db
            if tax not in taxid.ids:
                continue

            #Only accept refseq protein hits, since others will not be found in the refseq genomes anyway (probably? not 100% tested)
            if acc[2] != '_':
                continue

            #reduce hits by evalue
            if evalue == '0.0':
                evalue = 201
            elif 'e-' in evalue:
                evalue = int(evalue.split('e-')[1])
            else:
                evalue = 0

            if evalue < cutoff:
                continue

            accs.add(acc)
            mapping[acc] = tax

    with open('Promoters/hits_{}.tsv'.format(protein), 'w') as out:
        for acc in accs:
            out.write(acc+'\t'+mapping[acc]+'\n')



def loadGenomeFeatures(protein, taxid):
    """Download genome annotation tables for all hits of given protein from NCBI ftp-server

    :uses: Promoters/hits_{protein}.tsv
    :param protein: a protein name from the phylogenetics workflow
    :param taxid: an instance of the taxids class
    :creates: Promoters/featuretables/*.txt"""

    print('Downloading Genome annotation tables for hits of {:<10}'.format(protein), end='\r')

    accs = set()
    os.makedirs('Promoters/featuretables', exist_ok=True)

    with open('Promoters/hits_{}.tsv'.format(protein), 'r') as f:
        for line in f:
            acc, tax = line.rstrip().split('\t')
            accs.add((acc, tax))

    cmd = 'rsync --copy-links --times rsync{}_feature_table.txt.gz featuretables/'

    for acc, tax in accs:
        if os.path.isfile('Promoters/featuretables/{}_feature_table.txt'.format(tax)):
            continue

        print('Downloading Genome annotation table for species: {:<10}'.format(protein), end='\r')

        ftp = taxid.map[tax][1]
        path = ftp.split('/')[-1]
        ret = call(cmd.format(ftp[3:]+'/'+path).split()) #downlaod from ncbi, this always produces text output
        ret += call('gzip -d featuretables/{}_feature_table.txt.gz'.format(path).split()) #unpack
        ret += call('mv featuretables/{}_feature_table.txt featuretables/{}_feature_table.txt'.format(path, tax).split()) # rename
        if ret:
            print(ret) #should print errors


# Searching through genome annotation files may be faster with reg-ex (probably only when applying regex to whole file)
def extractPromoterPos(protein):
    """Extract start & strand of Promoter for all (usable) hits a of protein
    For each write Prot-Acc, TaxID, Genome-Acc, StartPos, Strand, GeneSym, Prot-Name into File

    :uses: Promoters/hits_{protein}.tsv, Promoters/featuretables/{tax}_feature_table.txt
    :param protein: a protein name from the phylogenetics workflow
    :creates: Promoters/{protein}_hitfeatures.tsv"""

    #sort hits by species, so feature tables have to be accessed only once per species
    tax_accs = dict()
    with open('Promoters/hits_{}.tsv'.format(protein), 'r') as f:
        for line in f:
            acc, tax = line.rstrip().split('\t')
            if tax in tax_accs:
                tax_accs[tax] += [acc]
            else:
                tax_accs[tax] = [acc]

    with open('Promoters/{}_hitfeatures.tsv'.format(protein), 'w') as out:
        out.write('#Prot-Acc\tTaxID\tGenome-Acc\tStartPos\tEndPos\tStrand\tGeneSym\tProt-Name\n')
        for tax, accs in tax_accs.items():

            print('Extracting promoter sequences for hits of {}, species: {:<10}'.format(protein, tax), end='\r')
            with open('Promoters/featuretables/{}_feature_table.txt'.format(tax), 'r') as f:
                for line in f:
                    #safe some time / only look for genes
                    if not line.startswith('CDS'):
                        continue

                    #without splitting first might also find 'related_accession', maybe not though after sorting for CDS
                    lline = line.rstrip().split('\t')
    #Columns in *_feature_table.txt
    #feature0	class1	assembly2	assembly_unit3	seq_type4	chromosome5	genomic_accession6	start7	end8	strand9
    #product_accession10	non-redundant_refseq11	related_accession12	name13	symbol14	GeneID15	locus_tag16
                    if lline[10] not in accs:
                        continue

                    accs.remove(lline[10])
                    outl = [lline[10], tax, lline[6], lline[7], lline[8], lline[9], lline[14], lline[13]]
                    out.write('\t'.join(outl)+'\n')

                    if len(accs) ==0:
                        break
                else:
                    print('Extracting promoter sequences for hits of {}, species: {}, features not found: {}'.format(protein, tax, accs))


def loadPromoterSeq(protein, lenght=1500):
    """Using the information of {protein}_hitfeatures.tsv download the associated promoter sequence from NCBI

    :uses: Promoters/{protein}_hitfeatures.tsv
    :param protein: a protein name from the phylogenetics workflow
    :param lenght: The number of NT in front of CDS start to pull
    :creates: Promoters/PromoterSeqs_{Protein}/*.fasta
    """

    entries = set()
    os.makedirs('Promoters/PromoterSeqs_{}'.format(protein),exist_ok=True)

    with open('Promoters/{}_hitfeatures.tsv'.format(protein), 'r') as f:
        next(f)
        for i, line in enumerate(f):
            lline = line.rstrip().split('\t')
            tax, genome_acc, start, end, strand = lline[1:6]

            #Multiple proteins may come from one gene
            entry = (tax, genome_acc, start, end)
            if entry in entries:
                continue
            else:
                entries.add(entry)

            #Dont load same sequence again - mostly for test use
            #if os.path.isfile('PromoterSeqs_{}/{}.fasta'.format(protein, lline[6]+'_'+tax)):
                #continue

            print('Downloading promoter sequences for hits of {}. Found: {} ({} multiple){:<8}'.format(protein, len(entries), i-len(entries), ''), end='\r')

            #start from the *_hitfeatures file is the CDS start, we want region downstream of that
            if strand == '+':
                end = int(start)
                start = int(start) - lenght
                strand = 1
            #If the gene is on the reverse strand, the Promoter is upstream of the gene
            else:
                start = int(end)
                end = int(end) + lenght
                strand = 2


            handle = Entrez.efetch(db="nucleotide",
                                   id=genome_acc,
                                   rettype="fasta",
                                   strand=strand,
                                   seq_start=start,
                                   seq_stop=end)
            record = SeqIO.read(handle, "fasta")
            record.id = 'ref|{}|{}-{}'.format(genome_acc, lline[3], lline[4])
            name = lline[7].replace(':', '').replace('(', '').replace(')', '')
            record.description = "promoter of "+name+"|protein "+lline[0]+"|tax-id "+tax

            with open('Promoters/PromoterSeqs_{}/{}.fasta'.format(protein, lline[6]+'_'+tax), 'w') as out:
                SeqIO.write(record, out, 'fasta')




if __name__ == '__main__':
    import argparse
    workflow = [getUniqueAccs, loadGenomeFeatures, extractPromoterPos, loadPromoterSeq]

    parser = argparse.ArgumentParser(
        description='''This module allows automatic download of promoter sequences starting from protein accession numbers.\n
                    Any number of protein names with blast results from the phylogenetics.py workflow may be used as arguments''')

    parser.add_argument('-i', '--init', action='store_true', help='Force new initialisation of helper files.')
    parser.add_argument('protein', nargs='+', help='Name(s) of protein(s) as used in the phylogenetics.py workflow.')
    #Possibly add support for giving a filename (with protein names) instead of the protein names themselves

    args = parser.parse_args()
    proteins = args.protein

    if args.init:
        initFiles()
        sys.exit()

    if not os.path.isfile('Promoters/refseq_taxids.txt'):
        initFiles()

    TX = taxids()

    for i, func in enumerate(workflow):
        print('')
        for p in proteins:
            if i < 2:
                func(p, TX)
            else:
                func(p)
    print('')