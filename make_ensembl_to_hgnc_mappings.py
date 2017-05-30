#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function
from berkley_hash_set import BerkleyHashSet 
from biomart_data import BiomartData 
from biomart_data import BiomartDataset
from tqdm import tqdm 

#BIOMART_URL = 'http://www.ensembl.org/biomart'
BIOMART_URL = 'http://grch37.ensembl.org/biomart'


def make_mappings(output_db_name='hgnc_by_ensembl.db'):

    hgnc_by_ensembl = BerkleyHashSet(output_db_name)
    hgnc_by_ensembl.reset()

    data = BiomartData(
        dataset=BiomartDataset(BIOMART_URL, name='hsapiens_gene_ensembl'),
        attributes=[
            'hgnc_symbol',
            'ensembl_transcript_id'     # one could use ensembl_gene_id instead
        ],
        filters={}
    )

    for entry in tqdm(data, total=data.count()):
        if entry.hgnc_symbol and entry.ensembl_transcript_id:
            hgnc_by_ensembl[entry.ensembl_transcript_id].add(entry.hgnc_symbol) 
    """

    with open('ensembl_to_hgnc.tsv') as f:
        for line in f:
            ensembl, hgnc = line.split()
            hgnc_by_ensembl[ensembl].add(hgnc)

    """
if __name__ == '__main__':
    make_mappings()
