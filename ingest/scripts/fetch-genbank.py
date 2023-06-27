#!/usr/bin/env python3

import json
import argparse
from Bio import SeqIO, Entrez

# to use the efetch API, the docs indicate only around 200 IDs should be provided per request
# https://www.ncbi.nlm.nih.gov/books/NBK25499/#chapter4.EFetch
BATCH_SIZE = 200

# As of mid 2023 there are ~11k records for HBV
MAX_IDS = 100_100

Entrez.email = "hello@nextstrain.org"

def parse_args():
    parser = argparse.ArgumentParser(description='Annotate sequences using a genbank reference')
    parser.add_argument('--term', help='Genbank search term. Replace spaces with "+"', default='Hepatitis+B+virus[All+Fields]complete+genome[All+Fields]')
    parser.add_argument('--output', required=True, type=str, help='Output file (Genbank)')
    return parser.parse_args()

def get_sequence_ids(term):
    handle = Entrez.esearch(db="nucleotide", term=term, retmode="json", retmax=MAX_IDS)
    data = json.loads(handle.read())
    ids = data['esearchresult']['idlist']
    print(f"Search term '{term}' returned {len(ids)} IDs. Processing in batches of n={BATCH_SIZE}.")
    idx = 0
    while idx < len(ids):
        yield ids[idx:idx+BATCH_SIZE]
        idx+=BATCH_SIZE

def fetch(id_list):
    handle = Entrez.efetch(db="nucleotide", id=','.join(id_list), rettype="gb", retmode="text")
    return SeqIO.parse(handle, "genbank")

if __name__=="__main__":
    args = parse_args()
    
    genomes = []
    for id_list in get_sequence_ids(args.term):
        genomes.extend(fetch(id_list))

    SeqIO.write(genomes, args.output, "genbank")
