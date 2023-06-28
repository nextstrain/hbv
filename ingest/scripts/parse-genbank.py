#!/usr/bin/env python3

"""
Most code was originally written by @kistlerk as part of `curate_genbank_to_fasta.ipynb`
[blab/adaptive-evolution](https://github.com/blab/adaptive-evolution)
"""

import argparse
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import dateutil
import dateutil.parser
from collections import defaultdict
import re
import csv

class MissingDate(Exception):
    pass

def parse_genotype(genbank_note):
    """
    Returns (genotype, subgenotype), each of which may be 'None' (string)
    """
    #try finding 'genotype:'
    regex_geno = (r'genotype:\s+(.*)')
    geno_match = re.findall(regex_geno, genbank_note)
    if len(geno_match)!=0:
        genotype = get_geno_from_match(geno_match[0])
        subgenotype = get_subgeno_from_match(genotype, geno_match[0])

    #could also be listed as 'genotype' 
    elif len(geno_match)==0:
        regex_geno2 = (r'genotype\s+(.*)')
        geno2_match = re.findall(regex_geno2, genbank_note)
        if len(geno2_match)!=0:
            genotype = get_geno_from_match(geno2_match[0])
            subgenotype = get_subgeno_from_match(genotype, geno2_match[0])
        
        #could be listed as 'genotype:' with no space after
        elif len(geno2_match)==0:
            regex_geno3 = (r'genotype:+(.*)')
            geno3_match = re.findall(regex_geno3, genbank_note)
            if len(geno3_match)!=0:
                genotype = get_geno_from_match(geno3_match[0])
                subgenotype = get_subgeno_from_match(genotype, geno3_match[0])
            
            #could be listed as 'subgenotype'
            elif len(geno3_match)==0:
                regex_subgeno = (r'subgenotype\s+(.*)')
                subgeno_match = re.findall(regex_subgeno, genbank_note)
                if len(subgeno_match)!=0:
                    genotype = get_geno_from_match(subgeno_match[0])
                    subgenotype = get_subgeno_from_match(genotype, subgeno_match[0])
                    
                elif len(subgeno_match)==0:
                    genotype = 'None'
                    subgenotype = 'None'
                    
    return (genotype, subgenotype)

def get_geno_from_match(match_result):
    if 'recombinant' in match_result or 'Recombinant' in match_result or '/' in match_result:
        #slash appears in a couple listings that are not recombinants
        if 'HBV/' in match_result:
            genotype = match_result.split('HBV/')[1].upper()
        else:
            genotype = 'recombinant'
    #deal with different ways that genotypes are specified in the genbank files
    elif 'HBV-' in match_result:
        genotype = match_result.split('HBV-')[1][0].upper()
    elif 'Hiroshima-' in match_result:
        genotype = match_result.split('Hiroshima-')[1][0].upper()
    elif 'not genotype' in match_result:
        genotype= 'None'
    elif 'subgenotype:' in match_result:
        genotype = match_result.split('subgenotype: ')[1][0].upper()
    elif 'subgenotype' in match_result:
        genotype = match_result.split('subgenotype ')[1][0].upper()
    else:
        genotype = match_result[0].upper()
    
    if genotype=="D4":
        genotype = "D"

    if genotype not in ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'recombinant']:
            genotype = 'None'

    return genotype

def get_subgeno_from_match(genotype, match_result):
    #going to miss a couple subgenotypes, but this catches most of them
    
    #going to ignore recombinants for the genotype-specific builds, so don't assign a subtype
    if len(match_result)==1 or genotype=='recombinant':
        subgenotype='None'
    else:
        #if just the subgenotype is listed, it will be max len(3) 
        if len(match_result)<=3:
            subgenotype=match_result.upper()
        elif 'ubgenotype:' in match_result:
            subgenotype = match_result.split('ubgenotype: ')[1]
            if ';' in subgenotype:
                subgenotype = subgenotype.split(';')[0].upper()
        elif 'ubgenotype' in match_result:
            subgenotype = match_result.split('ubgenotype ')[1]
            for x in ['.',')',';']:
                if x in subgenotype:
                    subgenotype = subgenotype.split(x)[0].upper()
        elif ';' in match_result:
            potential_sub = match_result.split(';')[0]
            if 1< len(potential_sub) <4:
                subgenotype = potential_sub.upper()
                if ')' in subgenotype:
                    subgenotype = subgenotype.split(')')[0].upper()
            else:
                subgenotype='None'
        elif ':' in match_result:
            potential_sub = match_result.split(':')[0]
            if 1< len(potential_sub) <4:
                subgenotype = potential_sub.upper()
            else:
                subgenotype = 'None'
        elif 'HBV/' in match_result:
            subgenotype = match_result.split('HBV/')[1]
        else:
            subgenotype = 'None'

    return subgenotype

def exclude_isolate(record):
    #exclude patent and synthetic sequences, that are not clinical isolates
    return record.annotations['data_file_division'] in ['PAT', 'SYN']


def summarise(metadata):
    counts = defaultdict(int)
    n = 0
    for data in metadata.values():
        n+=1
        for k,v in data.items():
            if v != 'None':
                counts[k]+=1
    print("\nSummary of parsed metadata:")
    for field, v in counts.items():
        stars = "*" * int(round(v/n*20, 0))
        print(f"{field:25}{stars:20} (n={v}/{n})")
    print()

def parse_collection_date(record, source):

    collection_date = extract(source, 'collection_date')
    if not collection_date:
        collection_date = record.annotations['date']
    
    if collection_date == '2018/2019':
        return False
    
    formatted_date = dateutil.parser.parse(collection_date).strftime('%Y-%m-%d')
    #dateutil parser will assign a day (today's date) to unknown days, and same for month, want XX instead
    if len(collection_date)==8:
        formatted_date = formatted_date[:-2] + 'XX'
    elif len(collection_date)==4:
        # NOTE: around 1/4 sequences only have a year!
        formatted_date = formatted_date[:5] + 'XX-XX'

    return formatted_date

def extract(feature, key):
    if key in feature.qualifiers:
         return feature.qualifiers[key][0]
    return False

def parse_metadata(record):
    accession = record.id.split('.')[0]
    source = [feature for feature in record.features if feature.type == 'source'][0]
    note = extract(source, 'note')

    metadata = {}
    metadata['name'] = accession # TODO XXX - change once new augur release is out with https://github.com/nextstrain/augur/pull/1240, 
    metadata['accession'] = accession
    metadata['strain_name'] = extract(source, 'strain') or extract(source, 'isolate') or 'None'
    metadata['country'] = extract(source, 'country') or 'None'
    metadata['host'] = extract(source, 'host') or 'None'
    metadata['genotype_genbank'], metadata['subgenotype_genbank'] = parse_genotype(note) if note else ("None", "None")
    metadata['date'] = parse_collection_date(record, source)
    if not metadata['date']:
        raise MissingDate
    return metadata

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--genbank", help="Input genbank file")
    parser.add_argument("--output-sequences")
    parser.add_argument("--output-metadata")
    args = parser.parse_args()

    n_records = 0
    seq_records = []
    exclude_counts = defaultdict(int)
    metadata = {}

    for record in SeqIO.parse(open(args.genbank,"r"), "genbank"):
        n_records+=1

        if exclude_isolate(record):
            exclude_counts['patient/synthetic sequence']+=1
            continue

        try:
            record_metadata = parse_metadata(record)
        except MissingDate:
            exclude_counts['missing date']+=1
            continue   

        seq_records.append(SeqRecord(record.seq, id=record_metadata['accession'], description=''))  
        metadata[record_metadata['accession']] = record_metadata

    print(f"{n_records} records parsed. Excluded "+', '.join([f"n={v} ({k})"for k,v in exclude_counts.items()]))
    summarise(metadata)

    SeqIO.write(seq_records, args.output_sequences, "fasta")
    with open(args.output_metadata, 'w') as fh:
        writer = csv.writer(fh, delimiter='\t', quotechar='"')
        header = list(next(iter(metadata.values())).keys())
        writer.writerow(header)
        for d in metadata.values():
            writer.writerow([d[k] for k in header])
