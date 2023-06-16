#!/usr/bin/env python3

"""
Most code was originally written by @kistlerk as part of `curate_genbank_to_fasta.ipynb`
[blab/adaptive-evolution](https://github.com/blab/adaptive-evolution)
"""

import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from dateutil.parser import parse
from collections import defaultdict, Counter
import dateutil
import re
import csv

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
                    
    return genotype, subgenotype

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
        if genotype not in ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I']:
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
                    subgenotyppe = subgenotype.split(')')[0].upper()
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

def main(input_fname, output_seqs_fname, output_meta_fname):
    """
    Takes a genbank file and
    1. Removes synthetic sequences and clinical isolates
    2. Extracts (interpets) metadata: collection_date, strain, country, host, genotype, subgenotype
    3. Drops sequences without a collection_date
    4. Writes out FASTA with header accession|strain_name|date|country|host|genotype|subgenotype
    """

    #store all edited sequence records
    seq_records = []
    counts = defaultdict(int)
    metadata = {}
    for record in SeqIO.parse(open(input_fname,"r"), "genbank"):
        counts['total']+=1
        accession = record.id.split('.')[0]
        collection_date, strain_name, country, host, genotype, subgenotype = 'None', 'None', 'None', 'None', 'None', 'None'
        #exclude patent and synthetic sequences, that are not clinical isolates
        if record.annotations['data_file_division'] not in ['PAT', 'SYN']:
            
            for feature in record.features:
                if feature.type == 'source':
                    if 'collection_date' in feature.qualifiers:
                        collection_date = feature.qualifiers['collection_date'][0]   
                    if 'strain' in feature.qualifiers:
                        strain_name = feature.qualifiers['strain'][0]
                    if 'country' in feature.qualifiers:
                        country = feature.qualifiers['country'][0]
                    if 'host' in feature.qualifiers:
                        host = feature.qualifiers['host'][0]
                    if 'note' in feature.qualifiers:
                        if 'genotype' in feature.qualifiers['note'][0]:
                            genotype, subgenotype = parse_genotype(feature.qualifiers['note'][0])
            
            if genotype == 'None':
                counts['no_genotype']+=1
            if subgenotype == 'None':
                counts['no_subgenotype']+=1

            if collection_date == 'None':
                collection_date = record.annotations['date']

            #only keep sequences with date
            
            if collection_date != 'None':
                if collection_date != '2018/2019':
                    formatted_date = dateutil.parser.parse(collection_date).strftime('%Y-%m-%d')
                    #dateutil parser will assign a day (today's date) to unknown days, and same for month, want XX instead
                    if len(collection_date)==8:
                        formatted_date = formatted_date[:-2] + 'XX'
                    elif len(collection_date)==4:
                        formatted_date = formatted_date[:5] + 'XX-XX'
                    
                    
                    counts['output']+=1
                    metadata[accession] = {
                        'name': accession, # TODO XXX - change once new augur release is out with https://github.com/nextstrain/augur/pull/1240,
                        'accession': accession,
                        'strain_name': strain_name,
                        'collection_date': formatted_date,
                        'country': country,
                        'host': host,
                        'genotype_genbank': genotype,
                        'subgenotype_genbank': subgenotype
                    }
                    # list_of_info = [accession, strain_name, formatted_date, country, host, genotype, subgenotype]
                    # new_record_info = '|'.join(list_of_info)
                    seq_records.append(SeqRecord(record.seq, id=accession, description=''))  
                else:
                    counts['excluded_missing_date']+=1
            else:
                counts['excluded_missing_date']+=1

    print(counts)

    #write fasta sequence file for ALL sequences
    SeqIO.write(seq_records, output_seqs_fname, "fasta")

    with open(output_meta_fname, 'w') as fh:
        writer = csv.writer(fh, delimiter='\t', quotechar='"')
        header = list(next(iter(metadata.values())).keys())
        writer.writerow(header)
        for d in metadata.values():
            writer.writerow([d[k] for k in header])

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--genbank", help="Input genbank file")
    parser.add_argument("--output-sequences")
    parser.add_argument("--output-metadata")
    args = parser.parse_args()
    main(args.genbank, args.output_sequences, args.output_metadata)
