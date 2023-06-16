#!/usr/bin/env python3

"""
HBV is a circular genome or ~3.2kb, and so a origin point has to be chosen to represent it in fasta/genbank
Upon alignment of all (NCBI) genomes to reference (JN182318) ~10% of (full-length) genomes had very good alignments
which terminated part-way through the reference as they were using a different (3') origin.

This script uses a simple seed-matching approach to find where the 3' end of the reference best matches,
and if it's suitably far into the genome we shift the genome accordingly.

Adds the 'circularise' field to the metadata TSV
                                                                                    @jameshadfield June 2023
"""

import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import csv

def analyse_ref(seq_fname, ref_name):
    ref = None
    records = {}
    for record in SeqIO.parse(open(seq_fname,"r"), "fasta"):
        records[record.name] = record
        if record.name == ref_name:
            ref = record
    if not ref:
        raise Exception("Ref not found")
    # print(ref.seq[0:20])
    seed_len = 30
    seeds = [
        {'start':0, 'seq': str(ref.seq[0:seed_len])},
        {'start':100, 'seq': str(ref.seq[100:100+seed_len])},
        {'start':200, 'seq': str(ref.seq[200:200+seed_len])},
        {'start':300, 'seq': str(ref.seq[300:300+seed_len])}
    ]
    # print(seeds)
    return (ref, records, seeds)

def seq_diff(a, b):
    return sum([x[0]!=x[1] for x in zip([*a], [*b])])

def print_match(a, b):
    if len(a)>100:
        a=a[0:100]
        b=b[0:100]
    print(f"\t{a}")
    print(f"\t{''.join(['|' if aa==b[i] else ' ' for i,aa in enumerate(a)])}")
    print(f"\t{b}")
# print_match("aaa", "aba")


def identify_origin(records, seeds, verbose=0):
    BAD_SEED_MISMATCH_COUNT = 10 # TODO XXX make argument
    START_BUFFER = 300 # TODO XXX make argument
    
    count = 0
    origins = {}
    for name, record in records.items():
        count+=1

        if count%1000==0:
            print(count)
        mismatches = 10000
        seed_used = None
        match_start = None
        try:
            for seed in seeds:
                l = len(seed['seq'])
                for i in range(0, len(record.seq)-l-seed['start']-1):
                    c = seq_diff(seed['seq'], str(record.seq[i:i+l]))
                    if c<mismatches:
                        mismatches = c
                        match_start = i
                        seed_used = seed
                    if c==0:
                        raise StopIteration
        except StopIteration:
            pass
                
        origin = match_start - seed_used['start']
        origins[record.name] = {
            'origin': origin,
            'start_pos_seed': match_start,
            'seed_used': seed,
            'mismatches': mismatches,
            'recut': origin > START_BUFFER and mismatches<BAD_SEED_MISMATCH_COUNT,
            'bad_match': mismatches>=BAD_SEED_MISMATCH_COUNT,
        }
        
        if origin < 0 and verbose>0:
            print(f"\t{record.name} origin is at 5' end!")
            
        if verbose>0:
            print(f"{record.name} start={origin} {mismatches} mismatches (seed@{seed_used['start']})" +
                  (" RECUT" if origins[name]['recut'] else "") +
                  (" BAD MATCH" if origins[name]['bad_match'] else ""))
            if verbose>1 and origins[name]['recut']:
                print_match(seed_used['seq'], record.seq[match_start:match_start+len(seed_used['seq'])])
                print("\n")
                # print_match(reference.seq[0:30], record.seq[origin:origin+30])
                # print("\n")
            
    print(f"Matches found:    {len([v for v in origins.values() if not v['bad_match']])}/{len(origins)}")
    print(f"Genomes to recut: {len([v for v in origins.values() if v['recut']])}/{len(origins)}")
    print(f"Bad Matches:      {len([v for v in origins.values() if v['bad_match']])}/{len(origins)}")
    return origins

def recircularise(records, origins, reference, verbose=False):
    records_to_write = []
    for name, record in records.items():
        if origins[name]['recut']:
            seed_offset = origins[name]['seed_used']['start']
            start = origins[name]['origin']
            new_record = SeqRecord(
                Seq(record.seq[start:] + record.seq[:start]),
                name = name,
                id = name,
                description = ''
            )
            assert(len(new_record.seq)==len(record.seq))
            records_to_write.append(new_record)
    
            if verbose:
                print(f"{name} - 3 prime pseudo-alignment (note: seed matched at {seed_offset})")
                print_match(reference.seq, new_record.seq)

                print(f"{name} - starting at {seed_offset}bp")
                print_match(reference.seq[seed_offset:], new_record.seq[seed_offset:])        

                print("\n")

        # elif origins[name]['bad_match']:
        #     if verbose:
        #         print(f"DROPPING {name}")

        else:
            records_to_write.append(record)
    return records_to_write

def append_to_metadata(origins, fname_in, fname_out):

    ## TODO -- shift to augur commands / stream processing
    with open(fname_in, 'r') as csvfile:
        reader = csv.DictReader(csvfile, delimiter="\t")
        header = reader.fieldnames

        with open(fname_out, 'w') as fh:
            writer = csv.writer(fh, delimiter='\t', quotechar='"')
            writer.writerow([*header, "circularise", "circularise_shift_bp"])

            for row in reader:
                o = origins.get(row['accession'], {})
                s = ['', '']
                if o.get('bad_match', False):
                    s[0] = 'no seed match'
                elif o.get('recut', False):
                    s[0] = 'shifted'
                    s[1] = o['origin']
                writer.writerow([*[row[field] for field in header], *s])

def main(args):
    reference, records, seeds = analyse_ref(args.seqs_in, args.reference)
    print("seeds:", seeds)

    origins = identify_origin(records, seeds, 0)

    records_to_write = recircularise(records, origins, reference, False)

    with open(args.seqs_out, 'w') as fh:
        SeqIO.write(records_to_write, fh, "fasta")
    
    append_to_metadata(origins, args.meta_in, args.meta_out)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--seqs-in")
    parser.add_argument("--meta-in")
    parser.add_argument("--seqs-out")
    parser.add_argument("--meta-out")
    parser.add_argument("--reference")
    args = parser.parse_args()
    main(args)