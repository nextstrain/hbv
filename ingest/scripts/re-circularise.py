#!/usr/bin/env python3

"""
HBV is a circular genome or ~3.2kb, and so a origin point has to be chosen to represent it in fasta/genbank
Some genomes terminate part-way through the reference because they use a different
origin.

This script uses a simple seed-matching approach to find where the 3' end of the reference best matches,
and if it's suitably far into the genome we shift the genome accordingly.

Adds the 'circularise' field to the metadata TSV
                                                                                    @jameshadfield June 2023
"""

import argparse
import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import csv

VERBOSE = os.environ.get("HBV_VERBOSE", "").lower() in {"1", "true", "yes", "on"}


def vprint(*args, **kwargs):
    if VERBOSE:
        print(*args, **kwargs)


def load_reference_from_genbank(ref_genbank, ref_name):
    for record in SeqIO.parse(open(ref_genbank, "r"), "genbank"):
        accession = record.id.split(".")[0]
        if accession == ref_name:
            return SeqRecord(
                record.seq,
                id=accession,
                name=accession,
                description="",
            )
    return None


def analyse_ref(seq_fname, ref_name, ref_genbank=None):
    ref = None
    records = {}
    for record in SeqIO.parse(open(seq_fname, "r"), "fasta"):
        records[record.name] = record
        if record.name == ref_name:
            ref = record
    if not ref:
        if ref_genbank:
            ref = load_reference_from_genbank(ref_genbank, ref_name)
        if not ref:
            source = ref_genbank or seq_fname
            raise Exception(f"Reference {ref_name!r} not found in {source}")
        vprint(f"Loaded reference {ref_name} from {ref_genbank}")
    seed_len = 30
    seeds = [
        {"start": 0, "seq": str(ref.seq[0:seed_len])},
        {"start": 100, "seq": str(ref.seq[100 : 100 + seed_len])},
        {"start": 200, "seq": str(ref.seq[200 : 200 + seed_len])},
        {"start": 300, "seq": str(ref.seq[300 : 300 + seed_len])},
    ]
    return (ref, records, seeds)


def seq_diff(a, b):
    return sum([x[0] != x[1] for x in zip([*a], [*b])])


def print_match(a, b):
    if len(a) > 100:
        a = a[0:100]
        b = b[0:100]
    vprint(f"\t{a}")
    vprint(f"\t{''.join(['|' if aa == b[i] else ' ' for i, aa in enumerate(a)])}")
    vprint(f"\t{b}")


def identify_origin(records, seeds, verbose=0):
    BAD_SEED_MISMATCH_COUNT = 10
    START_BUFFER = 300

    count = 0
    origins = {}
    skipped = 0

    vprint("Processing", len(records), "items...")

    for name, record in records.items():
        count += 1
        if count % 1000 == 0:
            vprint(count)

        mismatches = 10000
        seed_used = None
        match_start = None

        try:
            for seed in seeds:
                seed_len = len(seed["seq"])
                max_i = len(record.seq) - seed_len - seed["start"] - 1
                if max_i <= 0:
                    continue

                for i in range(0, max_i):
                    c = seq_diff(seed["seq"], str(record.seq[i : i + seed_len]))
                    if c < mismatches:
                        mismatches = c
                        match_start = i
                        seed_used = seed
                    if c == 0:
                        raise StopIteration
        except StopIteration:
            pass

        if seed_used is None or match_start is None:
            skipped += 1
            origins[record.name] = {
                "origin": None,
                "start_pos_seed": None,
                "seed_used": None,
                "mismatches": None,
                "recut": False,
                "bad_match": True,
            }
            continue

        origin = match_start - seed_used["start"]
        origins[record.name] = {
            "origin": origin,
            "start_pos_seed": match_start,
            "seed_used": seed_used,
            "mismatches": mismatches,
            "recut": origin > START_BUFFER and mismatches < BAD_SEED_MISMATCH_COUNT,
            "bad_match": mismatches >= BAD_SEED_MISMATCH_COUNT,
        }

        if origin < 0 and verbose > 0:
            print(f"\t{record.name} origin is at 5' end!")

        if verbose > 0:
            print(
                f"{record.name} start={origin} {mismatches} mismatches (seed@{seed_used['start']})"
                + (" RECUT" if origins[name]["recut"] else "")
                + (" BAD MATCH" if origins[name]["bad_match"] else "")
            )
            if verbose > 1 and origins[name]["recut"]:
                print_match(
                    seed_used["seq"],
                    record.seq[match_start : match_start + len(seed_used["seq"])],
                )
                print("\n")

    vprint(f"Skipped (no seed search possible): {skipped}/{len(records)}")
    vprint(
        f"Matches found:    {len([v for v in origins.values() if not v['bad_match']])}/{len(origins)}"
    )
    vprint(
        f"Genomes to recut: {len([v for v in origins.values() if v['recut']])}/{len(origins)}"
    )
    vprint(
        f"Bad Matches:      {len([v for v in origins.values() if v['bad_match']])}/{len(origins)}"
    )
    return origins


def recircularise(records, origins, reference, verbose=False):
    records_to_write = []
    for name, record in records.items():
        if origins[name]["recut"]:
            seed_offset = origins[name]["seed_used"]["start"]
            start = origins[name]["origin"]
            new_record = SeqRecord(
                Seq(record.seq[start:] + record.seq[:start]),
                name=name,
                id=name,
                description="",
            )
            assert len(new_record.seq) == len(record.seq)
            records_to_write.append(new_record)

            if verbose:
                print(
                    f"{name} - 3 prime pseudo-alignment (note: seed matched at {seed_offset})"
                )
                print_match(reference.seq, new_record.seq)

                print(f"{name} - starting at {seed_offset}bp")
                print_match(reference.seq[seed_offset:], new_record.seq[seed_offset:])

                print("\n")

        else:
            records_to_write.append(record)
    return records_to_write


def append_to_metadata(origins, fname_in, fname_out):

    # Future improvement: replace this with an augur-based or stream-oriented implementation.
    with open(fname_in, "r") as csvfile:
        reader = csv.DictReader(csvfile, delimiter="\t")
        header = reader.fieldnames

        with open(fname_out, "w") as fh:
            writer = csv.writer(fh, delimiter="\t", quotechar='"')
            writer.writerow([*header, "circularise", "circularise_shift_bp"])

            for row in reader:
                o = origins.get(row["accession"], {})
                s = ["", ""]
                if o.get("bad_match", False):
                    s[0] = "no seed match"
                elif o.get("recut", False):
                    s[0] = "shifted"
                    s[1] = o["origin"]
                writer.writerow([*[row[field] for field in header], *s])


def main(args):
    reference, records, seeds = analyse_ref(
        args.seqs_in,
        args.reference,
        args.reference_genbank,
    )
    vprint("seeds:", seeds)

    origins = identify_origin(records, seeds, 0)

    records_to_write = recircularise(records, origins, reference, False)

    with open(args.seqs_out, "w") as fh:
        SeqIO.write(records_to_write, fh, "fasta")

    append_to_metadata(origins, args.meta_in, args.meta_out)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--seqs-in")
    parser.add_argument("--meta-in")
    parser.add_argument("--seqs-out")
    parser.add_argument("--meta-out")
    parser.add_argument("--reference")
    parser.add_argument("--reference-genbank")
    args = parser.parse_args()
    main(args)
