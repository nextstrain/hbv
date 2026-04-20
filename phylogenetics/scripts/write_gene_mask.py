#!/usr/bin/env python3
import argparse
from pathlib import Path

from Bio import SeqIO


def read_regions(path):
    regions = {}
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue

            name, start, end = line.split()[:3]
            regions[name] = (int(start), int(end))

    return regions


def alignment_length(path, reference_id):
    first_record = None

    for record in SeqIO.parse(path, "fasta"):
        if first_record is None:
            first_record = record

        if record.id == reference_id or record.id.split()[0] == reference_id:
            return len(record.seq)

    if first_record is not None:
        return len(first_record.seq)

    raise ValueError(f"No sequences found in {path}")


def genbank_length(path):
    return len(SeqIO.read(path, "genbank").seq)


def mask_positions(start, end, length):
    if start < 1 or end < 1 or start > length or end > length:
        raise ValueError(
            f"Invalid coordinates: start={start}, end={end}, alignment length={length}"
        )

    if start <= end:
        yield from range(1, start)
        yield from range(end + 1, length + 1)
    else:
        yield from range(end + 1, start)


def main():
    parser = argparse.ArgumentParser(
        description="Write augur mask positions for all alignment sites outside a gene."
    )
    parser.add_argument("--regions", required=True)
    parser.add_argument("--alignment")
    parser.add_argument("--reference-genbank")
    parser.add_argument("--gene", required=True)
    parser.add_argument("--reference")
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    if args.alignment and not args.reference:
        parser.error("--reference is required when using --alignment")
    if not args.alignment and not args.reference_genbank:
        parser.error("one of --alignment or --reference-genbank is required")

    regions = read_regions(args.regions)
    if args.gene not in regions:
        raise ValueError(f"Gene {args.gene!r} not found in {args.regions}")

    if args.alignment:
        length = alignment_length(args.alignment, args.reference)
    else:
        length = genbank_length(args.reference_genbank)

    start, end = regions[args.gene]

    if args.output == "-":
        for position in mask_positions(start, end, length):
            print(position)
        return

    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)

    with open(output, "w") as fh:
        for position in mask_positions(start, end, length):
            fh.write(f"{position}\n")


if __name__ == "__main__":
    main()
