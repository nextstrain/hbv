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
    parser.add_argument("--reference-genbank", required=True)
    parser.add_argument("--gene", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    regions = read_regions(args.regions)
    if args.gene not in regions:
        raise ValueError(f"Gene {args.gene!r} not found in {args.regions}")

    length = genbank_length(args.reference_genbank)
    start, end = regions[args.gene]

    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)

    with open(output, "w") as fh:
        for position in mask_positions(start, end, length):
            fh.write(f"{position}\n")

if __name__ == "__main__":
    main()
