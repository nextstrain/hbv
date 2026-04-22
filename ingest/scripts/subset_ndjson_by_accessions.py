#!/usr/bin/env python3

import argparse
import json
import os


VERBOSE = os.environ.get("HBV_VERBOSE", "").lower() in {"1", "true", "yes", "on"}


def parse_args():
    parser = argparse.ArgumentParser(
        description="Filter NDJSON records to a supplied accession list.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--input", required=True, help="Input NDJSON file")
    parser.add_argument(
        "--accessions",
        required=True,
        help="Newline-delimited accession list used to select NDJSON records",
    )
    parser.add_argument("--output", required=True, help="Output NDJSON file")
    return parser.parse_args()


def load_accessions(path):
    with open(path, "r") as fh:
        accessions = {line.strip() for line in fh if line.strip()}
    if not accessions:
        raise RuntimeError("No accessions were available for active NCBI selection.")
    return accessions


def main():
    args = parse_args()
    keep = load_accessions(args.accessions)
    matched = set()
    written = 0

    with open(args.input, "r") as input_handle, open(args.output, "w") as output_handle:
        for line in input_handle:
            record = json.loads(line)
            accession = record.get("accession")
            if accession in keep:
                output_handle.write(line)
                matched.add(accession)
                written += 1

    missing = len(keep - matched)
    if VERBOSE:
        print(
            "Selected "
            f"{written} NCBI records from {len(keep)} Entrez-derived accessions "
            f"(missing_in_ncbi={missing})."
        )


if __name__ == "__main__":
    main()
