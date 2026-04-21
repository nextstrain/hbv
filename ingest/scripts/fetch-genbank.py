#!/usr/bin/env python3

import argparse
import json
import random
import time
from http.client import IncompleteRead
from urllib.error import HTTPError, URLError

from Bio import Entrez


Entrez.email = "hello@nextstrain.org"
BATCH_SIZE = 1000


def parse_args():
    parser = argparse.ArgumentParser(
        description="Fetch HBV GenBank records from Entrez.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--term", required=True, help="Entrez search term")
    parser.add_argument("--output", required=True, help="Output GenBank file")
    parser.add_argument(
        "--accessions",
        help=(
            "Optional newline-delimited accession list. When provided, fetch only "
            "these records instead of the full Entrez query result."
        ),
    )
    return parser.parse_args()


def get_esearch_history(term):
    handle = Entrez.esearch(
        db="nucleotide",
        term=term,
        retmode="json",
        usehistory="y",
        retmax=0,
    )
    esearch_result = json.loads(handle.read())["esearchresult"]
    print(f"Search term {term!r} returned {esearch_result['count']} IDs.")
    return {
        "count": int(esearch_result["count"]),
        "query_key": esearch_result["querykey"],
        "web_env": esearch_result["webenv"],
    }


def fetch_query_batch(query_key, web_env, start, retmax, tries=8):
    for attempt in range(tries):
        try:
            handle = Entrez.efetch(
                db="nucleotide",
                query_key=query_key,
                webenv=web_env,
                retstart=start,
                retmax=retmax,
                rettype="gb",
                retmode="text",
            )
            return handle.read()
        except (IncompleteRead, HTTPError, URLError, OSError):
            time.sleep(min(60, (2**attempt) + random.random()))
    raise RuntimeError(f"efetch failed after {tries} tries at retstart={start}")


def fetch_accession_batch(accessions, tries=8):
    for attempt in range(tries):
        try:
            handle = Entrez.efetch(
                db="nucleotide",
                id=",".join(accessions),
                rettype="gb",
                retmode="text",
            )
            return handle.read()
        except (IncompleteRead, HTTPError, URLError, OSError):
            time.sleep(min(60, (2**attempt) + random.random()))
    raise RuntimeError(
        "efetch failed after "
        f"{tries} tries for accession batch starting with {accessions[0]!r}"
    )


def fetch_from_accessions(accession_file, out_path):
    with open(accession_file, "r") as fh:
        accessions = [line.strip() for line in fh if line.strip()]

    if not accessions:
        raise RuntimeError(
            "Development-mode GenBank fetch needs accessions from the active "
            "NCBI NDJSON, but the accession list was empty."
        )

    print(
        "Development mode enabled; fetching GenBank records for "
        f"{len(accessions)} active NCBI accessions."
    )
    with open(out_path, "w") as output_handle:
        written = 0
        for start in range(0, len(accessions), BATCH_SIZE):
            batch_accessions = accessions[start : start + BATCH_SIZE]
            records = fetch_accession_batch(batch_accessions)
            output_handle.write(records)
            output_handle.flush()
            written += records.count("\nLOCUS")
            print(f"[batch] total_written={written}")
            time.sleep(0.4)


def fetch_from_query(term, out_path):
    history = get_esearch_history(term)
    count = history["count"]
    query_key = history["query_key"]
    web_env = history["web_env"]

    print(f"Fetching GenBank records in batches of n={BATCH_SIZE}")
    with open(out_path, "w") as output_handle:
        written = 0
        for start in range(0, count, BATCH_SIZE):
            batch_size = min(BATCH_SIZE, count - start)
            records = fetch_query_batch(query_key, web_env, start, batch_size)
            output_handle.write(records)
            output_handle.flush()
            written += records.count("\nLOCUS")
            print(f"[batch] total_written={written}")
            time.sleep(0.4)


def main():
    args = parse_args()
    if args.accessions:
        fetch_from_accessions(args.accessions, args.output)
    else:
        fetch_from_query(args.term, args.output)


if __name__ == "__main__":
    main()
