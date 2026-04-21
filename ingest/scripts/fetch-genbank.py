#!/usr/bin/env python3
"""Fetch HBV GenBank records from Entrez for either full or filtered ingest runs.

Support the main full-query mode plus optional development and complete-genome restrictions.
"""

import argparse
import json
import os
import random
import time
from http.client import IncompleteRead
from urllib.error import HTTPError, URLError

from Bio import Entrez


Entrez.email = "hello@nextstrain.org"
BATCH_SIZE = 1000
VERBOSE = os.environ.get("HBV_VERBOSE", "").lower() in {"1", "true", "yes", "on"}


def vprint(*args, **kwargs):
    if VERBOSE:
        print(*args, **kwargs)


def ensure_text(records):
    """Normalize Entrez reads to text before writing GenBank output."""
    if isinstance(records, bytes):
        return records.decode("utf-8", errors="replace")
    return records


def count_records(records):
    """Count GenBank records for progress logging."""
    return records.count("\nLOCUS") + records.startswith("LOCUS")


def add_complete_genome_clause(term):
    return f"({term}) AND (complete genome[All Fields])"


def parse_args():
    parser = argparse.ArgumentParser(
        description="Fetch HBV GenBank records from Entrez.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--term", required=True, help="Entrez search term")
    parser.add_argument("--output", required=True, help="Output GenBank file")
    parser.add_argument(
        "--limit",
        type=int,
        help=(
            "Optional limit on the number of Entrez records to fetch from the "
            "start of the query result."
        ),
    )
    parser.add_argument(
        "--complete-genomes",
        action="store_true",
        help="Restrict the Entrez query to records annotated as complete genomes.",
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
    vprint(f"Search term {term!r} returned {esearch_result['count']} IDs.")
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
            return ensure_text(handle.read())
        except (IncompleteRead, HTTPError, URLError, OSError):
            time.sleep(min(60, (2**attempt) + random.random()))
    raise RuntimeError(f"efetch failed after {tries} tries at retstart={start}")


def fetch_from_query(term, out_path, complete_genomes=False, limit=None):
    if complete_genomes:
        term = add_complete_genome_clause(term)

    history = get_esearch_history(term)
    count = history["count"]
    if limit is not None:
        if limit <= 0:
            raise RuntimeError("--limit must be a positive integer.")
        count = min(count, limit)
        vprint(f"Development mode enabled; fetching the first {count} Entrez records.")

    if complete_genomes:
        vprint("Restricting Entrez fetch to complete genomes.")

    query_key = history["query_key"]
    web_env = history["web_env"]

    vprint(f"Fetching GenBank records in batches of n={BATCH_SIZE}")
    with open(out_path, "w") as output_handle:
        written = 0
        for start in range(0, count, BATCH_SIZE):
            batch_size = min(BATCH_SIZE, count - start)
            records = fetch_query_batch(query_key, web_env, start, batch_size)
            output_handle.write(records)
            output_handle.flush()
            written += count_records(records)
            vprint(f"[batch] total_written={written}")
            time.sleep(0.4)


def main():
    args = parse_args()
    fetch_from_query(
        args.term,
        args.output,
        complete_genomes=args.complete_genomes,
        limit=args.limit,
    )


if __name__ == "__main__":
    main()
