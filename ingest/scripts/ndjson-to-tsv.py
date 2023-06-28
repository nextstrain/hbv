#!/usr/bin/env python3
"""
Parses NDJSON records from stdin to two different files: a metadata TSV and a
sequences FASTA.

Records that do not have an ID or sequence will be excluded from the output files.
"""

"""
Originally taken from https://github.com/nextstrain/monkeypox/blob/62fca491c6775573ad036eedf34b271b4952f2c2/ingest/bin/ndjson-to-tsv-and-fasta
and modified as the expected use case here doesn't (yet) involve sequences in the NDJSON
"""

import argparse
import csv
import json
from sys import stderr, stdin


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--metadata", metavar="TSV", default="data/metadata.tsv",
        help="The output metadata TSV file")
    parser.add_argument("--metadata-columns", nargs="+",
        help="List of fields from the NDJSON records to include as columns in the metadata TSV. " +
             "Metadata TSV columns will be in the order of the columns provided.")

    args = parser.parse_args()

    with open(args.metadata, 'wt') as metadata_output:
        metadata_csv = csv.DictWriter(
            metadata_output,
            args.metadata_columns,
            restval="",
            extrasaction='ignore',
            delimiter='\t'
        )
        metadata_csv.writeheader()

        for index, record in enumerate(stdin):
            record = json.loads(record)
            metadata_csv.writerow(record)