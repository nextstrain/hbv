#!/usr/bin/env python3
"""
Copied from https://github.com/nextstrain/monkeypox/blob/62fca491c6775573ad036eedf34b271b4952f2c2/ingest/bin/csv-to-ndjson
And modified to read TSV input rather than CSV

Convert TSV on stdin to NDJSON on stdout.
"""
import csv
import json
from sys import stdin, stdout

# 200 MiB; default is 128 KiB
csv.field_size_limit(200 * 1024 * 1024)

for row in csv.DictReader(stdin, delimiter="\t"):
    json.dump(row, stdout, allow_nan = False, indent = None, separators = ',:')
    print()