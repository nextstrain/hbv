#!/usr/bin/env python3
"""
Modifies country entries in the NDJSON records from stdin to split on the ':' character and discard any content after.
The modified records are output to stdout.
For instance, USA:Alaska -> USA
Also converts "None" to the empty string
"""

from sys import stdin, stdout
import argparse
import json

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    args = parser.parse_args()

for record in stdin:
        record = json.loads(record)
        record['country'] = record['country'].split(':')[0]
        if record['country'] == "None":
            record['country'] = ''
        json.dump(record, stdout, allow_nan=False, indent=None, separators=',:')
        print()