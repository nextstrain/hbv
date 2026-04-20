#!/usr/bin/env python3
"""
Modifies country entries in the NDJSON records from stdin to split on the ':' character and discard any content after.
The modified records are output to stdout.
For instance, USA:Alaska -> USA

Priority:
1) country
2) location
3) geo-location

Splits on ':' and drops trailing content.
Converts missing or "None" to empty string.
"""

from sys import stdin, stdout
import argparse
import json

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.parse_args()

    for line in stdin:
        record = json.loads(line)

        c = record.get("country") or record.get("location") or record.get("geo-location")
        if c:
            c = c.split(":")[0]
            if c == "None":
                c = ""
        else:
            c = ""

        record["country"] = c

        json.dump(record, stdout, allow_nan=False, indent=None, separators=(",", ":"))
        print()
