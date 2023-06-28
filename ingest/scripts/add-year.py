#!/usr/bin/env python3

from sys import stdin, stdout
import json

for record in stdin:
    record = json.loads(record)
    record['year'] = record['date'][0:4]
    json.dump(record, stdout, allow_nan=False, indent=None, separators=',:')
    print()