#!/usr/bin/env python3


import argparse
import csv

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--meta-in")
    parser.add_argument("--meta-out")
    parser.add_argument("--nextalign")
    args = parser.parse_args()

    nextalign = {}
    with open(args.nextalign) as csvfile:
        reader = csv.DictReader(csvfile, delimiter=",")
        for row in reader:
            nextalign[row['seqName']] = row

    with open(args.meta_in, 'r') as csvfile:
        reader = csv.DictReader(csvfile, delimiter="\t")
        header = reader.fieldnames

        with open(args.meta_out, 'w') as fh:
            writer = csv.writer(fh, delimiter='\t', quotechar='"')
            writer.writerow([*header, "nextalign_error"])

            for metadata_row in reader:
                s = ''
                err_str = nextalign.get(metadata_row['accession'], {}).get('errors', '')
                if err_str:
                    if "low seed matching rate" in err_str:
                        s = 'low seed matching rate'
                    elif 'seed matches suggest large indels or are ambiguous due to duplications' in err_str:
                        s = 'indels/duplications'
                    elif 'not enough matches' in err_str:
                        s = 'not enough matches'
                    else:
                        print("UKNONWN ERROR:", err_str)
                writer.writerow([*[metadata_row[field] for field in header], s])
