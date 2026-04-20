#!/usr/bin/env python3
import argparse
import pandas as pd

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--metadata", required=True)
    ap.add_argument("--metadata-genbank", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--accession-col", default="accession")
    args = ap.parse_args()

    meta = pd.read_csv(args.metadata, sep="\t", dtype=str)
    gb   = pd.read_csv(args.metadata_genbank, sep="\t", dtype=str)

    key = args.accession_col

    meta = meta[meta[key].notna() & (meta[key] != "")]
    gb   = gb[gb[key].notna() & (gb[key] != "")]

    extra_cols = [c for c in gb.columns if c != key and c not in meta.columns]

    out = meta.merge(gb[[key] + extra_cols], on=key, how="left")
    out["url"] = "https://www.ncbi.nlm.nih.gov/nuccore/" + out[key].astype(str)

    out.to_csv(args.out, sep="\t", index=False)

if __name__ == "__main__":
    main()