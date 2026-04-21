#!/usr/bin/env python3
import argparse
import sys

import pandas as pd


def format_examples(values):
    values = sorted(values)
    return ", ".join(values[:5]) if values else "none"


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

    meta_ids = set(meta[key].astype(str))
    gb_ids = set(gb[key].astype(str))
    overlap = meta_ids & gb_ids
    overlap_fraction = (len(overlap) / len(meta_ids)) if meta_ids else 0.0

    if not overlap:
        raise SystemExit(
            "GenBank metadata merge failed: no accession overlap between curated "
            "NCBI metadata and curated GenBank metadata. "
            f"NCBI accessions={len(meta_ids)}, GenBank accessions={len(gb_ids)}. "
            "Likely reasons: the GenBank fetch did not target the active NCBI "
            "accessions, or many targeted GenBank records were dropped during "
            "parsing/curation (for example due to missing dates or excluded "
            "PAT/SYN records). "
            f"Example NCBI-only accessions: {format_examples(meta_ids - gb_ids)}. "
            f"Example GenBank-only accessions: {format_examples(gb_ids - meta_ids)}."
        )

    if overlap_fraction < 0.5:
        print(
            "Warning: poor accession overlap between curated NCBI metadata and "
            "curated GenBank metadata. "
            f"Overlap={len(overlap)}/{len(meta_ids)} ({overlap_fraction:.1%} of "
            "NCBI accessions). "
            f"Example NCBI-only accessions: {format_examples(meta_ids - gb_ids)}. "
            f"Example GenBank-only accessions: {format_examples(gb_ids - meta_ids)}.",
            file=sys.stderr,
        )

    extra_cols = [c for c in gb.columns if c != key and c not in meta.columns]

    out = meta.merge(gb[[key] + extra_cols], on=key, how="left")
    out["url"] = "https://www.ncbi.nlm.nih.gov/nuccore/" + out[key].astype(str)

    out.to_csv(args.out, sep="\t", index=False)

if __name__ == "__main__":
    main()
