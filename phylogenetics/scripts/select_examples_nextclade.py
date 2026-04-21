"""Select a small HBV example panel for the Nextclade dataset.

Support curated genotype-balanced sampling or random sampling with recombinant enrichment.
"""

import argparse
import pandas as pd

GENOTYPES = list("ABCDEFGHI")
QC_LEVELS = ["good", "mediocre", "bad"]


def is_recombinant(df):
    return (
        df["clade_nextclade"].str.contains("_re", na=False)
        | df["genotype_genbank"].str.lower().eq("recombinant")
        | (
            (df["genotype_genbank"] != "")
            & (df["clade_nextclade"] != "")
            & (df["clade_nextclade"] != "unassigned")
            & (df["genotype_genbank"] != df["clade_nextclade"])
        )
    )


def sample_random(df, n):
    if n <= 0 or df.empty:
        return df.iloc[0:0].copy()
    return df.sample(n=min(n, len(df)), random_state=42).copy()


def take_evenly_across_qc(df, n, prefer_explicit_recombinant=False):
    if n <= 0 or df.empty:
        return df.iloc[0:0].copy()

    x = df.copy()
    if prefer_explicit_recombinant:
        x["explicit_recombinant"] = x["genotype_genbank"].str.lower().eq("recombinant")

    chunks = []
    base = n // len(QC_LEVELS)
    remainder = n % len(QC_LEVELS)

    for i, qc in enumerate(QC_LEVELS):
        n_take = base + (1 if i < remainder else 0)
        if n_take == 0:
            continue

        chunk = x[x["QC_overall_status"] == qc].copy()
        if prefer_explicit_recombinant and not chunk.empty:
            explicit = sample_random(chunk[chunk["explicit_recombinant"]], n_take)
            if len(explicit) < n_take:
                rest = sample_random(
                    chunk[~chunk["accession"].isin(explicit["accession"])],
                    n_take - len(explicit),
                )
                chunk = pd.concat([explicit, rest], ignore_index=True)
            else:
                chunk = explicit
        else:
            chunk = sample_random(chunk, n_take)

        chunks.append(chunk)

    out = pd.concat(chunks, ignore_index=True) if chunks else x.iloc[0:0].copy()

    if len(out) < n:
        used = set(out["accession"]) if not out.empty else set()
        fill = x[~x["accession"].isin(used)].copy()

        if prefer_explicit_recombinant:
            explicit = sample_random(fill[fill["explicit_recombinant"]], n - len(out))
            if len(explicit) < (n - len(out)):
                rest = sample_random(
                    fill[~fill["accession"].isin(explicit["accession"])],
                    (n - len(out)) - len(explicit),
                )
                fill = pd.concat([explicit, rest], ignore_index=True)
            else:
                fill = explicit
        else:
            fill = sample_random(fill, n - len(out))

        out = pd.concat([out, fill], ignore_index=True)

    return out.head(n).copy()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--metadata", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--n-total", type=int, required=True)
    parser.add_argument("--recomb-frac", type=float, required=True)
    parser.add_argument("--sampling", choices=["curated", "random"], required=True)
    args = parser.parse_args()

    metadata = pd.read_csv(args.metadata, sep="\t", dtype=str).fillna("")

    n_recomb = round(args.n_total * args.recomb_frac)
    n_nonrec = args.n_total - n_recomb

    nonrec = metadata[~is_recombinant(metadata)].copy()
    rec = metadata[is_recombinant(metadata)].copy()

    if args.sampling == "curated":
        per_gt = n_nonrec // len(GENOTYPES)
        remainder = n_nonrec % len(GENOTYPES)

        rows = []
        for i, gt in enumerate(GENOTYPES):
            n_take = per_gt + (1 if i < remainder else 0)
            if n_take == 0:
                continue
            x = nonrec[nonrec["clade_nextclade"] == gt].copy()
            x = take_evenly_across_qc(x, n_take, prefer_explicit_recombinant=False)
            x["example_reason"] = f"genotype_{gt}"
            rows.append(x)

        selected = pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()
        used = set(selected["accession"]) if not selected.empty else set()

        rec = rec[~rec["accession"].isin(used)].copy()
        rec_selected = take_evenly_across_qc(
            rec, n_recomb, prefer_explicit_recombinant=True
        )
        rec_selected["example_reason"] = "recombinant_or_discordant"

    else:  # random
        selected = sample_random(nonrec, n_nonrec)
        selected["example_reason"] = "nonrecombinant_random"

        used = set(selected["accession"]) if not selected.empty else set()
        rec = rec[~rec["accession"].isin(used)].copy()
        rec_selected = sample_random(rec, n_recomb)
        rec_selected["example_reason"] = "recombinant_or_discordant_random"

    out = pd.concat([selected, rec_selected], ignore_index=True)

    cols = [
        "accession",
        "example_reason",
        "clade_nextclade",
        "genotype_genbank",
        "subgenotype_genbank",
        "QC_overall_status",
        "QC_overall_score",
        "country",
        "date",
    ]
    cols = [c for c in cols if c in out.columns]

    out[cols].to_csv(args.output, sep="\t", index=False)


if __name__ == "__main__":
    main()
