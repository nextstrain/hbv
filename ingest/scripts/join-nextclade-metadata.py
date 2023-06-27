#!/usr/bin/env python3

# Based on https://github.com/nextstrain/ncov-ingest/blob/c77c39cf18ace98631834e10047ac2e09101b504/bin/join-metadata-and-clades

import argparse
import pandas as pd
import sys

METADATA_JOIN_COLUMN_NAME = 'name'
NEXTCLADE_JOIN_COLUMN_NAME = 'seqName'
VALUE_MISSING_DATA = '?'

column_map = {
    "clade": "clade_nextclade",
    "coverage": "coverage",
    "qc.overallScore": "QC_overall_score",
    "qc.overallStatus": "QC_overall_status",
    "qc.missingData.status": "QC_missing_data",
    "qc.mixedSites.status": "QC_mixed_sites",
    "qc.privateMutations.status": "QC_rare_mutations",
    "qc.frameShifts.status": "QC_frame_shifts",
    "qc.stopCodons.status": "QC_stop_codons",
    "alignmentScore": "alignment_score",
    "totalSubstitutions": "total_substitutions",
    "totalDeletions": "total_deletions",
    "totalInsertions": "total_insertions",
    "totalFrameShifts": "total_frame_shifts",
    "totalMissing": "total_missing",
}

def print_summary(name, n, total, fh=sys.stdout):
    stars = "*" * int(round(n/total*20))
    print(f"{name:15}{stars:20} {round(n/total*100, 1):5}% (n={n}/{total})", file=fh)

def summarise(series, title, empty='', fh=sys.stdout):
    print(title, file=fh)
    n = series.values.sum()
    for key in series.keys():
        print_summary(key if key else empty, series[key], n, fh)
    print(file=fh)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--metadata", required=True)
    parser.add_argument("--nextclade", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--summary", required=True)
    args = parser.parse_args()

    metadata = pd.read_csv(args.metadata, index_col=METADATA_JOIN_COLUMN_NAME,
                           sep='\t', low_memory=False, na_filter = False)

    # Read and rename clade column to be more descriptive
    clades = pd.read_csv(args.nextclade, index_col=NEXTCLADE_JOIN_COLUMN_NAME,
                         usecols=[NEXTCLADE_JOIN_COLUMN_NAME, *(set(column_map.keys()))],
                         sep='\t', low_memory=True, dtype="object", na_filter = False) \
            .rename(columns=column_map)

    # Concatenate on columns
    result = pd.merge(
        metadata, clades,
        left_index=True,
        right_index=True,
        how='left'
    )

    result.to_csv(args.output, index_label=METADATA_JOIN_COLUMN_NAME, sep='\t')

    ## Summarise statistics & print to screen (as well as summary file)
    with open(args.summary, 'w') as fh:
        
        qc = result.groupby(["QC_overall_status"])["QC_overall_status"].count()
        summarise(qc, "Overall Nextclade QC status", empty="Failed align", fh=fh)

        # drop the un-aligned sequences from remaining summaries
        aligned = result.drop(result[result["QC_overall_status"]==''].index, inplace=False)

        genotypes = aligned.groupby(["clade_nextclade"])["QC_overall_status"].count()
        summarise(genotypes, "Genotypes assigned via Nextclade inference (all aligned sequences)", fh=fh)

        genotypes = aligned.groupby(["genotype_genbank"])["genotype_genbank"].count()
        summarise(genotypes, "Genotypes assigned via GenBank annotation (all aligned sequences)", fh=fh)

        # compare genbank assigned genotypes to nextclade assigned ones
        print("Accuracy of inferred genotypes (nextclade) vs metadata-assigned genotypes (GenBank)", file=fh)
        g = aligned.groupby(["genotype_genbank", "clade_nextclade"])[["genotype_genbank", "clade_nextclade"]].size().reset_index()
        genotypes = set(g['genotype_genbank'])
        genotypes.remove('None')
        for genotype in sorted([*genotypes]):
            df = g[g['genotype_genbank']==genotype]
            print_summary(genotype, g[(g['genotype_genbank']==genotype) & (g['clade_nextclade']==genotype)][0].sum(), g[g['genotype_genbank']==genotype][0].sum(), fh)
