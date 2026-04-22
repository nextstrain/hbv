"""Subsample exploratory tree inputs while keeping the reference record when available.

Write matching FASTA and metadata subsets for fast dev-mode comparison runs.
"""

import argparse
import random
from Bio import SeqIO


def parse_args():
    parser = argparse.ArgumentParser(
        description="Subset exploratory tree inputs while retaining the reference record when present."
    )
    parser.add_argument("--sequences", required=True, help="Input sequences FASTA.")
    parser.add_argument("--alignment", required=True, help="Input alignment FASTA.")
    parser.add_argument("--metadata", required=True, help="Input metadata TSV.")
    parser.add_argument(
        "--out-sequences", required=True, help="Output sequences FASTA."
    )
    parser.add_argument(
        "--out-alignment", required=True, help="Output alignment FASTA."
    )
    parser.add_argument("--out-metadata", required=True, help="Output metadata TSV.")
    parser.add_argument("--reference-id", required=True, help="Reference ID to retain.")
    parser.add_argument(
        "--sample-size",
        type=int,
        required=True,
        help="Total number of sequences to keep.",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=1,
        help="Random seed for subsampling.",
    )
    return parser.parse_args()


def matches_reference(record_id, reference_id):
    short_id = record_id.split()[0]
    return short_id == reference_id or short_id.startswith(reference_id)


def write_fasta_subset(in_path, out_path, keep_ids):
    records = [
        record for record in SeqIO.parse(in_path, "fasta") if record.id in keep_ids
    ]
    SeqIO.write(records, out_path, "fasta")
    return {record.id for record in records}


def main():
    args = parse_args()

    if args.sample_size < 1:
        raise SystemExit("--sample-size must be at least 1.")

    sequence_records = list(SeqIO.parse(args.sequences, "fasta"))
    reference_record = next(
        (
            record
            for record in sequence_records
            if matches_reference(record.id, args.reference_id)
        ),
        None,
    )
    rng = random.Random(args.seed)

    if reference_record is None:
        all_ids = [record.id for record in sequence_records]
        n_keep = min(len(all_ids), args.sample_size)
        keep_ids = set(rng.sample(all_ids, n_keep))
    else:
        other_ids = [
            record.id for record in sequence_records if record.id != reference_record.id
        ]
        n_other = min(len(other_ids), args.sample_size - 1)
        sampled_other_ids = set(rng.sample(other_ids, n_other))
        keep_ids = {reference_record.id} | sampled_other_ids

    kept_sequence_ids = write_fasta_subset(args.sequences, args.out_sequences, keep_ids)
    kept_alignment_ids = write_fasta_subset(
        args.alignment, args.out_alignment, keep_ids
    )

    if reference_record is not None and reference_record.id not in kept_alignment_ids:
        raise SystemExit(
            f"Reference ID '{reference_record.id}' was not found in {args.alignment}."
        )

    with open(args.metadata) as in_fh, open(args.out_metadata, "w") as out_fh:
        header = in_fh.readline()
        out_fh.write(header)
        for line in in_fh:
            if not line.strip():
                continue
            if line.split("\t", 1)[0] in kept_sequence_ids:
                out_fh.write(line)


if __name__ == "__main__":
    main()
