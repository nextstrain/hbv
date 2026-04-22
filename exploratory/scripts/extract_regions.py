"""Extract regional alignment slices and matching GenBank records.

The script writes per-region FASTA and GenBank outputs from reference-based coordinates.
"""

import os
from Bio import AlignIO, SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
import argparse

ap = argparse.ArgumentParser(
    description="Extract alignment regions by reference coordinates; write per-region FASTA, GenBank, and optional metadata."
)

ap.add_argument("aln_path", help="Input multiple sequence alignment (FASTA).")
ap.add_argument(
    "regions_path",
    help="Regions file: 'name start end' (1-based, inclusive; start>end means circular wrap).",
)
ap.add_argument(
    "ref_id", help="Reference sequence ID in the alignment (exact or prefix match)."
)
ap.add_argument("gb_path", help="GenBank file for the reference sequence.")
ap.add_argument(
    "outdir",
    nargs="?",
    default=".",
    help="Output directory (default: current directory).",
)

ap.add_argument(
    "--min-cov",
    type=float,
    default=0.0,
    help="Minimum non-gap fraction required to keep a sequence per region (0.0–1.0).",
)
ap.add_argument(
    "--metadata", default=None, help="Optional metadata TSV file to subset per region."
)
ap.add_argument(
    "--keep-cols",
    nargs="*",
    default=None,
    help="Metadata columns to keep (default: all).",
)

args = ap.parse_args()

if not (0.0 <= args.min_cov <= 1.0):
    ap.error("--min-cov must be between 0.0 and 1.0")

aln_path = args.aln_path
regions_path = args.regions_path
ref_id = args.ref_id
gb_path = args.gb_path
outdir = args.outdir
min_cov = args.min_cov
metadata = args.metadata
keep_cols = args.keep_cols

aln = AlignIO.read(aln_path, "fasta")
gb_ref = SeqIO.read(gb_path, "genbank")


# --- pick reference row in alignment ---
def pick_ref_record(aln, ref_id):
    # exact match first
    for r in aln:
        if r.id == ref_id:
            return r
    # fallback: tolerate extra suffixes (spaces, pipes, etc.)
    for r in aln:
        rid = r.id.split()[0]
        if rid == ref_id or rid.startswith(ref_id):
            return r
    return None


ref = pick_ref_record(aln, ref_id)
gb_ref_len = len(gb_ref.seq)
aln_len = aln.get_alignment_length()

if ref is not None:
    ref_seq = str(ref.seq)

    # --- build mapping: reference ungapped pos (1-based) -> alignment column (0-based) ---
    refpos_to_col = {}
    pos = 0
    for col, c in enumerate(ref_seq):
        if c != "-":
            pos += 1
            refpos_to_col[pos] = col

    aln_ref_len = pos

    if aln_ref_len != gb_ref_len:
        raise ValueError(
            f"Reference length mismatch:\n"
            f"  alignment reference ({ref.id}): {aln_ref_len} bp\n"
            f"  GenBank reference ({gb_ref.id}): {gb_ref_len} bp\n"
            f"This means the alignment reference and GenBank are not the same origin or not the same reference."
        )
else:
    if aln_len != gb_ref_len:
        raise ValueError(
            f"Reference ID '{ref_id}' not found in alignment, and the alignment length "
            f"({aln_len} bp) does not match the GenBank reference length ({gb_ref_len} bp).\n"
            f"First few IDs: {[r.id for r in aln[:10]]}"
        )
    refpos_to_col = {pos: pos - 1 for pos in range(1, gb_ref_len + 1)}


def clip_features(gb_ref, seg0, seg1, out_offset):
    """
    Clip features to overlap with [seg0, seg1) on original reference (0-based, half-open),
    then shift into output coords by out_offset (0-based).
    Works for FeatureLocation and CompoundLocation (join).
    """
    out = []
    for feat in gb_ref.features:
        if feat.location is None:
            continue

        loc = feat.location
        parts = loc.parts if hasattr(loc, "parts") else [loc]

        new_parts = []
        for p in parts:
            # intersection of [p.start,p.end) with [seg0,seg1)
            a0, a1 = int(p.start), int(p.end)
            b0, b1 = seg0, seg1
            i0, i1 = max(a0, b0), min(a1, b1)
            if i1 <= i0:
                continue

            new_parts.append(
                FeatureLocation(
                    i0 - seg0 + out_offset,
                    i1 - seg0 + out_offset,
                    strand=p.strand,
                )
            )

        if not new_parts:
            continue

        # If original was a join and we kept 2 parts, Biopython will keep it as CompoundLocation automatically
        new_loc = (
            new_parts[0] if len(new_parts) == 1 else sum(new_parts[1:], new_parts[0])
        )

        out.append(
            SeqFeature(
                location=new_loc,
                type=feat.type,
                qualifiers=feat.qualifiers,
            )
        )

    return out


with open(regions_path) as f:
    for line in f:
        line = line.strip()
        if not line or line.startswith("#"):
            continue

        name, start, end = line.split()
        start, end = int(start), int(end)

        if start < 1 or end < 1 or start > gb_ref_len or end > gb_ref_len:
            raise ValueError(
                f"{name}: invalid interval start={start}, end={end}, ref_len={gb_ref_len}"
            )

        wraps = start > end  # circular interval spanning end of genome

        # --- aligned FASTA slice (wrap-aware) ---
        if not wraps:
            col_start = refpos_to_col[start]
            col_end = refpos_to_col[end]
            sub_aln = aln[:, col_start : col_end + 1]
        else:
            # [start..L] + [1..end]
            col_start1 = refpos_to_col[start]
            col_end1 = refpos_to_col[gb_ref_len]
            col_start2 = refpos_to_col[1]
            col_end2 = refpos_to_col[end]
            sub_aln = (
                aln[:, col_start1 : col_end1 + 1] + aln[:, col_start2 : col_end2 + 1]
            )

        # --- coverage filter (per sequence in this region) ---
        if min_cov > 0.0:
            print(
                f"Applying minimum coverage filter of {min_cov * 100}% for region {name} ({start}-{end})",
                end=": ",
            )
        region_len = sub_aln.get_alignment_length()
        kept = []
        for rec in sub_aln:
            s = str(rec.seq)
            cov = (region_len - s.count("-")) / region_len
            if cov >= min_cov:
                kept.append(rec)
        print(
            "Keeping {}/{} sequences for region {}".format(
                len(kept), len(sub_aln), name
            )
        )

        sub_aln = sub_aln.__class__(kept)  # MultipleSeqAlignment from kept records

        region_dir = f"{outdir}/{name}"
        os.makedirs(region_dir, exist_ok=True)
        AlignIO.write(sub_aln, f"{region_dir}/{name}_sub-alignment.fasta", "fasta")

        if metadata is not None:
            kept_ids = {r.id for r in sub_aln}
            subseq_len_by_id = {r.id: len(str(r.seq).replace("-", "")) for r in sub_aln}

            with open(metadata) as f:
                header = f.readline().rstrip("\n").split("\t")

                id_col = header[0]  # first column is the ID
                keep = (
                    header
                    if keep_cols is None
                    else [id_col] + [c for c in keep_cols if c in header]
                )
                keep_out = ["length_subsequence" if c == "length" else c for c in keep]

                idx = {c: header.index(c) for c in keep}

                out_meta = f"{region_dir}/{name}_metadata.tsv"
                with open(out_meta, "w") as out_fh:
                    out_fh.write("\t".join(keep_out) + "\n")
                    for line in f:
                        row = line.rstrip("\n").split("\t")
                        if not row:
                            continue
                        if row[0] in kept_ids:
                            seq_id = row[0]
                            out_fh.write(
                                "\t".join(
                                    str(subseq_len_by_id[seq_id])
                                    if c == "length"
                                    else row[idx[c]]
                                    for c in keep
                                )
                                + "\n"
                            )

        # --- GenBank sequence slice (wrap-aware) ---
        if not wraps:
            sub_seq = gb_ref.seq[start - 1 : end]
            seg0, seg1 = start - 1, end  # 0-based, half-open
            sub_features = clip_features(gb_ref, seg0, seg1, 0)
        else:
            # segment A: [start-1, L), segment B: [0, end)
            segA0, segA1 = start - 1, gb_ref_len
            segB0, segB1 = 0, end

            seqA = gb_ref.seq[segA0:segA1]
            seqB = gb_ref.seq[segB0:segB1]
            lenA = len(seqA)

            sub_seq = seqA + seqB
            sub_features = clip_features(gb_ref, segA0, segA1, 0) + clip_features(
                gb_ref, segB0, segB1, lenA
            )

        # optional: keep exactly one source feature spanning the extracted record
        sub_features = [ft for ft in sub_features if ft.type != "source"]
        sub_features.insert(
            0,
            SeqFeature(
                FeatureLocation(0, len(sub_seq)),
                type="source",
                qualifiers=(
                    gb_ref.features[0].qualifiers
                    if gb_ref.features and gb_ref.features[0].type == "source"
                    else {}
                ),
            ),
        )

        rec = SeqRecord(
            sub_seq,
            id=f"{gb_ref.id}_{name}",
            name=name,
            description=f"{name} {start}-{end} (circular; ref {gb_ref.id})",
        )
        rec.annotations = gb_ref.annotations
        rec.features = sub_features

        SeqIO.write(rec, f"{region_dir}/{name}.gb", "genbank")
