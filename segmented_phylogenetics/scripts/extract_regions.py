# This script extracts subsequences from a multiple sequence alignment based on specified regions
# in a reference sequence, and also writes the corresponding sliced GenBank record (features clipped).
import sys
import os
from Bio import AlignIO, SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation

aln_path = sys.argv[1]
regions_path = sys.argv[2]
ref_id = sys.argv[3]
gb_path = sys.argv[4]
outdir = sys.argv[5] if len(sys.argv) > 5 else "."

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
    raise ValueError(
        f"Reference ID '{ref_id}' not found in alignment. "
        f"First few IDs: {[r.id for r in aln[:10]]}"
    )

ref = pick_ref_record(aln, ref_id)
ref_seq = str(ref.seq)

# --- build mapping: reference ungapped pos (1-based) -> alignment column (0-based) ---
refpos_to_col = {}
pos = 0
for col, c in enumerate(ref_seq):
    if c != "-":
        pos += 1
        refpos_to_col[pos] = col

aln_ref_len = pos
gb_ref_len = len(gb_ref.seq)

if aln_ref_len != gb_ref_len:
    raise ValueError(
        f"Reference length mismatch:\n"
        f"  alignment reference ({ref.id}): {aln_ref_len} bp\n"
        f"  GenBank reference ({gb_ref.id}): {gb_ref_len} bp\n"
        f"This means the alignment reference and GenBank are not the same origin or not the same reference."
    )

def clip_features(gb_ref, seg0, seg1, out_offset):
    """
    Keep the portions of features overlapping [seg0, seg1) on the original reference (0-based, half-open),
    and shift them into the output coordinate system by out_offset (0-based).
    """
    out = []
    for feat in gb_ref.features:
        if feat.location is None:
            continue

        f0 = int(feat.location.start)
        f1 = int(feat.location.end)

        if f1 <= seg0 or f0 >= seg1:
            continue

        new0 = max(f0, seg0) - seg0 + out_offset
        new1 = min(f1, seg1) - seg0 + out_offset

        out.append(
            SeqFeature(
                location=FeatureLocation(new0, new1, strand=feat.location.strand),
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
            sub_aln = aln[:, col_start:col_end + 1]
        else:
            # [start..L] + [1..end]
            col_start1 = refpos_to_col[start]
            col_end1 = refpos_to_col[gb_ref_len]
            col_start2 = refpos_to_col[1]
            col_end2 = refpos_to_col[end]
            sub_aln = aln[:, col_start1:col_end1 + 1] + aln[:, col_start2:col_end2 + 1]

        region_dir = f"{outdir}/{name}"
        os.makedirs(region_dir, exist_ok=True)
        AlignIO.write(sub_aln, f"{region_dir}/{name}_subsequence.fasta", "fasta")
        
        # --- GenBank sequence slice (wrap-aware) ---
        if not wraps:
            sub_seq = gb_ref.seq[start - 1:end]
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
            sub_features = (
                clip_features(gb_ref, segA0, segA1, 0)
                + clip_features(gb_ref, segB0, segB1, lenA)
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
