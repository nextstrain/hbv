# This script lifts genomic region breakpoints from an old reference to a new reference using pairwise alignment and writes a segments file based on the lifted breakpoints.
#!/usr/bin/env python3
import argparse
from Bio import SeqIO
from Bio import pairwise2

def read_ref_seq_from_genbank(gb_path: str) -> str:
    rec = SeqIO.read(gb_path, "genbank")
    return str(rec.seq).upper()

def read_old_ref_sequence(old_ref: str) -> str:
    """
    Accepts either:
      - a path to a FASTA/GenBank file, or
      - an accession (e.g. X02763.1) and tries NCBI Entrez fetch.
    """
    # 1) treat as file path if it exists / has an extension we can parse
    for fmt in ("genbank", "fasta"):
        try:
            rec = SeqIO.read(old_ref, fmt)
            return str(rec.seq).upper()
        except Exception:
            pass

    # 2) otherwise treat as accession and try Entrez
    try:
        from Bio import Entrez  # only used if needed
        Entrez.email = "hello@nextstrain.org"
        with Entrez.efetch(db="nucleotide", id=old_ref, rettype="fasta", retmode="text") as h:
            rec = SeqIO.read(h, "fasta")
        return str(rec.seq).upper()
    except Exception as e:
        raise RuntimeError(
            f"Could not read old reference from file or fetch accession '{old_ref}'. "
            f"Provide a local FASTA/GenBank file path instead. Underlying error: {e}"
        )

def build_src_to_tgt_map(src_seq: str, tgt_seq: str) -> dict[int, int | None]:
    aln = pairwise2.align.globalms(src_seq, tgt_seq, 2, -1, -5, -1, one_alignment_only=True)[0]
    src_aln, tgt_aln = aln.seqA, aln.seqB

    src_pos = 0
    tgt_pos = 0
    src2tgt: dict[int, int | None] = {}


    for a, b in zip(src_aln, tgt_aln):
        if a != "-":
            src_pos += 1
        if b != "-":
            tgt_pos += 1
        if a != "-":
            src2tgt[src_pos] = (tgt_pos if b != "-" else None)

    return src2tgt

def lift_breakpoint(bp: int, src2tgt: dict[int, int | None], max_scan: int = 500) -> int:
    """
    Lift a breakpoint position. If it lands in a gap (None), scan outward for nearest mappable pos.
    """
    if bp in src2tgt and src2tgt[bp] is not None:
        return src2tgt[bp]  # type: ignore[return-value]

    for d in range(1, max_scan + 1):
        left = bp - d
        right = bp + d
        if left >= 1 and left in src2tgt and src2tgt[left] is not None:
            return src2tgt[left]  # type: ignore[return-value]
        if right in src2tgt and src2tgt[right] is not None:
            return src2tgt[right]  # type: ignore[return-value]

    raise ValueError(f"Breakpoint {bp} could not be lifted (gap/indel region too large).")



def shrink_interval(start, end, buf, L):
    # start/end are 1-based inclusive; may wrap if start > end
    if buf == 0:
        return start, end

    if start <= end:
        # linear
        s2 = start + buf
        e2 = end - buf
        if s2 > e2:
            raise ValueError(f"segment became invalid after buffer={buf}: {start}-{end} -> {s2}-{e2}")
        return s2, e2
    else:
        # wrap: [start..L] + [1..end]
        tail_len = L - start + 1
        head_len = end
        total = tail_len + head_len

        if 2 * buf >= total:
            raise ValueError(f"wrap segment too short for buffer={buf}: {start}-{end} (len {total})")

        # shrink both ends (stay wrap)
        s2 = start + buf
        e2 = end - buf

        # if we shrink past the origin, normalize by wrapping around
        if s2 > L:
            s2 = s2 - L
        if e2 < 1:
            e2 = e2 + L

        return s2, e2


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ref-gb", required=True, help="Target reference GenBank (e.g. ../ingest/defaults/NC_003977.gb)")
    ap.add_argument("--old-ref", required=True, help="Old reference accession (e.g. X02763.1) OR local fasta/genbank path")
    ap.add_argument("--breakpoints", required=True, nargs="+", type=int, help="Breakpoint list on old reference (1-based)")
    ap.add_argument("--buffer", type=int, default=0, help="Shrink each segment on both sides by this many bp")
    ap.add_argument("--out", required=True, help="Output txt path (segments file)")
    ap.add_argument("--map-out", default=None,
                help="Optional path to write source→target coordinate mapping (tab-separated)")

    args = ap.parse_args()

    buf = args.buffer
    if buf < 0:
        raise ValueError("--buffer must be >= 0")

    tgt_seq = read_ref_seq_from_genbank(args.ref_gb)
    src_seq = read_old_ref_sequence(args.old_ref)
    L = len(tgt_seq)

    src2tgt = build_src_to_tgt_map(src_seq, tgt_seq)
    if args.map_out:
        with open(args.map_out, "w") as fh:
            fh.write("# src_pos\ttgt_pos\n")
            for src_pos in sorted(src2tgt):
                tgt_pos = src2tgt[src_pos]
                fh.write(f"{src_pos}\t{tgt_pos if tgt_pos is not None else 'NA'}\n")


    bps_src = sorted(set(args.breakpoints))
    if not bps_src:
        raise ValueError("No breakpoints provided")
    if bps_src[0] != 1:
        bps_src = [1] + bps_src

    # lift breakpoints to target reference
    lifted = [lift_breakpoint(bp, src2tgt) for bp in bps_src]

    # unique + sorted
    bps_tgt = sorted(set(lifted))
    if bps_tgt[0] != 1:
        bps_tgt = [1] + bps_tgt

    with open(args.out, "w") as fh:
        fh.write("# start\tend\tname (1-based; inclusive; start>end means wrap)\n")

        # segments between consecutive breakpoints
        seg_idx = 1
        for i in range(len(bps_tgt) - 1):
            start = bps_tgt[i]
            end = bps_tgt[i + 1] - 1
            if end == 0:
                end = L
            start2, end2 = shrink_interval(start, end, buf, L)
            fh.write(f"{start2}\t{end2}\tsegment{seg_idx}\n")
            seg_idx += 1

        # closing circular segment: last -> first-1
        start = bps_tgt[-1]
        end = bps_tgt[0] - 1
        if end == 0:
            end = L
        start2, end2 = shrink_interval(start, end, buf, L)
        fh.write(f"{start2}\t{end2}\tsegment{seg_idx}\n")

if __name__ == "__main__":
    main()
