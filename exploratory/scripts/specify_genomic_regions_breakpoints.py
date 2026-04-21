#!/usr/bin/env python3
"""Lift genomic breakpoints from an old reference to a new reference.

Verify that the lifted breakpoints behave as intended before using them.
"""

import argparse
from Bio import SeqIO
from Bio.Align import PairwiseAligner


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
        with Entrez.efetch(
            db="nucleotide", id=old_ref, rettype="fasta", retmode="text"
        ) as h:
            rec = SeqIO.read(h, "fasta")
        return str(rec.seq).upper()
    except Exception as e:
        raise RuntimeError(
            f"Could not read old reference from file or fetch accession '{old_ref}'. "
            f"Provide a local FASTA/GenBank file path instead. Underlying error: {e}"
        )


def build_src_to_tgt_map(src_seq: str, tgt_seq: str) -> dict[int, int | None]:
    """
    Returns:
      cols[col] = (tgt_pos, src_pos)   # both 0-based, None for gaps
      src_pos is folded into [0, L)
    """
    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = 2
    aligner.mismatch_score = -1
    aligner.open_gap_score = -5
    aligner.extend_gap_score = -1

    circ = src_seq + src_seq
    aln = aligner.align(tgt_seq, circ)[0]

    coords = aln.coordinates
    t = coords[0]
    q = coords[1]
    L = len(src_seq)

    cols = []
    for j in range(len(t) - 1):
        t0, t1 = int(t[j]), int(t[j + 1])
        q0, q1 = int(q[j]), int(q[j + 1])

        dt = t1 - t0
        dq = q1 - q0
        steps = max(dt, dq)

        for s in range(steps):
            tgt_pos = (t0 + s) if dt else None

            if dq:
                circ_pos = q0 + s  # 0-based in src+src
                src_pos = circ_pos % L  # fold back to original src
            else:
                src_pos = None
            cols.append((tgt_pos, src_pos))
    return cols


def lift_breakpoint_from_cols(x, cols, L, max_scan=500):
    """
    x: source position (0-based, in [0, L))
    cols: list[(tgt_pos, src_pos)] with src_pos folded to [0, L)
    Returns tgt_pos (0-based)
    """
    # build an index once for speed: src_pos -> list of tgt_pos seen in alignment
    src2tgt = [[] for _ in range(L)]
    for tgt_pos, src_pos in cols:
        if src_pos is not None and tgt_pos is not None:
            src2tgt[src_pos].append(tgt_pos)

    if src2tgt[x]:
        return src2tgt[x][0]

    for d in range(1, max_scan + 1):
        left = (x - d) % L
        right = (x + d) % L

        if src2tgt[left]:
            return src2tgt[left][0]
        if src2tgt[right]:
            return src2tgt[right][0]

    raise ValueError(
        f"Breakpoint {x} could not be lifted (gap/indel region too large)."
    )


def shrink_interval(start, end, buf, L):
    # start/end are 1-based inclusive; may wrap if start > end
    if buf == 0:
        return start, end

    if start <= end:
        # linear
        s2 = start + buf
        e2 = end - buf
        if s2 > e2:
            raise ValueError(
                f"segment became invalid after buffer={buf}: {start}-{end} -> {s2}-{e2}"
            )
        return s2, e2
    else:
        # wrap: [start..L] + [1..end]
        tail_len = L - start + 1
        head_len = end
        total = tail_len + head_len

        if 2 * buf >= total:
            raise ValueError(
                f"wrap segment too short for buffer={buf}: {start}-{end} (len {total})"
            )

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
    ap.add_argument("--ref-gb", required=True, help="Target reference GenBank")
    ap.add_argument(
        "--old-ref",
        required=True,
        help="Old reference accession OR local fasta/genbank path",
    )
    ap.add_argument(
        "--breakpoints",
        required=True,
        nargs="+",
        type=int,
        help="Breakpoint list on old reference (1-based)",
    )
    ap.add_argument(
        "--buffer",
        type=int,
        default=0,
        help="Shrink each segment on both sides by this many bp",
    )
    ap.add_argument("--out", required=True, help="Output txt path (segments file)")
    ap.add_argument(
        "--map-out",
        default=None,
        help="Optional path to write alignment-column index (tab-separated)",
    )
    args = ap.parse_args()

    buf = args.buffer
    if buf < 0:
        raise ValueError("--buffer must be >= 0")

    tgt_seq = read_ref_seq_from_genbank(args.ref_gb)
    src_seq = read_old_ref_sequence(args.old_ref)

    L_tgt = len(tgt_seq)
    L_src = len(src_seq)

    bps_src_1based = sorted(set(args.breakpoints))
    if not bps_src_1based:
        raise ValueError("No breakpoints provided")
    if 1 not in bps_src_1based and L_src not in bps_src_1based:
        print(
            "No breakpoint at position 1 or at sequence end; "
            "assuming circular coverage.",
        )

    # Build alignment-column index using your new code
    # Expectation: build_src_to_tgt_map returns `cols` where:
    #   cols[col] = (tgt_pos_0based_or_None, src_pos_0based_or_None)   (src folded to [0, L_src))
    cols = build_src_to_tgt_map(src_seq, tgt_seq)

    if args.map_out:
        with open(args.map_out, "w") as fh:
            fh.write("# aln_col\ttgt_pos0\tsrc_pos0\n")
            for aln_col, (tgt_pos, src_pos) in enumerate(cols):
                fh.write(
                    f"{aln_col}\t"
                    f"{tgt_pos if tgt_pos is not None else 'NA'}\t"
                    f"{src_pos if src_pos is not None else 'NA'}\n"
                )

    # Lift breakpoints: old-ref (src) -> target (tgt)
    # breakpoints are 1-based in src; lift_breakpoint_from_cols expects 0-based in [0, L_src)
    lifted_tgt_0based = [
        lift_breakpoint_from_cols(bp - 1, cols, L_src) for bp in bps_src_1based
    ]
    bps_tgt_1based = sorted(set(p + 1 for p in lifted_tgt_0based))  # back to 1-based

    with open(args.out, "w") as fh:
        fh.write(
            "# name\tstart\tend\t(1-based, inclusive; start>end means wrap-around)\n"
        )

        k = len(bps_tgt_1based)
        for seg_idx in range(k):
            start = bps_tgt_1based[seg_idx]
            nxt = bps_tgt_1based[(seg_idx + 1) % k]
            end = nxt - 1
            if end == 0:
                end = L_tgt

            start2, end2 = shrink_interval(start, end, buf, L_tgt)
            fh.write(f"segment{seg_idx + 1}\t{start2}\t{end2}\n")


if __name__ == "__main__":
    main()
