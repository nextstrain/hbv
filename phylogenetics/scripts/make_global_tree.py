"""Stitch genotype-specific trees into one global tree using inferred subtree roots.

Remove duplicate reference tips, infer one root sequence per subtree, and graft the subtrees onto a root tree.
"""

import glob
import subprocess
import tempfile
from pathlib import Path

from Bio import Phylo, SeqRecord, Seq, SeqIO
import numpy as np
from treetime import TreeAnc

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Make a global tree")
    parser.add_argument(
        "--tree-dir", type=str, required=True, help="Directory of trees to combine"
    )
    parser.add_argument(
        "--alignment-dir", type=str, required=True, help="Directory of alignment files"
    )
    parser.add_argument("--output", type=str, required=True, help="output file")
    parser.add_argument(
        "--ref-id",
        type=str,
        required=True,
        help="reference accession to remove before stitching",
    )
    parser.add_argument(
        "--log",
        type=str,
        help="Optional destination for the temporary root-tree build log",
    )
    args = parser.parse_args()

    REF_ID = args.ref_id
    log_path = Path(args.log) if args.log else None

    root_sequences = {}  # subtype -> Seq
    subtrees = {}  # subtype -> Bio.Phylo tree

    with tempfile.TemporaryDirectory(prefix="make_global_tree_") as tempdir:
        tempdir = Path(tempdir)

        for t in glob.glob(str(Path(args.tree_dir) / "*.nwk")):
            subtype = Path(t).stem
            aln_fname = Path(args.alignment_dir) / f"{subtype}.fasta"

            sub_tree = Phylo.read(t, "newick")

            # prune duplicate reference tips (and any ref at all) from the subtree
            if REF_ID:
                while True:
                    hits = [
                        tip for tip in sub_tree.get_terminals() if tip.name == REF_ID
                    ]
                    if not hits:
                        break
                    sub_tree.prune(hits[0])

            sub_tree.root_at_midpoint()

            # TreeAnc needs a temporary alignment without the duplicate reference tip.
            aln_for_tt = str(aln_fname)
            if REF_ID:
                tmp_aln = tempdir / f"tmp_{subtype}.fasta"
                with open(tmp_aln, "w") as fh:
                    for rec in SeqIO.parse(aln_fname, "fasta"):
                        if rec.id != REF_ID:
                            SeqIO.write(rec, fh, "fasta")
                aln_for_tt = str(tmp_aln)

            tt = TreeAnc(sub_tree, aln=aln_for_tt, gtr="JC69")
            tt.infer_ancestral_sequences(marginal=True)

            ambig = np.where(np.sum(tt.tree.root.marginal_profile**2, axis=1) < 0.8)[0]
            root_array = tt.sequence(tt.tree.root, as_string=False)
            for pos in ambig:
                root_array[tt.data.compressed_to_full_sequence_map[pos]] = "N"

            root_sequences[subtype] = Seq.Seq("".join(root_array))

            # Prefix internal node names to avoid collisions across subtrees.
            for n in sub_tree.get_nonterminals():
                if n.name is None:
                    n.name = "NODE"
                n.name = f"{subtype}_{n.name}"

            subtrees[subtype] = sub_tree

        # Build the temporary root-sequence alignment from one inferred root per subtree.
        root_fasta = tempdir / "root_seq.fasta"
        SeqIO.write(
            [
                SeqRecord.SeqRecord(root_sequences[s], id=s, description="")
                for s in sorted(root_sequences)
            ],
            root_fasta,
            "fasta",
        )

        root_tree_raw = tempdir / "root_tree_raw.nwk"
        subprocess.run(
            [
                "augur",
                "tree",
                "--alignment",
                str(root_fasta),
                "--output",
                str(root_tree_raw),
            ],
            check=True,
        )

        if log_path:
            temp_logs = sorted(tempdir.glob("*.log"))
            if temp_logs:
                with open(log_path, "w") as out:
                    for i, temp_log in enumerate(temp_logs):
                        if i:
                            out.write("\n")
                        out.write(f"== {temp_log.name} ==\n")
                        out.write(temp_log.read_text())

        T = Phylo.read(root_tree_raw, "newick")
        T.root_at_midpoint()

        # Graft each subtree onto the corresponding tip in the root tree.
        for n in T.get_terminals():
            st = subtrees.get(n.name)
            if st is None:
                raise KeyError(
                    f"Tip '{n.name}' in root tree not found among subtrees: {sorted(subtrees)[:5]} ..."
                )
            n.clades = st.root.clades

        Phylo.write(T, args.output, "newick")
