from Bio import Phylo, SeqRecord, Seq, SeqIO
import numpy as np
import glob, os
from treetime import TreeAnc

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Make a global tree")
    parser.add_argument("--tree-dir", type=str, required=True, help="Directory of trees to combine")
    parser.add_argument("--alignment-dir", type=str, required=True, help="Directory of alignment files")
    parser.add_argument("--output", type=str, required=True, help="output file")
    parser.add_argument("--outdir", type=str, required=True, help="output directory")
    parser.add_argument("--ref-id", type=str, required=True, help="reference accession to remove before stitching")
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    REF_ID = args.ref_id

    root_sequences = {}  # subtype -> Seq
    subtrees = {}        # subtype -> Bio.Phylo tree

    for t in glob.glob(os.path.join(args.tree_dir, "*.nwk")):
        subtype = os.path.basename(t).replace(".nwk", "")
        aln_fname = os.path.join(args.alignment_dir, subtype + ".fasta")

        sub_tree = Phylo.read(t, "newick")

        # prune duplicate reference tips (and any ref at all) from the subtree
        if REF_ID:
            while True:
                hits = [tip for tip in sub_tree.get_terminals() if tip.name == REF_ID]
                if not hits:
                    break
                sub_tree.prune(hits[0])

        sub_tree.root_at_midpoint()

        # write a temp alignment without the reference (if requested)
        aln_for_tt = aln_fname
        if REF_ID:
            tmp_aln = os.path.join(args.outdir, f"tmp_{subtype}.fasta")
            with open(tmp_aln, "w") as fh:
                for rec in SeqIO.parse(aln_fname, "fasta"):
                    if rec.id != REF_ID:
                        SeqIO.write(rec, fh, "fasta")
            aln_for_tt = tmp_aln

        tt = TreeAnc(sub_tree, aln=aln_for_tt, gtr="JC69")
        tt.infer_ancestral_sequences(marginal=True)

        ambig = np.where(np.sum(tt.tree.root.marginal_profile**2, axis=1) < 0.8)[0]
        root_array = tt.sequence(tt.tree.root, as_string=False)
        for pos in ambig:
            root_array[tt.data.compressed_to_full_sequence_map[pos]] = "N"

        root_sequences[subtype] = Seq.Seq("".join(root_array))

        # prefix internal node names to avoid collisions across subtrees
        for n in sub_tree.get_nonterminals():
            if n.name is None:
                n.name = "NODE"
            n.name = f"{subtype}_{n.name}"

        subtrees[subtype] = sub_tree

    # build root-sequence alignment (IDs are subtype names, not file paths)
    root_fasta = os.path.join(args.outdir, "root_seq.fasta")
    SeqIO.write(
        [SeqRecord.SeqRecord(root_sequences[s], id=s, description="") for s in sorted(root_sequences)],
        root_fasta,
        "fasta",
    )

    root_tree_raw = os.path.join(args.outdir, "root_tree_raw.nwk")
    os.system(f'augur tree --alignment "{root_fasta}" --output "{root_tree_raw}"')

    T = Phylo.read(root_tree_raw, "newick")
    T.root_at_midpoint()

    # graft each subtree onto the corresponding tip in the root tree
    for n in T.get_terminals():
        st = subtrees.get(n.name)
        if st is None:
            raise KeyError(f"Tip '{n.name}' in root tree not found among subtrees: {sorted(subtrees)[:5]} ...")
        n.clades = st.root.clades

    Phylo.write(T, args.output, "newick")
