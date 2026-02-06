from Bio import Phylo, SeqRecord, Seq, SeqIO
import numpy as np
import glob, os
from treetime import TreeAnc

if __name__=="__main__":
    import argparse
    parser = argparse.ArgumentParser(description='Make a global tree')
    parser.add_argument('--tree-dir', type=str, help='Directory of trees to combine')
    parser.add_argument('--alignment-dir', type=str, help='Directory of alignment files')
    parser.add_argument('--output', type=str, help='output file')
    args = parser.parse_args()

    root_sequences = {}
    subtrees = {}
    for t in glob.glob(args.tree_dir + "/*.nwk"):
        aln_fname = args.alignment_dir + "/" + t.split("/")[-1].replace(".nwk", ".fasta")
        sub_tree = Phylo.read(t, "newick")
        sub_tree.root_at_midpoint()
        tt = TreeAnc(sub_tree, aln=aln_fname, gtr="JC69")
        tt.infer_ancestral_sequences(marginal=True)
        ambig = np.where(np.sum(tt.tree.root.marginal_profile**2, axis=1)<0.8)[0]
        root_array = tt.sequence(tt.tree.root, as_string=False)
        for pos in ambig:
            root_array[tt.data.compressed_to_full_sequence_map[pos]] = 'N'
        root_sequences[t] = Seq.Seq(''.join(root_array))

        subtype = t.split("/")[-1].replace(".nwk", "")
        for n in sub_tree.get_nonterminals():
            n.name = subtype + "_" + n.name
        subtrees[t] = sub_tree


    SeqIO.write([SeqRecord.SeqRecord(root_sequences[t], id=t) for t in root_sequences], "results/root_seq.fasta", "fasta")
    os.system("augur tree --alignment results/root_seq.fasta --output results/root_tree_raw.nwk")

    T = Phylo.read("results/root_tree_raw.nwk", "newick")
    T.root_at_midpoint()
    for n in T.get_terminals():
        n.clades = subtrees[n.name].root.clades

    # for t, sub_tree in subtrees.items():
    #     sub_tree.root.branch_length = 0.1
    #     T.root.clades.append(sub_tree.root)

    # T = Phylo.BaseTree.Tree()
    #     sub_tree.root.branch_length = 0.1
    #     T.root.clades.append(sub_tree.root)
    Phylo.write(T, args.output, "newick")
