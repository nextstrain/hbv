from Bio import Phylo
import argparse

def prune_side_of_tree(tree, tips_to_exclude):
    """
    Preorder traversal until we find the large clade with none of the provided names
    """
    for root_node in tree.find_clades():
        tips = {n.name for n in root_node.get_terminals()}
        if len(tips)<100:
            continue
        if len(tips_to_exclude & tips)==0:
            return root_node
    raise Exception("Tree pruning algorithm failed")
    
if __name__=="__main__":
    parser = argparse.ArgumentParser(
        description="Remove outgroup sequences, and those which fall alongside them on that side of the tree",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('--tree', type=str, required=True)
    parser.add_argument('--output', type=str, required=True)
    parser.add_argument('--names', type=str,  nargs='+', required=True, help="sequence names to remove")
    args = parser.parse_args()
    tree = prune_side_of_tree(Phylo.read(args.tree, "newick"), set(args.names))
    Phylo.write(tree, args.output, "newick")
