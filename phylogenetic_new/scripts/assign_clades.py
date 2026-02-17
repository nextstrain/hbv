from Bio import Phylo
import pandas as pd
import argparse
from collections import defaultdict
import json


def is_missing(v):
    return ( v in ("other", "NA", "") or v is None or (isinstance(v, float) and pd.isna(v) ))

def monophyletic_label_clades(tree, label, min_count, node_to_label):
    # annotate counts bottom-up
    for cl in tree.find_clades(order="postorder"):
        if cl.is_terminal():
            is_lab = (node_to_label.get(cl.name) == label)
            cl._n_label = 1 if is_lab else 0
            cl._n_other = 0 if is_lab else 1
        else:
            cl._n_label = sum(getattr(ch, "_n_label", 0) for ch in cl.clades)
            cl._n_other = sum(getattr(ch, "_n_other", 0) for ch in cl.clades)

    out = []
    for cl in tree.find_clades(order="preorder"):
        if cl._n_label > 0 and cl._n_other == 0 and cl._n_label >= min_count:
            parent = tree.get_path(cl)[-2] if len(tree.get_path(cl)) >= 2 else None
            parent_pure = (parent is not None and parent._n_label > 0 and parent._n_other == 0)
            if not parent_pure:
                out.append(cl)

    # cleanup (optional)
    for cl in tree.find_clades():
        if hasattr(cl, "_n_label"): del cl._n_label
        if hasattr(cl, "_n_other"): del cl._n_other

    return out



if __name__=="__main__":
    parser = argparse.ArgumentParser(description='Assign clades to a tree')
    parser.add_argument('--tree', type=str, help='Newick tree file')
    parser.add_argument('--metadata', type=str, help='metadata file with subtype information')
    parser.add_argument('--output', type=str, help='output file')
    parser.add_argument('--min-count', type=int, help='subtypes to ignore if they have fewer than this many samples', default=5)
    parser.add_argument('--subtype-col', type=str, help='column in metadata file to use for subtype information')
    parser.add_argument('--nextclade-label', type=str, help='label to use for nextclade clades')

    parser.add_argument('--fallback-genotype-col', type=str, help='column in metadata file to use as clade if subtype information is missing', required=False)
    parser.add_argument('--fallback-annotation-file', type=str, help='file with fallback annotations', required=False)
    parser.add_argument('--fallback-annotation-nextclade-label', type=str, help='label to use for nextclade clades in fallback annotation file', required=False)

    parser.add_argument('--branch_display', type=str, help='whether to enforce monophyly when assigning branch labels (options: "display_once", "display_all", "manually_curated")', default="display_once")
    parser.add_argument('--only_monophyletic', type="bool", help='only store monophyletic clades; split labels into multiple monophyletic clades if needed')

    
    args = parser.parse_args()

    tree = Phylo.read(args.tree, 'newick')
    metadata = pd.read_csv(args.metadata, sep='\t', index_col=0)

    node_to_subtype = {x.Index: getattr(x, args.subtype_col) for x in metadata.itertuples()}   # not sure thats genotype instead of
    nodes_by_subtype_all = defaultdict(list)

    for node in tree.get_terminals():
        nodes_by_subtype_all[node_to_subtype[node.name]].append(node)

    nodes_by_subtype = defaultdict(list)
    for subtype, nodes in nodes_by_subtype_all.items():
        if len(nodes) >= args.min_count:
            nodes_by_subtype[subtype] = nodes
        else:
            nodes_by_subtype['other'].extend(nodes)

    for node in nodes_by_subtype['other']:
        node_to_subtype[node.name] = 'other'

    clade_demarcations = {}
    for subtype, nodes in nodes_by_subtype.items():
        clade = tree.is_monophyletic(nodes)
        if clade:
            clade_demarcations[subtype] = clade
        else:
            clade_demarcations[subtype] = tree.common_ancestor(nodes)
            print("Subtype {} is not monophyletic".format(subtype))
        clade_demarcations[subtype].clade_label = subtype

    for node in tree.find_clades(order = 'preorder'):
        node.clade = 'other'

    for node in tree.find_clades(order = 'preorder'):
        if hasattr(node, "clade_label"):
            node.clade = node.clade_label
            for child in node.find_clades():
                child.clade = node.clade_label

    fallback_annotations = None
    if args.fallback_annotation_file is not None:
        with open(args.fallback_annotation_file, "r") as f:
            fallback_annotations = json.load(f)

        
    node_data = {'nodes': {}, 'branches': {}}
    for node in tree.find_clades(order = 'preorder'):
        if node.is_terminal():
            if node.clade != node_to_subtype[node.name]:
                print("Node {} has clade {} but subtype {}".format(node.name, node.clade, node_to_subtype[node.name]))
                node.clade = node_to_subtype[node.name]
        
        datum = {args.nextclade_label: node.clade}
    
        # If the node is missing subtype information and a fallback annotation file is provided, use the fallback annotation for this node if available
        if fallback_annotations is not None and is_missing(node.clade):
            fallback_node = fallback_annotations.get("nodes", {}).get(node.name, {})
            fallback_value = fallback_node.get(args.fallback_annotation_nextclade_label)

            if fallback_value is not None:
                datum[args.nextclade_label] = fallback_value

        # TODO add branch settings functionality here
        if hasattr(node, 'clade_label'):
            node_data['branches'][node.name] = {"labels": {args.nextclade_label: node.clade_label}}

        if fallback_annotations is not None:
            fallback_branch = fallback_annotations.get("branches", {}).get(node.name, {})
            fallback_value = fallback_branch.get("labels", {}).get(args.fallback_annotation_nextclade_label)
            if fallback_value is not None:
                node_data['branches'][node.name] = {"labels": {args.nextclade_label: fallback_value}}

        node_data['nodes'][node.name] = datum


    with open(args.output, 'w') as f:
        json.dump(node_data, f, indent=2)