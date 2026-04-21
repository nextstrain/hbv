from Bio import Phylo
import pandas as pd
import argparse
from collections import defaultdict
import json

def is_missing(v):
    return ( v in ("other", "NA", "") or v is None or (isinstance(v, float) and pd.isna(v) ))

def monophyletic_label_clades(tree, label, node_to_label, args):
    # parent map (Biopython clades don't store parent pointers)
    parent_of = {}
    for p in tree.find_clades(order="preorder"):
        for ch in p.clades:
            parent_of[ch] = p

    # annotate counts bottom-up
    for cl in tree.find_clades(order="postorder"):
        if cl.is_terminal():
            is_lab = (node_to_label.get(cl.name) == label)
            cl._n_label = 1 if is_lab else 0
            cl._n_other = 0 if is_lab else 1
        else:
            cl._n_label = sum(ch._n_label for ch in cl.clades)
            cl._n_other = sum(ch._n_other for ch in cl.clades)

    out = []
    for cl in tree.find_clades(order="preorder"):
        if cl._n_label > 0 and cl._n_other == 0 :
            parent = parent_of.get(cl, None)
            parent_pure = (parent is not None and parent._n_label > 0 and parent._n_other == 0)
            if not parent_pure:
                if args.min_count_mode =="use_for_annotation" and args.min_count>cl._n_label:
                    for leaf in cl.get_terminals():
                        node_to_label[leaf.name] = "other"
                else:
                    out.append(cl)

    # cleanup
    for cl in tree.find_clades(order="preorder"):
        del cl._n_label
        del cl._n_other

    return out

def assign_clade_info(tree , node_data, args, fallback_annotations, node_to_subtype):
    for node in tree.find_clades(order = 'preorder'):
        node.clade = 'other'

    for node in tree.find_clades(order = 'preorder'):
        if hasattr(node, "clade_label"):
            node.clade = node.clade_label
            for child in node.find_clades():
                child.clade = node.clade_label

    for node in tree.find_clades(order = 'preorder'):
        if node.is_terminal():
            if node.clade != node_to_subtype[node.name]:
                node.clade = node_to_subtype[node.name]

        datum = {args.nextclade_label: node.clade}

        # If the node is missing subtype information and a fallback annotation file is provided, use the fallback annotation for this node if available
        if fallback_annotations is not None and is_missing(node.clade):
            fallback_node = fallback_annotations.get("nodes", {}).get(node.name, {})
            fallback_value = fallback_node.get(args.fallback_annotation_nextclade_label)

            if fallback_value is not None:
                datum[args.nextclade_label] = fallback_value

        node_data['nodes'][node.name] = datum
    return tree, node_data

def assign_branch_clades(tree, node_data, args, fallback_annotations=None):

    best_node_for_label = {}
    best_size_for_label = {}

    for n in tree.find_clades(order="preorder"):
        if hasattr(n, "clade_label") and n.name not in (None, ""):
            sz = len(n.get_terminals())
            lab = n.clade_label
            if sz > best_size_for_label.get(lab, -1):
                best_size_for_label[lab] = sz
                best_node_for_label[lab] = n

    for node in tree.find_clades(order = 'preorder'):
        if hasattr(node, 'clade_label') and node.name not in (None, ""):
            # Optional: Hide branch labels for very small clades
            if args.min_count_mode == "use_for_branch_display" and len(node.get_terminals())  < args.min_count:
                continue
            # Optional: Only display label for largest clade
            if args.only_monophyletic and args.branch_display_only_largest:
                if best_node_for_label.get(node.clade_label) is not node:
                    continue
            node_data['branches'][node.name] = {"labels": {args.nextclade_label: node.clade_label}}

        if fallback_annotations is not None:
            fallback_branch = fallback_annotations.get("branches", {}).get(node.name, {})
            fallback_value = fallback_branch.get("labels", {}).get(args.fallback_annotation_nextclade_label)
            if fallback_value is not None:
                node_data['branches'][node.name] = {"labels": {args.nextclade_label: fallback_value}}

    return tree, node_data

if __name__=="__main__":

    # ______Read in arguments and data_____________________________________________________________________________________________________
    parser = argparse.ArgumentParser(description='Assign clades to a tree')
    parser.add_argument('--tree', type=str, help='Newick tree file')
    parser.add_argument('--metadata', type=str, help='metadata file with subtype information')
    parser.add_argument('--output', type=str, help='output file')
    parser.add_argument('--min-count', type=int, help='subtypes to ignore if they have fewer than this many samples', default=5)
    parser.add_argument('--min-count-mode', type=str, help='mode to use for subtypes with fewer than min-count samples (options: "use_for_branch_display", "use_for_annotation")') # "use_for_pruning" config option has no effect here but is used in earlier scirpt
    parser.add_argument('--subtype-col', type=str, help='column in metadata file to use for subtype information')
    parser.add_argument('--nextclade-label', type=str, help='label to use for nextclade clades')

    parser.add_argument('--fallback-genotype-col', type=str, help='column in metadata file to use as clade if subtype information is missing', required=False)
    parser.add_argument('--fallback-annotation-file', type=str, help='file with fallback annotations', required=False)
    parser.add_argument('--fallback-annotation-nextclade-label', type=str, help='label to use for nextclade clades in fallback annotation file', required=False)

    parser.add_argument('--only_monophyletic', action="store_true", help='only store monophyletic clades; split labels into multiple monophyletic clades if needed')
    parser.add_argument('--branch_display_only_largest', action="store_true", help='whether to display branch labels multiple times on tree (if only-monophyletic is set) or only once for largest clade')

    args = parser.parse_args()

    print(
        f"Assigning clades based on "
        f"{'monophyly' if args.only_monophyletic else 'MRCA'} "
        f"in {args.subtype_col} metadata annotation and new nextclade label "
        f"{args.nextclade_label}.\n"
        f"Clades with fewer members than {args.min_count} are "
        f"{'removed from annotation.' if args.min_count_mode == 'use_for_annotation' else 'not displayed on branches.' if args.min_count_mode == 'use_for_branch_display' else 'have already been pruned.'}\n"
        f"{'For multiple monophyletic clades per subtype only the largest is displayed.\n' if args.branch_display_only_largest else '\n'}"
        )

    tree = Phylo.read(args.tree, 'newick')
    metadata = pd.read_csv(args.metadata, sep='\t', index_col=0)

    fallback_annotations=None
    if args.fallback_annotation_file is not None:
        with open(args.fallback_annotation_file, "r") as f:
            fallback_annotations = json.load(f)

    # ______Match nodes and subtypes_____________________________________________________________________________________________________
    node_to_subtype = {x.Index: getattr(x, args.subtype_col) for x in metadata.itertuples()}
    nodes_by_subtype_all = defaultdict(list)

    for node in tree.get_terminals():
        nodes_by_subtype_all[node_to_subtype[node.name]].append(node)

    nodes_by_subtype = defaultdict(list)

    # Optional: Make tree annotation cleaner by dropping subtype information for subtypes with fewer than min-count samples
    for subtype, nodes in nodes_by_subtype_all.items():
        if args.min_count_mode == "use_for_annotation":
            if len(nodes) >= args.min_count:
                nodes_by_subtype[subtype] = nodes
            else:
                nodes_by_subtype['other'].extend(nodes)
        else:
            nodes_by_subtype[subtype] = nodes

    for node in nodes_by_subtype['other']:
        node_to_subtype[node.name] = 'other'

    # ______Assign clades to internal nodes and branches_____________________________________________________________________________________________________
    for subtype, nodes in nodes_by_subtype.items():
        if subtype == "other":
            continue

        clade = tree.is_monophyletic(nodes)
        if not clade:
            print("Subtype {} is not monophyletic".format(subtype))
            if not args.only_monophyletic:
                clade = tree.common_ancestor(nodes)
            else:
                # split into maximal monophyletic clades for this subtype
                clades = monophyletic_label_clades(tree, subtype, node_to_subtype, args)
                for c in clades:
                    c.clade_label = subtype
                continue
        clade.clade_label = subtype

    node_data = {'nodes': {}, 'branches': {}}

    tree, node_data = assign_clade_info   (tree, node_data, args, fallback_annotations, node_to_subtype)
    tree, node_data = assign_branch_clades(tree, node_data, args, fallback_annotations)

    with open(args.output, 'w') as f:
        json.dump(node_data, f, indent=2)
