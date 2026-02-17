import argparse
from pathlib import Path
from ete3 import Tree
from numpy import size


def suppress_unary_nodes(t):
    for n in list(t.traverse()):
        if not n.is_root() and len(n.children) == 1:
            child = n.children[0]
            child.dist += n.dist
            n.delete(prevent_nondicotomic=False)


def suppress_unary_root(t):
    while len(t.children) == 1:
        child = t.children[0]
        child.detach()
        t = child
        t.dist = 0.0
    return t


def keep_biggest_after_each_cut(t: Tree, cutoff_allbranches: float, excl_fh):
    removed = set()

    while True:
        candidates = [
            n for n in t.traverse()
            if (not n.is_root()) and (n.dist > cutoff_allbranches)
        ]
        if not candidates:
            return t, removed

        n = max(candidates, key=lambda z: z.dist)

        all_leaves = set(t.get_leaf_names())
        sub_leaves = set(n.get_leaf_names())
        rest_leaves = all_leaves - sub_leaves

        if len(sub_leaves) >= len(rest_leaves):
            discarded = rest_leaves
            excl_fh.write(" ".join(sorted(discarded)) + "\n")
            removed |= discarded

            n.detach()
            t = n
            t.dist = 0.0
        else:
            discarded = sub_leaves
            excl_fh.write(" ".join(sorted(discarded)) + "\n")
            removed |= discarded
            n.detach()

        suppress_unary_nodes(t)
        t = suppress_unary_root(t)

        assert all(len(node.children) != 1 for node in t.traverse()), \
            "Unary node detected after suppression"


def prune_long_tips_iteratively(t: Tree, cutoff_tips: float, excl_fh=None):
    """
    Iteratively remove tips whose terminal branch length exceeds cutoff.
    Writes removed tip names per round as one line to excl_fh (if provided).

    Returns:
        (Tree pruned_tree, set[str] removed_tips)
    """
    removed = set()

    while True:
        bad = [leaf for leaf in t.iter_leaves() if leaf.dist > cutoff_tips]
        if not bad:
            return t, removed

        discarded = set()
        for leaf in bad:
            if leaf.name not in (None, ""):
                discarded.add(leaf.name)
            leaf.detach()
            suppress_unary_nodes(t)
            t = suppress_unary_root(t)



        if excl_fh is not None and discarded:
            excl_fh.write(" ".join(sorted(discarded)) + "\n")
            
        removed |= discarded

        assert all(len(node.children) != 1 for node in t.traverse()), \
            "Unary node detected after tip pruning"

def prune_via_purity(t: Tree, maximal_monophyletic_fraction: int, minimal_monophyletic_purity: float, metadata_col: str, metadata_in: str, excl_fh=None ):
    # load metadata once
    max_phylogenetic_size_abs = maximal_monophyletic_fraction * len(t.get_leaf_names()) 
    with open(metadata_in) as fin:
        header = fin.readline().rstrip("\n").split("\t")
        idx = header.index(metadata_col)

        meta = dict(
            (f[0], f[idx])
            for f in (line.rstrip("\n").split("\t") for line in fin if line.strip())
        )

    removed = set()

    while True:
        clades = []
        for node in t.traverse():
            if node.is_leaf() or node.is_root():
                continue
            if (size(node.get_leaf_names()) <= max_phylogenetic_size_abs and (not node.is_root()) and size(node.up.get_leaf_names()) > max_phylogenetic_size_abs):
                clades.append(node)

        if not clades: return t, removed

        discarded = set()

        for clade in clades:
            leaf_names = clade.get_leaf_names()

            counts = {}
            total = 0
            for nm in leaf_names:
                v = meta.get(nm, "")
                if v in ("", None): continue
                counts[v] = counts.get(v, 0) + 1
                total += 1
            if total == 0: continue

            bad_types = {v for v, c in counts.items() if (c / total) < minimal_monophyletic_purity}
            if not bad_types: continue

            for nm in leaf_names:
                if meta.get(nm, "") in bad_types:
                    discarded.add(nm)

        if not discarded:
            return t, removed

        for leaf in list(t.iter_leaves()):
            if leaf.name in discarded:
                leaf.detach()
                suppress_unary_nodes(t)
                t = suppress_unary_root(t)

        if excl_fh is not None:
            excl_fh.write(" ".join(sorted(discarded)) + "\n")

        removed |= discarded

        assert all(len(node.children) != 1 for node in t.traverse()), \
            "Unary node detected after purity pruning"


def prune_metadata(metadata_in: str, metadata_out: str, removed: set):
    with open(metadata_in) as fin:
        header = fin.readline()
        out = [header]

        for line in fin:
            if not line.strip():
                continue
            key = line.split("\t", 1)[0]
            if key not in removed:
                out.append(line)

    with open(metadata_out, "w") as fout:
        fout.writelines(out)


def plot_log_tip_length_distribution(t, out_png, bins=200):
    import matplotlib.pyplot as plt
    import numpy as np

    vals = np.array([leaf.dist for leaf in t.iter_leaves() if leaf.dist > 0])
    if len(vals) == 0:
        return

    lx = np.log10(vals)

    plt.figure(figsize=(6, 4))
    plt.hist(lx, bins=bins)
    plt.yscale("log")
    plt.xlabel("log10(terminal branch length)")
    plt.ylabel("Number of tips (log scale)")
    plt.title("Terminal branch length distribution")
    plt.tight_layout()
    plt.savefig(out_png, dpi=300)
    plt.close()


# -------------------- main --------------------

parser = argparse.ArgumentParser()
parser.add_argument("--tree", required=True)
parser.add_argument("--metadata-in", required=True)
parser.add_argument("--metadata-out", required=True)
parser.add_argument("--out-tree", required=True)
parser.add_argument("--exclude", required=True)
parser.add_argument("--cutoff_allbranches", type=float, required=True)
parser.add_argument("--cutoff_tips", type=float, default=None)
parser.add_argument("--maximal_monophyletic_fraction", type=float, required=False)
parser.add_argument("--minimal_monophyletic_purity", type=float, required=False)
parser.add_argument("--purity_metadata_col", type=str, required=False)
args = parser.parse_args()

t = Tree(args.tree, format=1)
total = len(t.get_leaf_names())

Path(args.exclude).parent.mkdir(parents=True, exist_ok=True)
Path(args.out_tree).parent.mkdir(parents=True, exist_ok=True)
Path(args.metadata_out).parent.mkdir(parents=True, exist_ok=True)
Path("results/tip_length_distr").mkdir(parents=True, exist_ok=True)

# internal-branch pruning
with open(args.exclude, "w") as excl_fh:
    kept_tree, removed_via_all = keep_biggest_after_each_cut(
        t, args.cutoff_allbranches, excl_fh
    )

# terminal-tip pruning (optional)
removed_via_tips = set()
if args.cutoff_tips is not None:
    with open(args.exclude, "a") as excl_fh:
        kept_tree, removed_via_tips = prune_long_tips_iteratively(
            kept_tree, args.cutoff_tips, excl_fh
        )
removed_via_monophyly = set()
if args.maximal_monophyletic_fraction is not None and args.minimal_monophyletic_purity is not None:
    with open(args.exclude, "a") as excl_fh:
        kept_tree, removed_via_monophyly = prune_via_purity(
            kept_tree,
            args.maximal_monophyletic_fraction,
            args.minimal_monophyletic_purity,
            args.purity_metadata_col,
            args.metadata_in,
            excl_fh
        )

# final sanity
assert all(l.name not in (None, "") for l in kept_tree.iter_leaves()), \
    "Unnamed leaf present"

if kept_tree.name in (None, ""):
    kept_tree.name = "ROOT"

kept_tree.write(outfile=args.out_tree, format=1)

removed_all = removed_via_all | removed_via_tips | removed_via_monophyly
prune_metadata(args.metadata_in, args.metadata_out, removed_all)

excluded_via_all = len(removed_via_all)
excluded_via_tips = len(removed_via_tips)
excluded_via_monophyly = len(removed_via_monophyly)

print(
    f"[{Path(args.tree).name}] "
    f"Excluded: \n   {excluded_via_all} via long internal branches > {args.cutoff_allbranches},"
    + (f" \n + {excluded_via_tips} via terminal branches > {args.cutoff_tips}"
       if args.cutoff_tips is not None else "")
    + f"  \n + {excluded_via_monophyly} to enforce purity of {args.minimal_monophyletic_purity} of monophyletic clades with up to {int(args.maximal_monophyletic_fraction * total)} tips"
    + f"  \n out of {total} total tips"
)

plot_log_tip_length_distribution(
    kept_tree,
    f"results/tip_length_distr/{Path(args.tree).stem}.tip_length_loglog.png"
)
