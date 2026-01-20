import argparse
import os
import random
from Bio import Phylo


def read_tree(path: str):
    return Phylo.read(path, "newick")


def tip_names(tree):
    return {t.name for t in tree.get_terminals() if t.name is not None}


def prune_to_set(tree, keep):
    # prune all leaves not in keep
    to_prune = [t for t in tree.get_terminals() if t.name not in keep]
    for t in to_prune:
        tree.prune(t)


def select_random_pairs(names, n_pairs, rng):
    names = sorted(names)
    if len(names) < 2:
        raise ValueError("Need at least 2 tips to form pairs")
    return [tuple(rng.sample(names, 2)) for _ in range(n_pairs)]


def calculate_distance(trees, pairs):
    tree_distances = []
    for tree in trees:
        lookup = {t.name: t for t in tree.get_terminals() if t.name is not None}
        distances = []
        for n1, n2 in pairs:
            dist = tree.distance(lookup[n1], lookup[n2])
            distances.append((n1, n2, dist))
        tree_distances.append(distances)
    return tree_distances

def calculate_correlation(distances):
    # distances is a list of lists of (name1, name2, dist) for each tree
    import numpy as np
    from scipy.stats import pearsonr

    n_trees = len(distances)
    corr_matrix = np.zeros((n_trees, n_trees))

    for i in range(n_trees):
        dists_i = [d[2] for d in distances[i]]
        for j in range(i, n_trees):
            dists_j = [d[2] for d in distances[j]]
            corr, _ = pearsonr(dists_i, dists_j)
            corr_matrix[i, j] = corr
            corr_matrix[j, i] = corr

    return corr_matrix

def save_corr_heatmap(corr_matrix, labels, n_pairs, out_png):
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(1.2 * len(labels), 1.0 * len(labels)))

    im = ax.imshow(corr_matrix, cmap="RdBu_r", vmin=-1, vmax=1)
    ax.set_aspect("equal")

    ax.set_xticks(range(len(labels)))
    ax.set_yticks(range(len(labels)))
    ax.set_xticklabels(labels, rotation=45, ha="right")
    ax.set_yticklabels(labels)

    ax.tick_params(top=True, bottom=False,
                   labeltop=True, labelbottom=False)

    for i in range(len(labels)):
        for j in range(len(labels)):
            ax.text(j, i, f"{corr_matrix[i, j]:.2f}",
                    ha="center", va="center", fontsize=8)

    ax.set_title(f"Pairwise tip-distance correlation (n = {n_pairs} pairs)", pad=20)

    cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label("Pearson correlation")

    plt.tight_layout()
    plt.savefig(out_png, dpi=300)
    plt.close()


def main():
    ap = argparse.ArgumentParser(
        description="Keep only shared tips across Newick trees, sample random tip pairs, correlate distances between trees."
    )
    ap.add_argument("trees", nargs="+", help="Input .nwk tree files")
    ap.add_argument("sample_pairs", type=int, help="Number of random tip pairs to sample")
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()

    trees = [read_tree(p) for p in args.trees]
    shared = set.intersection(*(tip_names(t) for t in trees))
    if not shared:
        raise SystemExit("No shared tips across all trees.")

    if args.sample_pairs <= 0:
        raise SystemExit("sample_pairs must be a positive integer.")

    # prune each tree to the shared tip set (in place)
    for path, tree in zip(args.trees, trees):
        prune_to_set(tree, shared)
        out = f"{os.path.splitext(path)[0]}.shared.nwk"
        Phylo.write(tree, out, "newick")

    rng = random.Random(args.seed)
    pairs = select_random_pairs(shared, args.sample_pairs, rng)

    tree_distances = calculate_distance(trees, pairs)
    corr_matrix = calculate_correlation(tree_distances)

    with open("results/pairwise_tip_distance_correlation_matrix.tsv", "w") as f:
        header = "\t" + "\t".join(os.path.basename(p) for p in args.trees)
        print(header, file=f)
        for i, row in enumerate(corr_matrix):
            line = os.path.basename(args.trees[i]) + "\t" + "\t".join(f"{v:.4f}" for v in row)
            print(line, file=f)

    labels = [os.path.basename(p).replace(".tree.nwk", "") for p in args.trees]
    save_corr_heatmap(
    corr_matrix,
    labels,
    args.sample_pairs,
    "results/pairwise_tip_distance_correlation_matrix.png",
)







if __name__ == "__main__":
    main()

