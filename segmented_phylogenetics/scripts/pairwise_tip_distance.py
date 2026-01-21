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



def calculate_distances_both(trees, pairs):
    pat_all, topo_all = [], []
    for tree in trees:
        lookup = {t.name: t for t in tree.get_terminals() if t.name is not None}
        pat, topo = [], []
        for n1, n2 in pairs:
            a, b = lookup[n1], lookup[n2]

            # patristic
            pat.append(tree.distance(a, b))

            # topological (edges)
            m = tree.common_ancestor(a, b)
            topo.append(len(tree.get_path(a)) + len(tree.get_path(b)) - 2 * len(tree.get_path(m)))

        pat_all.append(pat)
        topo_all.append(topo)
    return pat_all, topo_all


import matplotlib.pyplot as plt



def calculate_correlation(values_per_tree):
    import numpy as np
    from scipy.stats import pearsonr

    n = len(values_per_tree)
    corr_matrix = np.zeros((n, n))

    for i in range(n):
        for j in range(i, n):
            x, y = values_per_tree[i], values_per_tree[j]
            corr = 1.0 if i == j else pearsonr(x, y)[0]
            if corr != corr:  # NaN check
                corr = 0.0

            corr_matrix[i, j] = corr
            corr_matrix[j, i] = corr

    return corr_matrix

def save_corr_heatmap(corr_matrix, labels, n_pairs, out_png, metric):
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
    metric="Patristic" if metric=="patristic" else "Topological"
    ax.set_title(f"Pairwise tip-distance correlation (n = {n_pairs} pairs), metric={metric}", pad=20)

    cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label("Pearson correlation")

    plt.tight_layout()
    plt.savefig(out_png, dpi=300)
    plt.close()


def plot_correlation_grid(values, labels, metric, out_png):
    import matplotlib.pyplot as plt

    n = len(values)
    fig, axes = plt.subplots(n, n, figsize=(2*n, 2*n))

    for i in range(n):
        for j in range(n):
            ax = axes[i, j]
            ax.scatter(values[i], values[j], s=2, alpha=0.3)
            ax.set_xticks([])
            ax.set_yticks([])
            if i == n - 1:
                ax.set_xlabel(labels[j], rotation=90)
            if j == 0:
                ax.set_ylabel(labels[i])

    fig.suptitle(f"{metric} distance scatter matrix")
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

    # compute metrics (per tree, per pair)
    pat_dists, topo_dists = calculate_distances_both(trees, pairs)

    labels = [os.path.basename(p).replace(".tree.nwk", "") for p in args.trees]

    for metric, values_per_tree in (("patristic", pat_dists), ("topo", topo_dists)):
        corr_matrix = calculate_correlation(values_per_tree)

        tsv_out = f"results/pairwise_tip_distance_correlation_matrix.{metric}.tsv"
        png_out = f"results/pairwise_tip_distance_correlation_matrix.{metric}.png"

        with open(tsv_out, "w") as f:
            header = "\t" + "\t".join(os.path.basename(p) for p in args.trees)
            print(header, file=f)
            for i, row in enumerate(corr_matrix):
                line = os.path.basename(args.trees[i]) + "\t" + "\t".join(f"{v:.4f}" for v in row)
                print(line, file=f)

        save_corr_heatmap(corr_matrix, labels, args.sample_pairs, png_out, metric)

        # ---- ADD THIS ----
        grid_out = f"results/pairwise_tip_distance_scatter_grid.{metric}.png"
        plot_correlation_grid(
            values_per_tree,
            labels,
            metric=metric,
            out_png=grid_out,
        )



if __name__ == "__main__":
    main()

