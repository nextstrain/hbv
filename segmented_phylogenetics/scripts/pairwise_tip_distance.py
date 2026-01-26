import argparse
import os
import random
import csv
import matplotlib.pyplot as plt
import numpy as np
from Bio import Phylo
from math import ceil
from scipy.stats import pearsonr
from mpl_toolkits.axes_grid1.inset_locator import inset_axes


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





def calculate_correlation(values_per_tree, pairs, out_csv, labels, metric, top_frac=0.05):

    n = len(values_per_tree)
    if n == 0: raise ValueError("values_per_tree is empty")

    for k, v in enumerate(values_per_tree): 
        if len(v) != len(pairs): raise ValueError(f"Tree {k}: expected {len(pairs)} values (len(pairs)), got {len(v)}")

    # --- correlation matrix ---
    corr_matrix = np.zeros((n, n), dtype=float)
    for i in range(n):
        for j in range(i, n):
            if i == j: corr = 1.0
            else:
                r = pearsonr(values_per_tree[i], values_per_tree[j])[0]
                corr = 0.0 if (r != r) else float(r)  # NaN -> 0
            corr_matrix[i, j] = corr
            corr_matrix[j, i] = corr

    if metric=="topo": 
        return corr_matrix

    # --- outliers for EACH tree pair (i<j), one CSV per pair, only for patristic distance ---
    for i in range(n):
        x = np.asarray(values_per_tree[i], dtype=float)
        for j in range(i+1, n):
            y = np.asarray(values_per_tree[j], dtype=float)

            # best-fit line y = a*x + b
            a, b = np.polyfit(x, y, 1)

            # orthogonal distances
            denom = (a * a + 1.0) ** 0.5
            orth = np.abs(a * x - y + b) / denom

            k = max(1, int(ceil(top_frac * len(pairs))))
            idx = np.argsort(-orth)[:k]

            out_csv_ij = os.path.join(
                out_csv,
                f"{labels[i]}_vs_{labels[j]}.csv"
            )

            rows = []
            for t in idx:
                n1, n2 = pairs[t]
                rows.append((n1, n2, round(x[t], 2), round(y[t], 2), round((orth[t]),2)))

            rows.sort(key=lambda r: r[-1], reverse=True)

            with open(out_csv_ij, "w", newline="") as f:
                w = csv.writer(f)
                w.writerow(["tip1", "tip2", f"dist_{labels[i]}", f"dist_{labels[j]}", "orthogonal_distance"])
                w.writerows(rows)

    return corr_matrix




def save_corr_heatmap(corr_matrix, labels, n_pairs, out_png, metric):

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


def plot_correlation_grid(values, labels, metric, out_pdf):
    n = len(values)

    out_base = os.path.splitext(out_pdf)[0]
    out_dir_scatter = out_base + "_scatter_individual"
    out_dir_hex = out_base + "_hexbin_individual"
    os.makedirs(out_dir_scatter, exist_ok=True)
    os.makedirs(out_dir_hex, exist_ok=True)

    # --- scatter grid ---
    fig, axes = plt.subplots(n, n, figsize=(2 * n, 2 * n), constrained_layout=True)

    # --- hexbin grid (reserve right margin for shared colorbar) ---
    fig_hex, axes_hex = plt.subplots(n, n, figsize=(2 * n, 2 * n))
    fig_hex.subplots_adjust(left=0.06, right=0.86, bottom=0.06, top=0.92, wspace=0.35, hspace=0.35)

    hbs = []

    for i in range(n):
        for j in range(n):
            x = values[i]
            y = values[j]

            # ------------------ scatter grid ------------------
            ax = axes[i, j]
            ax.scatter(x, y, s=2, alpha=0.3)
            ax.tick_params(labelsize=8)
            if i == n - 1:
                ax.set_xlabel(labels[j], rotation=90)
            if j == 0:
                ax.set_ylabel(labels[i])

            # ---- individual scatter ----
            fig_ij, ax_ij = plt.subplots(figsize=(4, 4))
            ax_ij.scatter(x, y, s=2, alpha=0.3)
            ax_ij.set_xlabel(labels[j])
            ax_ij.set_ylabel(labels[i])
            ax_ij.tick_params(labelsize=10)
            fig_ij.savefig(
                os.path.join(out_dir_scatter, f"{metric}_{labels[i]}_vs_{labels[j]}.pdf"),
                bbox_inches="tight"
            )
            plt.close(fig_ij)

            # ------------------ hexbin grid ------------------
            axh = axes_hex[i, j]
            hb = axh.hexbin(
                x, y,
                gridsize=40,
                bins="log",
                mincnt=1
            )
            hbs.append(hb)

            axh.tick_params(labelsize=8)
            if i == n - 1:
                axh.set_xlabel(labels[j], rotation=90)
            if j == 0:
                axh.set_ylabel(labels[i])

            # ---- individual hexbin plot with its own colorbar ----
            fig_h, ax_h = plt.subplots(figsize=(4, 4))
            hb2 = ax_h.hexbin(
                x, y,
                gridsize=60,
                bins="log",
                mincnt=1
            )
            ax_h.set_xlabel(labels[j])
            ax_h.set_ylabel(labels[i])
            ax_h.tick_params(labelsize=10)

            cb2 = fig_h.colorbar(hb2, ax=ax_h)
            cb2.set_label("count (log)")
            cb2.ax.tick_params(labelsize=9)

            fig_h.savefig(
                os.path.join(out_dir_hex, f"{metric}_{labels[i]}_vs_{labels[j]}.pdf"),
                bbox_inches="tight"
            )
            plt.close(fig_h)

    # --- shared color scale for hex grid ---
    vmax = max((h.get_array().max() for h in hbs if h.get_array().size), default=1.0)
    for h in hbs:
        if h.get_array().size:
            h.set_clim(1, vmax)

    # --- one shared colorbar (outside grid) ---
    cax = fig_hex.add_axes([0.88, 0.15, 0.02, 0.70])
    cbar = fig_hex.colorbar(hbs[0], cax=cax)
    cbar.set_label("count (log)")

    # --- save scatter grid ---
    fig.suptitle(f"{metric} distance scatter matrix")
    fig.savefig(out_base + "_scatter.pdf", bbox_inches="tight")
    plt.close(fig)

    # --- save hexbin grid (NO tight_layout) ---
    fig_hex.suptitle(f"{metric} distance hexbin matrix (log density)")
    fig_hex.savefig(out_base + "_hexbin.pdf", bbox_inches="tight")
    plt.close(fig_hex)


def main():
    ap = argparse.ArgumentParser(
        description="Keep only shared tips across Newick trees, sample random tip pairs, correlate distances between trees."
    )
    ap.add_argument("trees", nargs="+", help="Input .nwk tree files")
    ap.add_argument("sample_pairs", type=int, help="Number of random tip pairs to sample")
    ap.add_argument("outdir", default="results/correlation_analysis", help="output directory")
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
    outdir = args.outdir
    os.makedirs(os.path.join(outdir, "corr_matrix"), exist_ok=True)
    os.makedirs(os.path.join(outdir, "corr_grid"), exist_ok=True)
    os.makedirs(os.path.join(outdir, "outliers"), exist_ok=True)


    for metric, values_per_tree in (("patristic", pat_dists), ("topo", topo_dists)):

        corr_matrix = calculate_correlation(values_per_tree,pairs, out_csv=os.path.join(outdir, "outliers"), top_frac=0.05, labels=labels, metric=metric)
        
        tsv_out = os.path.join(outdir, "corr_matrix", f"pw_tip_distance_correlation_matrix.{metric}.tsv")
        png_out = os.path.join(outdir, "corr_matrix", f"pw_tip_distance_correlation_matrix.{metric}.png")


        with open(tsv_out, "w") as f:
            header = "\t" + "\t".join(os.path.basename(p) for p in args.trees)
            print(header, file=f)
            for i, row in enumerate(corr_matrix):
                line = os.path.basename(args.trees[i]) + "\t" + "\t".join(f"{v:.4f}" for v in row)
                print(line, file=f)


        save_corr_heatmap(corr_matrix, labels, args.sample_pairs, png_out, metric)

        grid_out = grid_out = os.path.join(outdir, "corr_grid", f"grid_{metric}.pdf")

        plot_correlation_grid(
            values_per_tree,
            labels,
            metric=metric,
            out_pdf=grid_out,
        )



if __name__ == "__main__":
    main()

