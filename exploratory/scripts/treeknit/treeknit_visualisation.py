"""Render TreeKnit comparison matrices from the aggregated summary CSV.

Write one PNG heatmap per metric into the requested output directory.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import argparse

def short_name(x):
    return Path(x).stem.replace(".tree", "")

def make_sym_matrix(df_sub, value_col):
    labels = sorted(set(df_sub["name1"]) | set(df_sub["name2"]))
    mat = pd.DataFrame(np.nan, index=labels, columns=labels, dtype=float)

    for _, row in df_sub.iterrows():
        a, b = row["name1"], row["name2"]
        v = row[value_col]
        mat.loc[a, b] = v
        mat.loc[b, a] = v

    if value_col == "sum_mcc_frac_K":
        np.fill_diagonal(mat.values, 100.0)
    elif value_col in {"largest_mcc_frac_K", "singleton_frac_K", "non_singleton_leaf_frac_K"}:
        np.fill_diagonal(mat.values, np.nan)
    else:
        np.fill_diagonal(mat.values, 0.0)

    return mat

def plot_matrix(mat, title, out, cmap="viridis", vmin=None, vmax=None, fmt=".2f"):
    fig, ax = plt.subplots(figsize=(1.2 * len(mat.columns), 1.0 * len(mat.index)))
    im = ax.imshow(mat.values, cmap=cmap, vmin=vmin, vmax=vmax)
    ax.set_xticks(range(len(mat.columns)))
    ax.set_yticks(range(len(mat.index)))
    ax.set_xticklabels(mat.columns, rotation=45, ha="right")
    ax.set_yticklabels(mat.index)
    ax.tick_params(top=True, bottom=False, labeltop=True, labelbottom=False)

    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]):
            val = mat.iat[i, j]
            if not np.isnan(val):
                ax.text(j, i, format(val, fmt), ha="center", va="center", fontsize=8)

    ax.set_title(title, pad=20)
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    plt.tight_layout()
    fig.savefig(out, dpi=300)
    plt.close(fig)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("in_csv")
    ap.add_argument("--outdir", default=None)
    args = ap.parse_args()

    in_csv = args.in_csv
    outdir = Path(args.outdir) if args.outdir else Path(in_csv).parent / "plots"
    outdir.mkdir(parents=True, exist_ok=True)
    df = pd.read_csv(in_csv)

    df["name1"] = df["tree1"].map(short_name)
    df["name2"] = df["tree2"].map(short_name)

    metrics = {
        "n_mcc": dict(title="Number of MCCs", cmap="viridis", vmin=0, vmax=None, fmt=".0f"),
        "largest_mcc_frac_K": dict(title="Largest MCC (% of K)", cmap="magma", vmin=0, vmax=None, fmt=".2f"),
        "singleton_frac_K": dict(title="Singleton leaves (% of K)", cmap="magma", vmin=0, vmax=100, fmt=".2f"),
        "non_singleton_leaf_frac_K": dict(title="Leaves in non-singleton MCCs (% of K)", cmap="magma", vmin=0, vmax=100, fmt=".2f"),
        "mcc_max": dict(title="Largest MCC size (n leaves)", cmap="viridis", vmin=0, vmax=None, fmt=".0f"),
    }

    for metric, opts in metrics.items():
        mat = make_sym_matrix(df, metric)
        plot_matrix(
            mat,
            title=f"TreeKnit summary: {opts['title']}",
            out=outdir / f"{metric}.png",
            cmap=opts["cmap"],
            vmin=opts["vmin"],
            vmax=opts["vmax"],
            fmt=opts["fmt"],
        )

if __name__ == "__main__":
    main()
