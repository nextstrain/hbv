"""Render a Robinson-Foulds distance matrix from the aggregated comparison table.

Write one PDF heatmap covering all region pairs in the current exploratory run.
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import sys

tsv = sys.argv[1]
out_pdf = sys.argv[2]

# read
df = pd.read_csv(tsv, sep="\t")

# extract segment names
segs = sorted(s.split(".")[0] for s in set(df["tree1"]) | set(df["tree2"]))
idx = {s: i for i, s in enumerate(segs)}
n = len(segs)

# matrix
M = np.full((n, n), np.nan)

for _, r in df.iterrows():
    i, j = idx[r.tree1.split(".")[0]], idx[r.tree2.split(".")[0]]
    M[i, j] = r.rf_normalized
    M[j, i] = r.rf_normalized  # mirror

fig, ax = plt.subplots(figsize=(0.8*n, 0.8*n))
im = ax.imshow(M, vmin=0, vmax=1)

# ticks
ax.set_xticks(range(n))
ax.set_yticks(range(n))
ax.set_xticklabels(segs, rotation=90)
ax.set_yticklabels(segs)

# overlay numbers (only upper triangle)
for i in range(n):
    for j in range(i+1, n):
        if not np.isnan(M[i, j]):
            ax.text(j, i, f"{M[i,j]:.2f}", ha="center", va="center", fontsize=8)

plt.title("Normalised Robinson-Foulds-Distance: \n (RF/RF_max)",fontsize=7)
plt.colorbar(im, ax=ax, label="RF Normalised")
plt.tight_layout()
plt.savefig(out_pdf)
