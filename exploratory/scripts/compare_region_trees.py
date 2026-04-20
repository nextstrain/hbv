from ete3 import Tree
import sys, os

t1 = Tree(sys.argv[1])
t2 = Tree(sys.argv[2])
out_path = sys.argv[3]

# restrict to shared leaves
common = set(t1.get_leaf_names()) & set(t2.get_leaf_names())
t1.prune(common, preserve_branch_length=True)
t2.prune(common, preserve_branch_length=True)

res = t1.robinson_foulds(t2, unrooted_trees=True)
rf = res[0]
max_rf = res[1]
common_leaves = res[2]

rf_norm = rf / max_rf if max_rf else 0.0
print("RF:", rf, "max:", max_rf, "normalized:", rf_norm, "n_leaves:", len(common))

header = (
    "tree1\t"
    "tree2\t"
    "n_shared_leaves\t"
    "rf\t"
    "max_rf\t"
    "rf_normalized\n"
)

row = (
    f"{os.path.basename(sys.argv[1])}\t"
    f"{os.path.basename(sys.argv[2])}\t"
    f"{len(common)}\t"
    f"{rf}\t"
    f"{max_rf}\t"
    f"{rf_norm}\n"
)

write_header = (not os.path.exists(out_path)) or os.path.getsize(out_path) == 0

with open(out_path, "a") as f:
    if write_header:
        f.write(header)
    f.write(row)
