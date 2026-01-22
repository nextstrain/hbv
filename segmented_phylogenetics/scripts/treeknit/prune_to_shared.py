# prune_to_shared.py
# Usage:
#   python prune_to_shared.py tree1.nwk tree2.nwk out1.nwk out2.nwk [n] [seed]
#
# If n is provided: keep a random subset of size n from the shared tips.

from Bio import Phylo
import sys, random

t1_path, t2_path, out1, out2 = sys.argv[1:5]
n = int(sys.argv[5]) if len(sys.argv) >= 6 and sys.argv[5] != "None" else None
seed = int(sys.argv[6]) if len(sys.argv) >= 7 and sys.argv[6] != "None" else None

t1 = Phylo.read(t1_path, "newick")
t2 = Phylo.read(t2_path, "newick")

def norm(x):
    if x is None:
        return None
    x = x.strip()
    if len(x) >= 2 and x[0] == x[-1] == '"':
        x = x[1:-1]
    return x

for term in t1.get_terminals():
    term.name = norm(term.name)
for term in t2.get_terminals():
    term.name = norm(term.name)

tips1 = {t.name for t in t1.get_terminals()}
tips2 = {t.name for t in t2.get_terminals()}
shared = tips1 & tips2
if not shared:
    raise SystemExit("No shared leaves between the two trees.")

# optional downsample
# optional downsample
if n is not None:
    if n <= 0:
        raise SystemExit("n must be > 0")
    if n > len(shared):
        print(f"Requested n={n} but only #shared={len(shared)}; keeping all shared leaves (no subsetting).")
        # keep shared as-is
    else:
        rng = random.Random(seed)
        shared = set(rng.sample(sorted(shared), n))

for tip in list(tips1 - shared):
    t1.prune(target=tip)
for tip in list(tips2 - shared):
    t2.prune(target=tip)

tips1p = {t.name for t in t1.get_terminals()}
tips2p = {t.name for t in t2.get_terminals()}
if tips1p != tips2p:
    only1 = sorted(tips1p - tips2p)[:10]
    only2 = sorted(tips2p - tips1p)[:10]
    raise SystemExit(f"Pruning failed.\nOnly in tree1: {only1}\nOnly in tree2: {only2}")

Phylo.write(t1, out1, "newick")
Phylo.write(t2, out2, "newick")

print(f"Shared leaves kept: {len(shared)}")
