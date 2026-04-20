# treeknit_summarize.py
# Usage: python treeknit_summarize.py tree1.nwk tree2.nwk outdir [n_subset|None]
import json, os, sys

tree1, tree2, outdir = sys.argv[1:4]
n_subset_arg = sys.argv[4] if len(sys.argv) >= 5 else None

t1_name = os.path.basename(tree1)
t2_name = os.path.basename(tree2)

mcc_path     = os.path.join(outdir, "MCCs.json")
param_path   = os.path.join(outdir, "parameters.json")
log_path     = os.path.join(outdir, "log.txt")
summary_path = os.path.join(outdir, "results_summary.txt")

# parsed MCC stats
n_mcc = None
mcc_sizes = []
largest_mcc = None
n_singletons = None
n_non_singletons = None

# denominators
n_leaves_tree1 = None
n_leaves_tree2 = None
n_shared_leaves = None

# subset denominator (what you want to normalize to)
subset_k = None

# normalized metrics (relative to subset_k)
largest_mcc_of_subset = None       # largest_mcc / subset_k
mcc_leaf_coverage = None           # sum(mcc_sizes) / subset_k
singleton_leaf_frac = None         # n_singletons / subset_k
nonsingleton_leaf_coverage = None  # sum(size>1) / subset_k

note = ""

def median_int(xs):
    xs = sorted(xs)
    if not xs:
        return None
    mid = len(xs) // 2
    return xs[mid] if len(xs) % 2 == 1 else xs[mid - 1]

def read_newick_leaves(path):
    s = open(path).read()
    labels = set()
    token = []
    in_quote = False
    quote_char = None

    def flush():
        nonlocal token
        if not token:
            return
        lab = "".join(token).strip()
        token = []
        if not lab:
            return
        try:
            float(lab)
            return
        except Exception:
            pass
        labels.add(lab)

    for ch in s:
        if in_quote:
            if ch == quote_char:
                in_quote = False
            else:
                token.append(ch)
            continue

        if ch in ("'", '"'):
            in_quote = True
            quote_char = ch
            continue

        if ch in "(),:;":
            flush()
        else:
            token.append(ch)
    flush()

    # filter common internal node labels in your trees
    labels = {x for x in labels if not x.startswith("NODE_")}
    return labels

def _collect_mcc_sizes_from_json(data):
    sizes = []
    mcc_lists = []

    if isinstance(data, dict) and "MCC_dict" in data:
        mcc_dict = data.get("MCC_dict")
        if isinstance(mcc_dict, dict):
            for entry in mcc_dict.values():
                if not isinstance(entry, dict):
                    continue
                mccs = entry.get("mccs", [])
                if isinstance(mccs, list):
                    for m in mccs:
                        if isinstance(m, list):
                            sizes.append(len(m))
                            mcc_lists.append(m)
                        elif isinstance(m, dict):
                            if isinstance(m.get("leaves"), list):
                                sizes.append(len(m["leaves"]))
                                mcc_lists.append(m["leaves"])
                            elif isinstance(m.get("tips"), list):
                                sizes.append(len(m["tips"]))
                                mcc_lists.append(m["tips"])
        return sizes, mcc_lists

    if isinstance(data, dict) and "MCCs" in data:
        mccs = data.get("MCCs")
        if isinstance(mccs, list):
            for m in mccs:
                if isinstance(m, dict):
                    if isinstance(m.get("leaves"), list):
                        sizes.append(len(m["leaves"]))
                        mcc_lists.append(m["leaves"])
                    elif isinstance(m.get("tips"), list):
                        sizes.append(len(m["tips"]))
                        mcc_lists.append(m["tips"])
                elif isinstance(m, list):
                    sizes.append(len(m))
                    mcc_lists.append(m)
        return sizes, mcc_lists

    return sizes, mcc_lists

def parse_subset_k(n_subset_arg, n_shared_leaves):
    """
    If n_subset was used for pruning, the effective denominator is:
      K = min(n_subset, n_shared_leaves)   (because you can't sample more than shared)
    If n_subset is None/"None"/"" -> fall back to n_shared_leaves.
    """
    if n_subset_arg is None:
        return n_shared_leaves
    s = str(n_subset_arg).strip()
    if s == "" or s.lower() == "none":
        return n_shared_leaves
    try:
        k = int(s)
        if k <= 0:
            return n_shared_leaves
        if n_shared_leaves is None:
            return k
        return min(k, n_shared_leaves)
    except Exception:
        return n_shared_leaves

# --- parse MCCs.json ---
mcc_lists = []
if os.path.isfile(mcc_path):
    parsed_ok = False
    try:
        with open(mcc_path) as f:
            data = json.load(f)

        mcc_sizes, mcc_lists = _collect_mcc_sizes_from_json(data)

        if mcc_sizes:
            parsed_ok = True
            n_mcc = len(mcc_sizes)
            largest_mcc = max(mcc_sizes)
        else:
            keys = list(data.keys()) if isinstance(data, dict) else None
            note += "Parsed MCCs.json but found no MCCs in a recognized schema."
            if keys is not None:
                note += f" Top-level keys: {keys}\n"
            else:
                note += "\n"

    except Exception as e:
        note += f"Failed to parse MCCs.json as JSON ({type(e).__name__}: {e}).\n"

    if not parsed_ok and n_mcc is None:
        try:
            txt = open(mcc_path).read()
            n_mcc = txt.count('"leaves"') or txt.count('"mccs"')
            note += 'Fallback heuristic: counted occurrences of \'"leaves"\' or \'"mccs"\' in MCCs.json.\n'
        except Exception:
            note += "Failed to read MCCs.json.\n"
else:
    note += "Missing MCCs.json; TreeKnit may have failed or wrote elsewhere.\n"

# --- compute denominators from trees ---
try:
    leaves1 = read_newick_leaves(tree1)
    leaves2 = read_newick_leaves(tree2)
    n_leaves_tree1 = len(leaves1)
    n_leaves_tree2 = len(leaves2)
    shared = leaves1 & leaves2
    n_shared_leaves = len(shared)
except Exception as e:
    note += f"Failed to parse Newick leaves ({type(e).__name__}: {e}).\n"

subset_k = parse_subset_k(n_subset_arg, n_shared_leaves)

# --- derived metrics relative to subset_k ---
if mcc_sizes:
    n_singletons = sum(1 for s in mcc_sizes if s == 1)
    n_non_singletons = sum(1 for s in mcc_sizes if s > 1)
    sum_sizes = sum(mcc_sizes)
    sum_nonsingletons = sum(s for s in mcc_sizes if s > 1)

    if subset_k and subset_k > 0:
        largest_mcc_of_subset = (largest_mcc / subset_k) if largest_mcc is not None else None
        mcc_leaf_coverage = (sum_sizes / subset_k)
        singleton_leaf_frac = (n_singletons / subset_k)
        nonsingleton_leaf_coverage = (sum_nonsingletons / subset_k)
    else:
        note += "Could not determine subset denominator K; set n_subset (4th CLI arg) or ensure shared leaves > 0.\n"

# ---- print to stdout ----
print("TreeKnit results")
print(f"  tree1: {t1_name}")
print(f"  tree2: {t2_name}")
print(f"  outdir: {outdir}")
print(f"  MCCs.json: {mcc_path}")

if n_leaves_tree1 is not None and n_leaves_tree2 is not None:
    print(f"  leaves: tree1={n_leaves_tree1} tree2={n_leaves_tree2} shared={n_shared_leaves}")

print(f"  subset K (denominator): {subset_k}  (from n_subset={n_subset_arg})")

if n_mcc is not None:
    print(f"  MCCs: {n_mcc}")

if mcc_sizes:
    med = median_int(mcc_sizes)
    print(f"  MCC sizes: min={min(mcc_sizes)} median={med} max={max(mcc_sizes)}")
    print(f"  MCC composition: singletons={n_singletons} non-singletons={n_non_singletons}")

    if largest_mcc_of_subset is not None:
        print(f"  largest MCC / K: {largest_mcc_of_subset:.2%} ({largest_mcc}/{subset_k})")
    if mcc_leaf_coverage is not None:
        print(f"  sum(MCC sizes) / K: {mcc_leaf_coverage:.2f}")
    if singleton_leaf_frac is not None:
        print(f"  singleton MCCs / K: {singleton_leaf_frac:.2%} ({n_singletons}/{subset_k})")
    if nonsingleton_leaf_coverage is not None:
        print(f"  non-singleton leaf coverage / K: {nonsingleton_leaf_coverage:.2%}")
if note.strip():
    print("  note:", note.strip())

# ---- write human-readable summary ----
with open(summary_path, "w") as io:
    io.write("TreeKnit results summary\n")
    io.write("=======================\n")
    io.write(f"tree1:\t{tree1}\n")
    io.write(f"tree2:\t{tree2}\n")
    io.write(f"outdir:\t{outdir}\n\n")

    io.write("Key outputs:\n")
    io.write(f"  MCCs.json:\t{mcc_path}\n")
    io.write(f"  parameters.json:\t{param_path}\n")
    io.write(f"  log.txt:\t{log_path}\n\n")

    if n_leaves_tree1 is not None and n_leaves_tree2 is not None:
        io.write("Leaves:\n")
        io.write(f"  tree1:\t{n_leaves_tree1}\n")
        io.write(f"  tree2:\t{n_leaves_tree2}\n")
        io.write(f"  shared:\t{n_shared_leaves}\n\n")

    io.write("Subset normalization:\n")
    io.write(f"  n_subset arg:\t{n_subset_arg}\n")
    io.write(f"  K used:\t{subset_k}\n")
    io.write("  (K = min(n_subset, #shared leaves); if n_subset=None -> K=#shared leaves)\n\n")

    if n_mcc is not None:
        io.write(f"MCCs:\t{n_mcc}\n")

    if mcc_sizes:
        mcc_sizes_sorted = sorted(mcc_sizes)
        med = median_int(mcc_sizes_sorted)
        io.write("MCC size stats (by #leaves listed per MCC):\n")
        io.write(f"  min:\t{mcc_sizes_sorted[0]}\n")
        io.write(f"  median:\t{med}\n")
        io.write(f"  max:\t{mcc_sizes_sorted[-1]}\n\n")

        io.write("Normalized metrics (relative to K):\n")
        if largest_mcc_of_subset is not None:
            io.write(f"  largest MCC / K:\t{largest_mcc_of_subset:.2%}\t({largest_mcc}/{subset_k})\n")
        if mcc_leaf_coverage is not None:
            io.write(f"  sum(MCC sizes) / K:\t{mcc_leaf_coverage:.2f}\n")
        if singleton_leaf_frac is not None:
            io.write(f"  singleton MCCs / K:\t{singleton_leaf_frac:.2%}\t({n_singletons}/{subset_k})\n")
        if nonsingleton_leaf_coverage is not None:
            io.write(f"  non-singleton leaf coverage / K:\t{nonsingleton_leaf_coverage:.2%}\n")

        io.write("\nRaw MCC composition:\n")
        io.write(f"  singletons:\t{n_singletons}\n")
        io.write(f"  non-singletons:\t{n_non_singletons}\n")

    if note.strip():
        io.write("\nNotes:\n")
        io.write(note if note.endswith("\n") else note + "\n")
