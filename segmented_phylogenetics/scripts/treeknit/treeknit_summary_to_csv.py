import os
import csv
import re

def parse_treeknit_summary(summary_path):
    with open(summary_path) as f:
        txt = f.read()

    def grab(pattern, cast=str, default=None):
        m = re.search(pattern, txt)
        return cast(m.group(1)) if m else default

    row = {
        "tree1": grab(r"tree1:\s*(.+)"),
        "tree2": grab(r"tree2:\s*(.+)"),
        "outdir": grab(r"outdir:\s*(.+)"),

        "leaves_tree1": grab(r"tree1:\s*\d+\n\s*tree2:\s*\d+\n\s*shared:\s*(\d+)", int),
        "leaves_tree2": grab(r"tree2:\s*(\d+)", int),
        "shared_leaves": grab(r"shared:\s*(\d+)", int),

        "n_subset_arg": grab(r"n_subset arg:\s*(\d+|None)", lambda x: None if x == "None" else int(x)),
        "K_used": grab(r"K used:\s*(\d+)", int),

        "n_mcc": grab(r"MCCs:\s*(\d+)", int),

        "mcc_min": grab(r"min:\s*(\d+)", int),
        "mcc_median": grab(r"median:\s*(\d+)", int),
        "mcc_max": grab(r"max:\s*(\d+)", int),

        "largest_mcc_frac_K": grab(r"largest MCC / K:\s*([\d.]+)%", float),
        "sum_mcc_frac_K": grab(r"sum\(MCC sizes\) / K:\s*([\d.]+)", float),
        "singleton_frac_K": grab(r"singleton MCCs / K:\s*([\d.]+)%", float),
        "non_singleton_leaf_frac_K": grab(r"non-singleton leaf coverage / K:\s*([\d.]+)%", float),

        "singleton_mccs": grab(r"singletons:\s*(\d+)", int),
        "non_singleton_mccs": grab(r"non-singletons:\s*(\d+)", int),
    }

    return row


def append_summary_csv(summary_path, out_csv):
    row = parse_treeknit_summary(summary_path)

    write_header = not os.path.exists(out_csv)

    with open(out_csv, "a", newline="") as f:
        w = csv.DictWriter(f, fieldnames=row.keys())
        if write_header:
            w.writeheader()
        w.writerow(row)

if __name__ == "__main__":
    import sys
    if len(sys.argv) != 3:
        sys.exit("Usage: python treeknit_summary_to_csv.py <results_summary.txt> <out.csv>")
    append_summary_csv(sys.argv[1], sys.argv[2])
