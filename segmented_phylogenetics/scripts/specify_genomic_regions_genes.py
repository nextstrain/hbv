from Bio import SeqIO

import sys
from Bio import SeqIO

gb_in = sys.argv[1]
out = sys.argv[2]
# gb_in = "NC_003977.gb"
# out = "genomic_regions_genes.txt"



genes_wanted = {"C", "P", "S", "X"}

rec = SeqIO.read(gb_in, "genbank")
L = len(rec.seq)

# collect all "gene" features by gene name
genes = {}
for feat in rec.features:
    if feat.type != "gene":
        continue
    name = (feat.qualifiers.get("gene") or [None])[0]
    if name not in genes_wanted:
        continue

    loc = feat.location
    parts = list(getattr(loc, "parts", [loc]))  # CompoundLocation -> parts, else single

    # Biopython locations are 0-based, end-exclusive
    starts0 = [int(p.start) for p in parts]
    ends0 = [int(p.end) for p in parts]

    # If gene is split (wrap-around), represent as start>end using first part start and last part end
    starts0_sorted = sorted(starts0)
    ends0_sorted = sorted(ends0)

    # heuristic: if there are 2+ parts and one starts near 0, treat as wrap
    wraps = len(parts) > 1 and min(starts0) == 0

    if not wraps:
        start_1based = min(starts0) + 1
        end_1based_inclusive = max(ends0)
    else:
        # part near end of genome + part near start; represent as start>end
        # choose start from the part with the largest start (tail), end from the part with the largest end (head)
        tail_start0 = max(starts0)
        head_end0 = max(e for s, e in zip(starts0, ends0) if s == 0)
        start_1based = tail_start0 + 1
        end_1based_inclusive = head_end0

    # keep the widest interval if multiple gene annotations exist
    prev = genes.get(name)
    span = (end_1based_inclusive - start_1based + 1) if start_1based <= end_1based_inclusive else (L - start_1based + 1 + end_1based_inclusive)
    if prev is None or span > prev["span"]:
        genes[name] = {"start": start_1based, "end": end_1based_inclusive, "span": span}

with open(out, "w") as fh:
    fh.write("# name\tstart\tend\t(1-based, inclusive; start>end means wrap-around)\n")
    for name in ["C", "P", "S", "X"]:
        g = genes.get(name)
        if g is None:
            raise ValueError(f"Gene {name} not found in {gb_in}")
        fh.write(f"{name}\t{g['start']}\t{g['end']}\n")



