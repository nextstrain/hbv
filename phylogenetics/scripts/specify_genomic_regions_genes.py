"""Write phylogenetics gene intervals from the reference GenBank annotation.

Handle wrapped coordinates on the circular HBV genome for the C, P, S, and X genes.
"""

from Bio import SeqIO
import sys

gb_in, out = sys.argv[1], sys.argv[2]
genes_wanted = ["C", "P", "S", "X"]

rec = SeqIO.read(gb_in, "genbank")

genes = {}
for f in rec.features:
    if f.type != "gene":
        continue
    name = (f.qualifiers.get("gene") or [None])[0]

    loc = f.location
    parts = loc.parts if hasattr(loc, "parts") else [loc]

    if len(parts) == 1:
        start = int(parts[0].start) + 1
        end   = int(parts[0].end)
    else:
        # exactly 2 parts, one starts at 0
        p0, p1 = parts
        tail = p0 if int(p0.start) > int(p1.start) else p1
        head = p1 if tail is p0 else p0
        start = int(tail.start) + 1
        end   = int(head.end)

    genes[name] = (start, end)

with open(out, "w") as fh:
    fh.write("# name\tstart\tend\t(1-based, inclusive; start>end means wrap-around)\n")
    for name in genes_wanted:
        if name not in genes:
            raise ValueError(f"Gene {name} not found in {gb_in}")
        s, e = genes[name]
        fh.write(f"{name}\t{s}\t{e}\n")
