from Bio import AlignIO, SeqIO
import sys

msa = sys.argv[1]
gb = sys.argv[2]

# read alignment and genbank
aln = AlignIO.read(msa, "fasta")
gb_rec = next(SeqIO.parse(gb, "genbank"))

# find genbank sequence in MSA
aln_rec = None

print("gbrec_id", gb_rec.id)
for r in aln:
    if r.id == gb_rec.id.split(".", 1)[0]:
        aln_rec = r
        break

if aln_rec is None:
    raise SystemExit("GenBank sequence not found in MSA")

aln_seq = str(aln_rec.seq)
gb_seq = str(gb_rec.seq).lower()

# check it starts the same (ungapped)
if aln_seq.replace("-", "")[:50] != gb_seq[:50]:
    print(aln_seq.replace("-", "")[:50])
    print(gb_seq[:50])
    raise SystemExit("Start of GenBank sequence does not match MSA")

# map alignment pos -> genbank pos
gb_pos = -1
print("aln_pos\tgb_pos")
for aln_pos, c in enumerate(aln_seq, start=0):
    if c == "-":
        print(f"{aln_pos}\t")
    else:
        gb_pos += 1
        print(f"{aln_pos}\t{gb_pos}")
