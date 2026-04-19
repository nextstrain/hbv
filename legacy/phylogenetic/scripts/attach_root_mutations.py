"""
Copied from https://github.com/nextstrain/nextclade_dataset_template/blob/d7e16ab67a784eedb3bfe4327dfe1974e1b9e705/scripts/attach_root_mutations.py
with modifications
Original author: Richard Neher
"""
import argparse
import json
from Bio import SeqIO, Phylo

if __name__=="__main__":
    parser = argparse.ArgumentParser(
        description="add root mutations via the nuc + aa reconstruction produced by augur ancestral",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('--tree', type=str, required=True)
    parser.add_argument('--reference-name', type=str, required=True)
    parser.add_argument('--translations', type=str, metavar='FILE_PATTERN', required=True, help="The translations of the nodes in the tree, from augur ancestral. The reference must be in here! Format: filename with %%GENE")
    parser.add_argument('--genes', type=str, nargs='+', required=True, help="gene names")
    parser.add_argument('--mutations-in',  metavar="NODE_DATA_JSON", type=str, required=True, help="the node-data JSON from augur ancestral")
    parser.add_argument('--mutations-out', metavar="NODE_DATA_JSON",  type=str, required=True, help="the node-data JSON to pass to augur export")
    args = parser.parse_args()
    genes = args.genes if type(args.genes)==list else [args.genes]

    # parse tree, get root name
    T = Phylo.read(args.tree, 'newick')
    root_name = T.root.name

    # Get the AA _reference_ and _root_ sequences. Note that the reference _must_ be in the tree (and not pruned out!)
    references = {}
    root_seqs = {}
    for gene_name in genes:
        for s in SeqIO.parse(args.translations.replace("%GENE", gene_name), 'fasta'):
            if s.id == args.reference_name:
                references[gene_name] = str(s.seq)
            elif s.id == root_name:
                root_seqs[gene_name] = str(s.seq)
        if gene_name not in references:
            raise Exception(f"Reference {args.reference_name} not found in the provided translations for gene {gene_name}")
        if gene_name not in root_seqs:
            raise Exception(f"Root {root_name} not found in the provided translations for gene {gene_name}")

    # Get the same from the node-data file for nuc
    with open(args.mutations_in, 'r') as fh:
        node_data = json.load(fh)
    try:
        references['nuc'] = node_data['nodes'][args.reference_name]['sequence']
        root_seqs['nuc'] = node_data['nodes'][root_name]['sequence']
    except KeyError:
        raise Exception(f"Root {root_name} and/or reference {args.reference_name} nuc sequence not found in {args.mutations_in}")

    # create the arrays where we'll store mutations on the root (in node-data format)
    # import pdb; pdb.set_trace()
    node_data['nodes'][root_name]['muts'] = []
    node_data['nodes'][root_name]['aa_muts'] = {}

    # for each gene, compare the root sequence to reference & store mutations
    for gene_name in genes:
        node_data['nodes'][root_name]['aa_muts'][gene_name] = []
        for pos, (ref, root) in enumerate(zip(references[gene_name], root_seqs[gene_name])):
            if ref!=root:
                node_data['nodes'][root_name]['aa_muts'][gene_name].append(f"{ref}{pos+1}{root}")
    # repeat for nucleotide mutations
    for pos, (ref, root) in enumerate(zip(references['nuc'], root_seqs['nuc'])):
        if ref!=root:
            node_data['nodes'][root_name]['muts'].append(f"{ref}{pos+1}{root}")

    # output files
    with open(args.mutations_out, 'w') as fh:
        json.dump(node_data, fh, indent=2)
