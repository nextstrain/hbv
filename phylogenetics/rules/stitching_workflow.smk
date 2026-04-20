#______________________________________________________________________________________________________________________________________________________________________________________________
# Rules for stitching together trees from different genotypes into a single tree 
#______________________________________________________________________________________________________________________________________________________________________________________________

rule collect_trees_for_stitching:
    input:
        trees=lambda wc: expand(
            RESULTS + "/stitched/{g}/{gene}_masked/{gene}_masked_refined.tree.nwk",
            g=ALL_GTS, gene=[wc.gene]
        ),
        alns=lambda wc: expand(
            RESULTS + "/stitched/{g}/{gene}_masked/{gene}_masked_aln.fasta",
            g=ALL_GTS, gene=[wc.gene]
        ),
    output:
        trees_dir=directory(RESULTS + "/stitched/{gene}_global/trees"),
        aln_dir=directory(RESULTS + "/stitched/{gene}_global/aln"),
        touch=RESULTS + "/stitched/{gene}_global/collect.done",
    params:
        g_list=" ".join(ALL_GTS)
    shell:
        r"""
        mkdir -p {output.trees_dir} {output.aln_dir}

        for g in {params.g_list}; do
            cp {RESULTS}/stitched/$g/{wildcards.gene}_masked/{wildcards.gene}_masked_refined.tree.nwk "{output.trees_dir}/$g.nwk"
            cp {RESULTS}/stitched/$g/{wildcards.gene}_masked/{wildcards.gene}_masked_aln.fasta "{output.aln_dir}/$g.fasta"
        done

        touch {output.touch}
        """
rule group_trees:
    input:
        done=RESULTS + "/stitched/{gene}_global/collect.done"
    output:
        stitched_tree=RESULTS + "/stitched/{gene}_global/tree_raw.nwk"
    params:
        trees=RESULTS + "/stitched/{gene}_global/trees",
        aln=RESULTS + "/stitched/{gene}_global/aln",
        outdir=directory(RESULTS + "/stitched/{gene}_global"),
        ref_id = config["reference"]["id"]
    shell:
        r"""
        python3 scripts/make_global_tree.py \
          --tree-dir {params.trees} \
          --alignment-dir {params.aln} \
          --output {output.stitched_tree} \
          --outdir {params.outdir} \
          --ref-id {params.ref_id}
        """


rule refine_stitched:
    input:
        tree = STITCHED_DIR + "/tree_raw.nwk",
        metadata="data/filtered/metadata.len_filtered.tsv",
    output:
        tree = STITCHED_DIR + "/{gene}_tree.nwk",
        node_data = STITCHED_DIR + "/node_data/{gene}_branch_lengths.json",
    shell:
        r"""
        augur refine --tree {input.tree} --metadata {input.metadata} \
                     --keep-root \
                     --output-tree {output.tree} --output-node-data {output.node_data}
        """




rule ancestral_stitched:
    input:
        tree = STITCHED_DIR + "/{gene}_tree.nwk",
        aln=expand(RESULTS + "/stitched/{key}/filtered.fasta", key=ALL_GTS),  # Using non-masked alignments for ancestral reconstruction
        annotation= config["reference"]["gff"],
        translations=expand("../ingest/data/nextclade/cds_{g}.fasta", g=config["genes"]), 
        root = "../nextclade_datasets/references/NC_003977/versions/2023-08-22/reference.fasta",   # Mutations are relative to the reference sequence

    output:
        node_data = STITCHED_DIR + "/node_data/{gene}_muts.json",
    params:
        genes=" ".join(ANCESTRAL_GENES),
        translation_pattern="../ingest/data/nextclade/cds_%GENE.fasta",
        outdir=RESULTS + "/stitched/{gene}_global",
        ref_id= config["reference"]["id"],
 
    shell:
        r"""
        set -euo pipefail

        outdir="{params.outdir}"

        combined="$outdir/combined_nonmasked_aln.fasta"

        awk -v ref="{params.ref_id}" '
        /^>/ {{
            is_ref = ($0 == ">" ref) 
            if (is_ref && seen_ref) skip=1
            else {{
                skip=0
                if (is_ref) seen_ref=1
            }}
        }}
        !skip {{ print }}
        ' {input.aln} > "$combined"

        augur ancestral --tree {input.tree} --alignment "$combined" \
                        --translations {params.translation_pattern} --genes {params.genes} \
                        --output-node-data {output.node_data} \
                        --annotation {input.annotation}  \
                        --root-sequence {input.root} 
        """
