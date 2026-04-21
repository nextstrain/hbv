"""
Rules for stitching genotype-specific trees into a single global tree.
"""

rule collect_trees_for_stitching:
    """Collect genotype-specific trees and alignments into one staging directory per gene."""
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
    """Build an initial stitched tree from the collected genotype-specific trees."""
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
    """Refine the stitched tree against the full length-filtered metadata."""
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

rule combine_stitched_alignments:
    """Combine genotype-specific alignments while keeping only one copy of the reference."""
    input:
        done=RESULTS + "/stitched/{gene}_global/collect.done",
        alns=expand(RESULTS + "/stitched/{key}/filtered.fasta", key=ALL_GTS),
    output:
        alignment=temp(STITCHED_DIR + "/combined_nonmasked_aln.fasta"),
    params:
        ref_id=config["reference"]["id"],
    shell:
        r"""
        # Track the stitched ancestral alignment explicitly so that
        # ancestral_stitched reruns whenever the genotype-specific inputs change.
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
        ' {input.alns} > "{output.alignment}"
        """

rule ancestral_stitched:
    """Infer ancestral states on the stitched tree using the combined non-masked alignment. Mutations are relative to the configured reference sequence."""
    input:
        tree = STITCHED_DIR + "/{gene}_tree.nwk",
        alignment = STITCHED_DIR + "/combined_nonmasked_aln.fasta",
        annotation= config["reference"]["gff"],
        translations=expand("data/ancestral_translations/cds_{g}.fasta", g=config["genes"]),
        root = config["reference"]["fasta"],

    output:
        node_data = STITCHED_DIR + "/node_data/{gene}_muts.json",
    params:
        genes=" ".join(ANCESTRAL_GENES),
        translation_pattern="data/ancestral_translations/cds_%GENE.fasta",

    shell:
        r"""
        augur ancestral --tree {input.tree} --alignment {input.alignment} \
                        --translations {params.translation_pattern} --genes {params.genes} \
                        --output-node-data {output.node_data} \
                        --annotation {input.annotation}  \
                        --root-sequence {input.root}
        """
