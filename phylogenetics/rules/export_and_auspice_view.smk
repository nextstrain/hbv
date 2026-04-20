rule augur_export:
    wildcard_constraints:
        gene="P|S|C|X"
    input:
        tree = RESULTS + "/{mode}/{key}/{gene}_masked/{gene}_masked_refined.tree.nwk",
        node_data = [
            RESULTS + "/{mode}/{key}/{gene}_masked/ancestral/{gene}.json",
            RESULTS + "/{mode}/{key}/{gene}_masked/{gene}_subclades.json",
            RESULTS + "/{mode}/{key}/{gene}_masked/{gene}_clades.json",
        ],
        metadata = RESULTS + "/{mode}/{key}/{gene}_masked/{gene}_metadata.pruned.tsv",
        config = config["auspice"]["auspice_config"],
    output:
        auspice_json = RESULTS + "/{mode}/{key}/{gene}_masked/{gene}.json",
    params:
        colours = " \\\n            ".join(config["auspice"]["color_by_metadata"])
    threads: 1
    shell:
        """
        augur export v2 \
        --tree "{input.tree}" \
        --node-data {input.node_data} \
        --metadata "{input.metadata}" \
        --metadata-id-columns accession \
        --auspice-config "{input.config}" \
        --color-by-metadata {params.colours}\
        --minify-json \
        --output "{output.auspice_json}"
        """

rule export_stitched:
    wildcard_constraints:
        gene="P|S|C|X"
    input:
        tree = STITCHED_DIR + "/{gene}_tree.nwk",
        metadata = "data/filtered/metadata.len_filtered.tsv",
        config = config["auspice"]["auspice_config"],
        node_data = [
            STITCHED_DIR + "/node_data/{gene}_branch_lengths.json",
            STITCHED_DIR + "/node_data/{gene}_muts.json",
            STITCHED_DIR + "/{gene}_subclades.json",
            STITCHED_DIR + "/{gene}_clades.json",]
            #STITCHED_DIR + "/{gene}_node_metadata.json",

    output:
        auspice_json = STITCHED_DIR + "/{gene}.json",
    params:
        colours = " \\\n            ".join(config["auspice"]["color_by_metadata"])
    threads: 1
    shell:
        """
        augur export v2 \
            --tree "{input.tree}" \
            --node-data {input.node_data} \
            --metadata "{input.metadata}" \
            --metadata-id-columns accession \
            --output "{output.auspice_json}" \
            --auspice-config "{input.config}" \
            --color-by-metadata {params.colours}\
            --minify-json \
            --include-root-sequence-inline
        """

#______________________________________________________________________________________________________________________________________________________________________________________________

rule create_auspice_view_main_clades:
    input:
        main_clades= RESULTS + "/stitched/{gene}_global/{gene}.json",
    output:
        main_clades="auspice_datasets/{gene}_masked/main-clades.json",
    wildcard_constraints:
        gene="|".join(config["gene_mask"]),
    shell:
        r"""
        cp {input.main_clades} {output.main_clades}
        """

rule create_auspice_view_full_tree:
    input:
        full_tree= RESULTS + "/basic/all/{gene}_masked/{gene}.json",
    output:
        full_tree="auspice_datasets/{gene}_masked/full-tree.json",
    wildcard_constraints:
        gene="|".join(config["gene_mask"]),
    shell:
        r"""
        cp {input.full_tree} {output.full_tree}
        """

rule create_auspice_view_single:
    input:
        json=RESULTS + "/single-clade/{key}/{gene}_masked/{gene}.json"
    output:
        single_clades="auspice_datasets/{gene}_masked/clade_{key}.json",
    wildcard_constraints:
        gene="|".join(config["gene_mask"]),
        key="|".join(SINGLE_BUILD_GENOTYPES),
    shell:
        r"""
        cp {input.json} {output.single_clades}
        """
