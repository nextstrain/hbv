"""
Rules for assigning clades and subclades to stitched and non-stitched trees.
"""

rule clades_stitched:
    """Assign major clades on the stitched global trees."""
    input:
        tree = STITCHED_DIR + "/{gene}_tree.nwk",
        metadata = "data/filtered/metadata.len_filtered.tsv",
        script = "scripts/assign_clades.py"
    output:
        clades = STITCHED_DIR + "/{gene}_clades.json",
    params:
         min_count = config["clade_settings"]["min_count_per_clade"],
         min_count_mode = config["clade_settings"]["min_count_mode"],
         only_monophyletic_flag = "--only_monophyletic" if config["clade_settings"].get("only_monophyletic", True) else "" ,
         branch_display_flag = "--branch_display_only_largest" if config["clade_settings"].get("branch_display_only_largest",False) else "",
         genotype_col = "genotype_genbank",
         nextclade_label = "HBV Type",
    shell:
        r"""
        python3 {input.script} \
          --tree {input.tree} \
          --metadata {input.metadata} \
          --output {output.clades} \
          --min-count {params.min_count} \
          --min-count-mode {params.min_count_mode} \
          --subtype-col {params.genotype_col} \
          --nextclade-label "{params.nextclade_label}" \
          {params.only_monophyletic_flag} \
          {params.branch_display_flag}
         """

rule sub_clades_stitched:
    """Assign subclades on the stitched global trees, with major clades as fallback annotations."""
    input:
        tree = STITCHED_DIR + "/{gene}_tree.nwk",
        metadata = "data/filtered/metadata.len_filtered.tsv",
        fallback_clades = STITCHED_DIR + "/{gene}_clades.json",
        script = "scripts/assign_clades.py"
    output:
        clades = STITCHED_DIR + "/{gene}_subclades.json",
    params:
         min_count = config["clade_settings"]["min_count_per_subclade"],
         min_count_mode = config["clade_settings"]["min_count_mode"],
         only_monophyletic_flag = "--only_monophyletic" if config["clade_settings"].get("only_monophyletic", True) else "" ,
         branch_display_flag = "--branch_display_only_largest" if config["clade_settings"].get("branch_display_only_largest",False) else "",
         genotype_col = "subgenotype_genbank",
         fallback_genotype_col = "genotype_genbank", # in case subgenotype information is missing, use genotype information to assign to major clade
         nextclade_label = "HBV Subtype",
         fallback_annotation_nextclade_label = "HBV Type", # if subtype information is missing and genotype information is used as fallback, use the major clade as annotation for nextclade clades as well
    shell:
        r"""
        python3 {input.script} \
          --tree {input.tree} \
          --metadata {input.metadata} \
          --output {output.clades} \
          --min-count {params.min_count} \
          --min-count-mode {params.min_count_mode} \
          --subtype-col {params.genotype_col} \
          --nextclade-label "{params.nextclade_label}" \
          --fallback-genotype-col {params.fallback_genotype_col} \
          --fallback-annotation-file {input.fallback_clades} \
          --fallback-annotation-nextclade-label "{params.fallback_annotation_nextclade_label}" \
          {params.only_monophyletic_flag} \
          {params.branch_display_flag}
        """

rule clades_non_stitched:
    """Assign major clades for non-stitched builds."""
    input:
        tree = RESULTS + "/{mode}/{key}/{gene}_masked/{gene}_masked_refined.tree.nwk",
        metadata = "data/filtered/metadata.len_filtered.tsv",
        script = "scripts/assign_clades.py"
    output:
        clades = RESULTS + "/{mode}/{key}/{gene}_masked/{gene}_clades.json"
    params:
        min_count = config["clade_settings"]["min_count_per_clade"],
        min_count_mode = config["clade_settings"]["min_count_mode"],
        only_monophyletic_flag = "--only_monophyletic" if config["clade_settings"].get("only_monophyletic", True) else "",
        branch_display_flag = "--branch_display_only_largest" if config["clade_settings"].get("branch_display_only_largest", False) else "",
        genotype_col = "genotype_genbank",
        nextclade_label = "HBV Type",
    shell:
        r"""
        python3 {input.script} \
          --tree {input.tree} \
          --metadata {input.metadata} \
          --output {output.clades} \
          --min-count {params.min_count} \
          --min-count-mode {params.min_count_mode} \
          --subtype-col {params.genotype_col} \
          --nextclade-label "{params.nextclade_label}" \
          {params.only_monophyletic_flag} \
          {params.branch_display_flag}
        """

rule sub_clades_non_stitched:
    """Assign subclades for non-stitched builds, with major clades as fallback annotations."""
    input:
        tree = RESULTS + "/{mode}/{key}/{gene}_masked/{gene}_masked_refined.tree.nwk",
        metadata = "data/filtered/metadata.len_filtered.tsv",
        fallback_clades = RESULTS + "/{mode}/{key}/{gene}_masked/{gene}_clades.json",
        script = "scripts/assign_clades.py"
    output:
        clades = RESULTS + "/{mode}/{key}/{gene}_masked/{gene}_subclades.json"
    params:
        min_count = config["clade_settings"]["min_count_per_subclade"],
        min_count_mode = config["clade_settings"]["min_count_mode"],
        only_monophyletic_flag = "--only_monophyletic" if config["clade_settings"].get("only_monophyletic", True) else "",
        branch_display_flag = "--branch_display_only_largest" if config["clade_settings"].get("branch_display_only_largest", False) else "",
        genotype_col = "subgenotype_genbank",
        fallback_genotype_col = "genotype_genbank",
        nextclade_label = "HBV Subtype",
        fallback_annotation_nextclade_label = "HBV Type",
    shell:
        r"""
        python3 {input.script} \
          --tree {input.tree} \
          --metadata {input.metadata} \
          --output {output.clades} \
          --min-count {params.min_count} \
          --min-count-mode {params.min_count_mode} \
          --subtype-col {params.genotype_col} \
          --nextclade-label "{params.nextclade_label}" \
          --fallback-genotype-col {params.fallback_genotype_col} \
          --fallback-annotation-file {input.fallback_clades} \
          --fallback-annotation-nextclade-label "{params.fallback_annotation_nextclade_label}" \
          {params.only_monophyletic_flag} \
          {params.branch_display_flag}
        """
