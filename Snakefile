
include: "ingest/ingest.smk"

BUILDS = ["all", "dev", "nextclade-tree"]

GENOTYPES = ["A", "B", "C", "D"]

ROOT = {
    "all": "HQ603073", # NHP-HBV isolate
    # genotype roots chosen by examining the entire tree and picking a suitably close isolate
    "A": "MK534669", # root is genotype I (I is A/C/G recombinant)
    "B": "MK534669", # root is genotype I (I is A/C/G recombinant)
    "C": "MK534669", # root is genotype I (I is A/C/G recombinant)
    "D": "KX186584", # root is genotype E
}

# rule all:
#     input:
#         auspice_json = expand("auspice/hepatitisB_{lineage}.json", lineage=LINEAGES)

# rule files:
#     params:
#         reference = "config/reference_hepatitisB_{lineage}.gb",
#         auspice_config = "config/auspice_config_all.json",
#         dropped_strains = "config/dropped_strains_hepatitisB_{lineage}.txt"


# files = rules.files.params



# rule filter:
#     message:
#         """
#         Filtering to
#           - {params.sequences_per_group} sequence(s) per {params.group_by!s}
#           - minimum genome length of {params.min_length}
#         """
#     input:
#         sequences = rules.parse.output.sequences,
#         metadata = rules.parse.output.metadata,
#         exclude = files.dropped_strains,
#         include = "config/include.txt", # TESTING ONLY
#     output:
#         sequences = "results/filtered_hepatitisB_{lineage}.fasta"
#     params:
#         group_by = "country year",
#         sequences_per_group = 50,
#         min_length = 3000
#     shell:
#         """
#         augur filter \
#             --sequences {input.sequences} \
#             --metadata {input.metadata} \
#             --exclude {input.exclude} \
#             --exclude-all --include {input.include} \
#             --output {output.sequences}
#         """

# --group-by {params.group_by} \
# --sequences-per-group {params.sequences_per_group} \
# --min-length {params.min_length}

# rule align:
#     message:
#         """
#         Aligning sequences to {input.reference}
#           - filling gaps with N
#         """
#     input:
#         sequences = rules.filter.output.sequences,
#         reference = files.reference
#     threads: 8
#     output:
#         alignment = "results/aligned_hepatitisB_{lineage}.fasta"
#     shell:
#         """
#         augur align \
#             --sequences {input.sequences} \
#             --output {output.alignment} \
#             --reference-name JN182318 \
#             --nthreads {threads} \
#             --fill-gaps
#         """

#                   #  --reference-sequence {input.reference} \
# #

def get_root(wildcards):
    if wildcards.build in ROOT:
        return ROOT[wildcards.build]
    return ROOT['all']

rule include_file:
    output:
        file = "results/{build}/include.txt",
    params:
        root = get_root
    shell:
        """
        echo {params.root:q} > {output.file}
        """

def filter_params(wildcards):
    if wildcards.build == "all":
        return ""
    elif wildcards.build == "dev":
        return "--group-by genotype_genbank --subsample-max-sequences 500"
    elif wildcards.build == "nextclade-tree":
        return "--group-by genotype_genbank --subsample-max-sequences 2000"
    elif wildcards.build == "nextclade-sequences":
        return "--group-by genotype_genbank --subsample-max-sequences 25"
    elif wildcards.build in GENOTYPES:
        if wildcards.build == "C":
            query = f"--query \"(clade_nextclade=='C') | (clade_nextclade=='C_re')\""
        else:
            query = f"--query \"clade_nextclade=='{wildcards.build}'\""
        return f"{query} --group-by year --subsample-max-sequences 500"
    raise Exception("Unknown build parameter")

rule filter:
    input:
        sequences = "ingest/results/aligned.fasta",
        metadata = "ingest/results/metadata.tsv",
        include = "results/{build}/include.txt",
        exclude = "config/exclude.txt",
    output:
        sequences = "results/{build}/filtered.fasta",
        metadata = "results/{build}/filtered.tsv",
    params:
        args = filter_params
    shell:
        """
        augur filter \
            --sequences {input.sequences} --metadata {input.metadata} \
            --include {input.include} --exclude {input.exclude} \
            {params.args} \
            --output-sequences {output.sequences} --output-metadata {output.metadata}
        """


rule tree:
    message: "Building tree"
    input:
        alignment = "results/{build}/filtered.fasta"
    output:
        tree = "results/{build}/tree_raw.nwk"
    threads: 8
    shell:
        """
        augur tree \
            --alignment {input.alignment} \
            --nthreads {threads} \
            --output {output.tree}
        """

def refine_parameters(wildcards):
    params = []
    ## Following code will be useful if we infer a timetree, which _may_
    ## be possible at the subgenotype level?
    # params.append("--timetree")
    # params.append("--date-inference marginal")
    # params.append("--date-confidence")
    # params.append("--stochastic-resolve")
    # params.append("--coalescent opt")
    # params.append("--clock-filter-iqd 4")
    # params.append("--root best")
    params.append(f"--root {get_root(wildcards)}")
    return " ".join(params)

rule refine:
    message: "Pruning the outgroup from the tree"
    input:
        tree = "results/{build}/tree_raw.nwk",
        alignment = "results/{build}/filtered.fasta",
        metadata = "ingest/results/metadata.tsv"
    output:
        tree = "results/{build}/tree.refined.nwk",
        node_data = "results/{build}/branch_lengths.json"
    params:
        refine = refine_parameters,
    shell:
        """
        augur refine \
            --tree {input.tree} --alignment {input.alignment} --metadata {input.metadata} \
            {params.refine} \
            --output-tree {output.tree} --output-node-data {output.node_data}
        """


## see https://github.com/nextstrain/augur/issues/340#issuecomment-545184212
rule prune_outgroup:
    input:
        tree = "results/{build}/tree.refined.nwk"
    output:
        tree = "results/{build}/tree.nwk"
    params:
        root = get_root
    run:
        from Bio import Phylo
        T = Phylo.read(input[0], "newick")
        outgroup = [c for c in T.find_clades() if str(c.name) == params[0]][0]
        T.prune(outgroup)
        Phylo.write(T, output[0], "newick")


rule ancestral:
    input:
        tree = "results/{build}/tree.nwk",
        alignment = "results/{build}/filtered.fasta",
    output:
        node_data = "results/{build}/nt_muts.json"
    params:
        inference = "joint"
    shell:
        """
        augur ancestral \
            --tree {input.tree} \
            --alignment {input.alignment} \
            --output-node-data {output.node_data} \
            --inference {params.inference}
        """

# rule translate:
#     message: "Translating amino acid sequences"
#     input:
#         tree = rules.refine.output.tree,
#         node_data = rules.ancestral.output.node_data,
#         reference = files.reference
#     output:
#         node_data = "results/aa_muts_hepatitisB_{lineage}.json"
#     shell:
#         """
#         augur translate \
#             --tree {input.tree} \
#             --ancestral-sequences {input.node_data} \
#             --reference-sequence {input.reference} \
#             --output {output.node_data} \
#         """


## Clades are hard to predict using mutational signatures for HBV due to the amount
## of recombination and the inherit difference in tree reconstruction.
## They are only run for the nextclade dataset build, which we manually check
rule clades:
    message: "Labeling clades as specified in config/clades.tsv"
    input:
        tree =  "results/{build}/tree.nwk",
        # aa_muts = rules.translate.output.node_data,
        nuc_muts = "results/{build}/nt_muts.json",
        clades = "config/clades-genotypes.tsv"
    output:
        clade_data = "results/{build}/clades-genotypes.json",
    shell:
        """
        augur clades \
            --tree {input.tree} \
            --mutations {input.nuc_muts} \
            --membership-name clade_membership \
            --label-name augur_clades \
            --clades {input.clades} \
            --output {output.clade_data}
        """

def node_data_files(wildcards):
    patterns = [
        "results/{build}/branch_lengths.json",
        "results/{build}/nt_muts.json",
        ## TODO XXX translated AA json
    ]

    # Only infer genotypes via `augur clades` for manually curated nextclade dataset builds
    if wildcards.build == "nextclade-tree" or wildcards.build == "dev":
        patterns.append("results/{build}/clades-genotypes.json")

    inputs = [f.format(**dict(wildcards)) for f in patterns]
    return inputs

rule colors:
    input:
        ordering="config/color_ordering.tsv",
        color_schemes="config/color_schemes.tsv",
        metadata="results/{build}/filtered.tsv",
    output:
        colors="results/{build}/colors.tsv",
    shell:
        """
        python3 scripts/assign-colors.py \
            --metadata {input.metadata} \
            --ordering {input.ordering} \
            --color-schemes {input.color_schemes} \
            --output {output.colors}
        """

rule export:
    input:
        tree = "results/{build}/tree.nwk",
        metadata = "ingest/results/metadata.tsv",
        node_data = node_data_files,
        auspice_config = "config/auspice_config_all.json", ### TODO XXX parameterise when necessary
        colors = "results/{build}/colors.tsv",
        lat_longs = "config/lat-longs.tsv",
    output:
        auspice_json = "auspice/hbv_{build}.json"
    shell:
        """
        augur export v2 \
            --tree {input.tree} \
            --metadata {input.metadata} \
            --node-data {input.node_data} \
            --auspice-config {input.auspice_config} \
            --colors {input.colors} --lat-longs {input.lat_longs} \
            --output {output.auspice_json} \
            --validation-mode skip
        """

## AUGUR BUG TODO XXX
## we skip validation because
##   .display_defaults.branch_label "genotype_inferred" failed pattern validation for "^(none|[a-zA-Z0-9]+)$"

# rule clean:
#     message: "Removing directories: {params}"
#     params:
#         "results ",
#         "auspice"
#     shell:
#         "rm -rfv {params}"
