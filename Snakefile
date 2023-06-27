
include: "ingest/ingest.smk"

BUILDS = ["all", "dev", "nextclade-tree"]

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


### TODO XXX
# iqtree writes log files etc to the input path + some suffix


def filter_params(wildcards):
    if wildcards.build == "all":
        return ""
    elif wildcards.build == "dev":
        return "--group-by genotype_genbank --subsample-max-sequences 500"
    elif wildcards.build == "nextclade-tree":
        return "--group-by genotype_genbank --subsample-max-sequences 2000"
    elif wildcards.build == "nextclade-sequences":
        return "--group-by genotype_genbank --subsample-max-sequences 25"
    raise Exception("Unknown build parameter")


rule filter:
    input:
        sequences = "ingest/results/aligned.fasta",
        metadata = "ingest/results/metadata.tsv",
        include = "config/include.txt",
    output:
        sequences = "results/{build}/filtered.fasta"
    params:
        args = filter_params
    shell:
        """
        augur filter \
            --sequences {input.sequences} --metadata {input.metadata} \
            --include {input.include} \
            {params.args} \
            --output {output.sequences}
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

rule refine:
    message:
        """
        Refining tree
        NO TIMETREE CURRENTLY
        """
    input:
        tree = "results/{build}/tree_raw.nwk",
        alignment = "results/{build}/filtered.fasta",
        metadata = "ingest/results/metadata.tsv"
    output:
        tree = "results/{build}/tree.nwk",
        node_data = "results/{build}/branch_lengths.json"
    params:
        coalescent = "opt",
        date_inference = "marginal",
        clock_filter_iqd= 4
    shell:
        """
        augur refine \
            --tree {input.tree} \
            --alignment {input.alignment} \
            --metadata {input.metadata} \
            --output-tree {output.tree} \
            --output-node-data {output.node_data} \
            --root HQ603073 \
            --coalescent {params.coalescent} \
            --date-confidence \
            --date-inference {params.date_inference} \
            --clock-filter-iqd {params.clock_filter_iqd}
        """


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

rule export:
    input:
        tree = "results/{build}/tree.nwk",
        metadata = "ingest/results/metadata.tsv",
        node_data = node_data_files,
        auspice_config = "config/auspice_config_all.json", ### TODO XXX parameterise when necessary
    output:
        auspice_json = "auspice/hbv_{build}.json"
    shell:
        """
        export AUGUR_RECURSION_LIMIT=10000;
        augur export v2 \
            --tree {input.tree} \
            --metadata {input.metadata} \
            --node-data {input.node_data} \
            --auspice-config {input.auspice_config} \
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
