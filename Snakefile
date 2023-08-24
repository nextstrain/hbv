configfile: "config/config.yaml"
include: "ingest/ingest.smk"

BUILDS = ["all", "dev", "nextclade-tree"]

DEFAULT_TARGETS = ["all", *config['genotypes']]

rule all:
    input:
        auspice_json = expand("auspice/hbv_{target}.json", target=DEFAULT_TARGETS)

def get_root(wildcards):
    if wildcards.build in config['roots']:
        return config['roots'][wildcards.build]
    return config['roots']['all']

## For "all-HBV" like builds - i.e. those which don't downsample to a single genotype
## We deliberatly include some outliers so we can prune out non-human sequences
## And we include the reference genome
rule include_file:
    output:
        file = "results/{build}/include.txt",
    params:
        root = get_root,
        outgroups = lambda w: config['outgroups'] if w.build in BUILDS else [],
        ref = lambda w: [config['reference_accession']] if w.build in BUILDS else [],
    run:
        with open(output[0], 'w') as fh:
            print(params[0], file=fh)
            for name in params[1]:
                print(name, file=fh)
            for name in params[2]:
                print(name, file=fh)

## TODO - there are a number of nextclade QC status' we can filter on here.
## Currently the settings in the nextclade dataset need to be looked at as 
## around 40% of all sequences (including the entirety of some genotypes)
## have QC=bad mainly due to frameshifts and stop codons.
def filter_params(wildcards):
    if wildcards.build == "all":
        return ""
    elif wildcards.build == "dev":
        return "--group-by genotype_genbank --subsample-max-sequences 500"
    elif wildcards.build == "nextclade-tree":
        return "--group-by genotype_genbank --subsample-max-sequences 2000"
    elif wildcards.build == "nextclade-sequences":
        return "--group-by genotype_genbank --subsample-max-sequences 25"
    elif wildcards.build in config['genotypes']:
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

rule prune_outgroup:
    input:
        tree = "results/{build}/tree.refined.nwk"
    output:
        tree = "results/{build}/tree.nwk"
    params:
        outgroups = lambda w: [get_root(w), *(config['outgroups'] if w.build in BUILDS else [])]
    shell:
        """
        python scripts/remove_outgroups.py --names {params.outgroups} --tree {input.tree} --output {output.tree}
        """


rule ancestral:
    input:
        tree = "results/{build}/tree.nwk",
        alignment = "results/{build}/filtered.fasta",
        translations = rules.align_everything.output.translations,
    output:
        node_data = "results/{build}/mutations.json",
        translations = expand("results/{{build}}/mutations_{gene}.translation.fasta", gene=config['genes'])
    params:
        inference = "joint",
        genes = config['genes'],
        annotation = config['temporary_genemap_for_augur_ancestral'],
        translation_pattern = 'results/{build}/mutations_%GENE.translation.fasta'
    shell:
        """
        augur ancestral \
            --tree {input.tree} \
            --alignment {input.alignment} \
            --inference {params.inference} \
            --annotation {params.annotation} \
            --translations ingest/results/nextclade_gene_%GENE.translation.fasta \
            --genes {params.genes} \
            --output-node-data {output.node_data} \
            --output-translations {params.translation_pattern}
        """

## Clades are hard to predict using mutational signatures for HBV due to the amount
## of recombination and the inherit difference in tree reconstruction.
## They are only run for the nextclade dataset build, which we manually check
rule clades:
    message: "Labeling clades as specified in config/clades.tsv"
    input:
        tree =  "results/{build}/tree.nwk",
        # aa_muts = rules.translate.output.node_data,
        nuc_muts = "results/{build}/mutations.json",
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

# make sure all differences between the alignment reference and the root are attached as mutations to the root.
# This is almost possible via `augur ancestral --root-sequence` however that functionality does not work for
# CDSs which wrap the origin.
rule attach_root_mutations:
    input:
        mutations = "results/{build}/mutations.json",
        tree = "results/{build}/tree.nwk",
        translations = expand("results/{{build}}/mutations_{gene}.translation.fasta", gene=config['genes'])
    output:
        mutations = "results/{build}/mutations.reference-mutations-on-root.json",
    params:
        genes = config['genes'],
        reference = config['reference_accession'],
        translation_pattern = 'results/{build}/mutations_%GENE.translation.fasta',
    shell:
        """
        python3 scripts/attach_root_mutations.py \
            --tree {input.tree} \
            --reference-name {params.reference} \
            --translations {params.translation_pattern} \
            --genes {params.genes} \
            --mutations-in {input.mutations} \
            --mutations-out {output.mutations}
        """

def node_data_files(wildcards):
    patterns = [
        "results/{build}/branch_lengths.json",
    ]

    # For nextclade dataset trees we need to annotate all differences between the root & the nextclade reference
    # But for everything else we don't want this!
    if wildcards.build == "nextclade-tree":
        patterns.append("results/{build}/mutations.reference-mutations-on-root.json")
    else:
        patterns.append("results/{build}/mutations.json")

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

rule auspice_config:
    output:
        "results/{build}/auspice_config.json"
    run:
        build_name = wildcards.build
        title = "Genetic diversity of Hepatitis B virus"
        if build_name in config['genotypes']:
            title = f"Genetic diversity of Hepatitis B virus genotype {build_name}" 

        dev_colorings = [
            { "key": "circularise", "title": "genome rotated?", "type": "categorical" },
            { "key": "circularise_shift_bp", "title": "rotation offset", "type": "continuous" },
            { "key": "QC_overall_score",     "title": "Nextclade QC overall score", "type": "categorical" },
            { "key": "QC_overall_status",    "title": "Nextclade QC overall status", "type": "categorical" },
            { "key": "total_substitutions",  "title": "Nextclade QC total substitutions", "type": "continuous" },
            { "key": "total_deletions",      "title": "Nextclade QC total deletions", "type": "continuous" },
            { "key": "total_insertions",     "title": "Nextclade QC total insertions", "type": "continuous" },
            { "key": "total_frame_shifts",   "title": "Nextclade QC total frame shifts", "type": "continuous" },
            { "key": "total_missing",        "title": "Nextclade QC total missing", "type": "continuous" },
            { "key": "alignment_score",      "title": "Nextclade QC aln score", "type": "continuous" },
            { "key": "coverage",             "title": "Nextclade QC coverage", "type": "continuous" },
            { "key": "QC_missing_data",      "title": "Nextclade QC status missing data", "type": "categorical" },
            { "key": "QC_mixed_sites",       "title": "Nextclade QC status mixed sites", "type": "categorical" },
            { "key": "QC_rare_mutations",    "title": "Nextclade QC status rare muts", "type": "categorical" },
            { "key": "QC_frame_shifts",      "title": "Nextclade QC status frame shifts", "type": "categorical" },
            { "key": "QC_stop_codons",       "title": "Nextclade QC status stop codons", "type": "categorical" },
        ] if build_name in ['dev', 'all'] else []

        augur_clades =     { "key": "clade_membership", "title": "HBV Genotype (clades TSV)", "type": "categorical" }
        nextclade_clades = { "key": "clade_nextclade",  "title": "HBV Genotype (Nextclade)",  "type": "categorical" }


        default_color_by = "region" # genotype builds are the ~same genotype...
        if build_name=='nextclade-tree':
            default_color_by = "clade_membership" # augur clades TSV
        elif build_name in ['dev', 'all']:
            default_color_by = "clade_nextclade"

        data = {
            "title": title,
            "build_url": "https://github.com/nextstrain/hepatitisB/",
            "data_provenance": [
                {"name": "GenBank", "url": "https://www.ncbi.nlm.nih.gov/nucleotide/"}
            ],
            "maintainers": [
                {"name": "Katie Kistler"},
                {"name": "James Hadfield"},
                {"name": "Nextstrain"},
            ],
            "colorings": [
                {"key": "gt", "title": "Genotype (by position)", "type": "categorical"},
                {"key": "region", "title": "region", "type": "categorical"},
                { "key": "country", "title": "country", "type": "categorical" },
                { "key": "host", "title": "host", "type": "categorical" },
                (augur_clades if build_name in ['dev', "nextclade-tree"] else None),
                (nextclade_clades if build_name != 'nextclade-tree' else None),
                { "key": "genotype_genbank",    "title": "HBV Genotype (Genbank)",    "type": "categorical" },
                { "key": "subgenotype_genbank", "title": "HBV Subgenotype (Genbank)", "type": "categorical" },
                { "key": "year", "title": "collection year", "type": "continuous"},
                { "key": "date", "title": "collection date", "type": "categorical" },
                *dev_colorings,
                { "key": "strain_name", "title": "Strain (name)", "type": "categorical" },
            ],
            "display_defaults": {
                "color_by": default_color_by,
                "branch_label": "augur_clades", # This won't be there for most builds...
                "map_triplicate": True
            },
            "filters": [
                "host",
                "country",
                "strain_name",
                (augur_clades['key'] if build_name in ['dev', "nextclade-tree"] else None),
                "clade_nextclade",
                "genotype_genbank",
                "subgenotype_genbank",
                *[d['key'] for d in dev_colorings if d['type']=='categorical'],
            ],
            "geo_resolutions": [
                "region",
                "country"
            ],
            "panels": [
                "tree",
                "map",
                "entropy"
            ]
        }

        ## Prune out empty values
        data['colorings'] = [c for c in data['colorings'] if c]
        data['filters'] = [f for f in data['filters'] if f]

        with open(output[0], 'w') as fh:
            json.dump(data, fh, indent=2)

# Note that `--validation-mode skip` is used because the JSON schema is overly restrictive
# and rejects our branch labels as they have underscores
rule export:
    input:
        tree = "results/{build}/tree.nwk",
        metadata = "ingest/results/metadata.tsv",
        node_data = node_data_files,
        auspice_config = "results/{build}/auspice_config.json",
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

rule deploy:
    input:
        *rules.all.input,
    output:
        touch(expand("results/{target}/deploy.done", target=DEFAULT_TARGETS))
    params:
        deploy_destination = config["deploy_destination"]
    shell:
        """
        nextstrain remote upload {params.deploy_destination} {input}
        """
