"""
This part of the workflow handles running Nextclade on the curated metadata
and sequences.

REQUIRED INPUTS:

    sequences = "data/circularised.fasta",
    metadata = "data/circularised.tsv",

OUTPUTS:

    metadata = "results/metadata.tsv",
    sequences = "results/sequences.fasta",
    aligned = "results/aligned.fasta"
    summary = "data/metadata.summary.txt",
    tree = "data/nextclade/nextclade.json",



See Nextclade docs for more details on usage, inputs, and outputs if you would
like to customize the rules:
https://docs.nextstrain.org/projects/nextclade/page/user/nextclade-cli.html
"""





rule nextclade:
    """
    Nextclade v3 is used to align all genomes using a reference dataset and infer genotypes ("clade_nextclade")
    Note that the minimum seed match rate is specified in the dataset itself.
    """
    input:
        sequences = "data/circularised.fasta",
    output:
        alignment = "data/nextclade/aligned.fasta",
        tree = "data/nextclade/nextclade.json",
        translations_snakemake = expand("data/nextclade/cds_{gene}.fasta", gene=config['genes']),
        metadata = "data/nextclade/metadata.tsv",
    params:
        dataset = config['nextclade_dataset'],
        translations_pattern = lambda w: "data/nextclade/cds_{cds}.fasta",
    threads: 4
    shell:
        """
        nextclade run \
            -j {threads} --silent --replace-unknown \
            --input-dataset {params.dataset} \
            --output-fasta {output.alignment} \
            --output-translations {params.translations_pattern} \
            --output-tsv {output.metadata} \
            --output-tree {output.tree} \
            {input.sequences}
        """

rule join_nextclade_metadata:
    input:
        metadata = "data/circularised.tsv",
        nextclade = "data/nextclade/metadata.tsv"
    output:
        metadata = "results/metadata.tsv",
        summary = "data/metadata.summary.txt",
    shell:
        """
        scripts/join-nextclade-metadata.py \
             --metadata {input.metadata} --nextclade {input.nextclade} \
             --output {output.metadata} --summary {output.summary}
        """

rule copy_ingest_files:
    input:
        sequences = "data/circularised.fasta",
        aligned = "data/nextclade/aligned.fasta",
    output:
        sequences = "results/sequences.fasta",
        aligned = "results/aligned.fasta"
    shell:
        """
        cp {input.sequences} {output.sequences}
        cp {input.aligned} {output.aligned}
        """

