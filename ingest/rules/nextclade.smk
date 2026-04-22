"""
Rules for running Nextclade on curated ingest sequences and merging its outputs.

Required external inputs:
- `results/sequences.fasta`
- `data/circularised/circularised_metadata.tsv`
- the local Nextclade dataset configured in `defaults/config.yaml`

Key outputs:
- `data/nextclade/translations/`
- `results/metadata.tsv`
- `results/aligned.fasta`
- `data/qc/metadata_summary.txt`

See the Nextclade CLI docs for details on configurable inputs and outputs:
https://docs.nextstrain.org/projects/nextclade/page/user/nextclade-cli.html
"""

# Future improvement: fetch the Nextclade dataset directly instead of relying on a local dataset path.

rule nextclade:
    """
    Nextclade v3 is used to align all genomes using a reference dataset, perform QC and infer genotypes ("clade_nextclade")
    We can output a preliminary tree here but do not need to for this pipeline.
    Note that the minimum seed match rate is specified in the dataset itself.
    Note that QC metrics are not used for filtering yet.
    """
    input:
        sequences = "results/sequences.fasta",
    output:
        alignment = temp("data/nextclade/aligned_sequences.fasta"),
        translations_snakemake = expand("data/nextclade/translations/cds_{gene}.fasta", gene=config['genes']),
        metadata = "data/nextclade/nextclade_metadata.tsv",
    params:
        dataset = config['nextclade_dataset'],
        translations_pattern = lambda w: "data/nextclade/translations/cds_{cds}.fasta",
    threads: 4
    shell:
        """
        nextclade run \
            -j {threads} --silent --replace-unknown \
            --input-dataset {params.dataset} \
            --output-fasta {output.alignment} \
            --output-translations {params.translations_pattern} \
            --output-tsv {output.metadata} \
            {input.sequences}
        """

rule join_nextclade_metadata:
    """Join Nextclade outputs back onto the circularised ingest metadata."""
    input:
        metadata = "data/circularised/circularised_metadata.tsv",
        nextclade = "data/nextclade/nextclade_metadata.tsv"
    output:
        metadata = "results/metadata.tsv",
        summary = "data/qc/metadata_summary.txt",
    shell:
        """
        scripts/join-nextclade-metadata.py \
             --metadata {input.metadata} --nextclade {input.nextclade} \
             --output {output.metadata} --summary {output.summary}
        """

rule copy_ingest_alignment:
    """Promote the Nextclade alignment to the final ingest results directory."""
    input:
        aligned = "data/nextclade/aligned_sequences.fasta",
    output:
        aligned = "results/aligned.fasta"
    shell:
        """
        cp {input.aligned} {output.aligned}
        """
