rule fetch_genbank:
    params:
        term = 'Hepatitis+B+virus[All+Fields]complete+genome[All+Fields]'
    output:
        genbank = "ingest/data/genbank.gb",
    retries: 1  # Requires snakemake 7.7.0 or later
    shell:
        """
        ingest/scripts/fetch-genbank.py --term {params.term:q} --output {output.genbank}
        """

rule parse_genbank:
    input:
        genbank = "ingest/data/genbank.gb",
    output:
        sequences = "ingest/results/genbank.fasta",
        metadata = "ingest/results/genbank.tsv",
    shell:
        """
        ingest/scripts/parse-genbank.py --genbank {input.genbank} \
            --output-sequences {output.sequences} --output-metadata {output.metadata}
        """

rule re_circularise:
    input:
        metadata = "ingest/results/genbank.tsv",
        sequences = "ingest/results/genbank.fasta",
    output:
        sequences = "ingest/results/genbank.recircular.fasta",
        metadata = "ingest/results/genbank.recircular.tsv",
    params:
        reference = "JN182318" # should be the same as the one used for the big alignment
    shell:
        """
        ingest/scripts/re-circularise.py \
            --seqs-in {input.sequences} --meta-in {input.metadata} \
            --seqs-out {output.sequences} --meta-out {output.metadata} \
            --reference {params.reference}
        """

"""
We use nextclade to align to reference JN182318 and also to call "genotype" ("clade" in the output CSV)
See early commits (e.g. https://github.com/nextstrain/hepatitisB/tree/93cbb184d329c12835d0ffedec61e0fbd8cb53db)
for a version of this rule which used nextalign instead.
Note that the minimum seed match rate is specified in the dataset itself
"""
rule align_everything:
    input:
        sequences = "ingest/results/genbank.recircular.fasta",
        metadata = "ingest/results/genbank.recircular.tsv",
        reference = "ingest/config/JN182318.fasta"
    output:
        alignment = "ingest/results/nextclade.aligned.fasta",
        summary = "ingest/results/nextclade.tsv"
    threads: 4
    shell:
        """
        nextclade run \
            -j {threads} --silent --replace-unknown \
            --input-dataset nextclade_datasets/references/JN182318/versions/2023-06-26 \
            --output-all ingest/results \
            --output-basename nextclade \
            {input.sequences}
        """

rule join_nextclade_metadata:
    input:
        metadata = "ingest/results/genbank.recircular.tsv",
        nextclade = "ingest/results/nextclade.tsv"
    output:
        metadata = "ingest/results/metadata_nextclade.tsv",
        summary = "ingest/results/metadata.summary.txt",
    shell:
        """
        ingest/scripts/join-nextclade-metadata.py \
             --metadata {input.metadata} --nextclade {input.nextclade} \
             --output {output.metadata} --summary {output.summary}
        """

rule transform_metadata:
    input:
        metadata = "ingest/results/metadata_nextclade.tsv"
    output:
        metadata = "ingest/results/metadata.tsv"
    params:
        metadata_columns = ['name', 'accesion', "strain_name", "region", "country", "host", "genotype_genbank", "subgenotype_genbank", "collection_date", \
        "circularise", "circularise_shift_bp","clade_nextclade","QC_overall_score","QC_overall_status","total_substitutions","total_deletions", \
        "total_insertions","total_frame_shifts","total_missing","alignment_score","coverage","QC_missing_data","QC_mixed_sites","QC_rare_mutations", \
        "QC_frame_shifts","QC_stop_codons"]
    shell:
        """
        ingest/scripts/tsv-to-ndjson.py < {input.metadata} |
            ingest/scripts/fix_country_field.py |
            ingest/scripts/apply-geolocation-rules.py --geolocation-rules ingest/config/geoLocationRules.tsv |
            ingest/scripts/ndjson-to-tsv.py --metadata-columns {params.metadata_columns} --metadata {output.metadata}
        """

rule provision:
    input:
        sequences = "ingest/results/genbank.recircular.fasta",
        aligned = "ingest/results/nextclade.aligned.fasta"
    output:
        sequences = "ingest/results/sequences.fasta",
        aligned = "ingest/results/aligned.fasta"
    shell:
        """
        cp {input.sequences} {output.sequences}
        cp {input.aligned} {output.aligned}
        """