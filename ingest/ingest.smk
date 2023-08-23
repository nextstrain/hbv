rule vendor_nextclade3_x86:
    """
    Nextclade v3 is unreleased. This repo includes an arch64 version but for GitHub
    actions we need an x86 version. 
    This executable is produced as part of Nextclade's CI - e.g. see artifacts from
    https://github.com/nextstrain/nextclade/actions/runs/5829996974

    This rule is temporary as once Nextclde v3 is released we can either use the
    version in the Nextstrain runtime or use an approach similar to ncov-ingest:
    https://github.com/nextstrain/ncov-ingest/blob/6c7ae5d7d259d62b8a0f016f02be51a11090592a/workflow/snakemake_rules/nextclade.smk#L130
    """
    params:
        path = config['nextclade_binary']
    shell:
        """
        curl -fsSL "https://nextstrain-scratch.s3.amazonaws.com/james/nextclade-x86_64-unknown-linux-gnu" -o {params[0]:q}
        chmod +x {params[0]:q}
        NEXTCLADE_VERSION="$(./{params[0]:q} --version)"
        echo "[ INFO] Nextclade version: $NEXTCLADE_VERSION" 
        """

rule fetch_genbank:
    params:
        term = 'Hepatitis+B+virus[All+Fields]complete+genome[All+Fields]'
    output:
        genbank = "ingest/data/genbank.gb",
    retries: 1  # Requires snakemake 7.7.0 or later
    shell:
        """
        ingest/vendored/fetch-from-ncbi-entrez --term {params.term:q} --output {output.genbank}
        """

rule add_extra_genomes:
    """
    This step shouldn't be necessary but the NCBI reference genome, NC_003977,
    is not returned via the ENTREZ query. Even if we change the reference it's good
    to add this genome.
    """
    input:
        entrez = "ingest/data/genbank.gb",
        ref = "ingest/config/NC_003977.gb"
    output:
        genbank = "ingest/data/genbank.extended.gb",
    shell:
        """
        cat {input.ref:q} {input.entrez:q} > {output.genbank:q}
        """

rule parse_genbank:
    input:
        genbank = "ingest/data/genbank.extended.gb",
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
        reference = config['reference_accession']
    shell:
        """
        ingest/scripts/re-circularise.py \
            --seqs-in {input.sequences} --meta-in {input.metadata} \
            --seqs-out {output.sequences} --meta-out {output.metadata} \
            --reference {params.reference}
        """

rule align_everything:
    """
    Nextclade v3 is used to align all genomes using a reference dataset.
    Note that the minimum seed match rate is specified in the dataset itself.
    """
    input:
        sequences = "ingest/results/genbank.recircular.fasta",
        metadata = "ingest/results/genbank.recircular.tsv",
    output:
        # Note that these outputs require `--output-basename nextclade`
        alignment = "ingest/results/nextclade.aligned.fasta",
        summary = "ingest/results/nextclade.tsv",
        translations = expand("ingest/results/nextclade_gene_{gene}.translation.fasta", gene=config['genes'])
    params:
        dataset = config['nextclade_dataset'],
        nextclade = config['nextclade_binary']
    threads: 4
    shell:
        """
        {params.nextclade} run \
            -j {threads} --silent --replace-unknown \
            --input-dataset {params.dataset} \
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
        metadata_columns = ['name', 'accesion', "strain_name", "date", "year", "region", "country", "host", "genotype_genbank", "subgenotype_genbank", \
        "circularise", "circularise_shift_bp","clade_nextclade","QC_overall_score","QC_overall_status","total_substitutions","total_deletions", \
        "total_insertions","total_frame_shifts","total_missing","alignment_score","coverage","QC_missing_data","QC_mixed_sites","QC_rare_mutations", \
        "QC_frame_shifts","QC_stop_codons"]
    shell:
        """
        ingest/scripts/tsv-to-ndjson.py < {input.metadata} |
            ingest/scripts/fix_country_field.py |
            ingest/vendored/apply-geolocation-rules --geolocation-rules ingest/config/geoLocationRules.tsv |
            ingest/scripts/add-year.py |
            ingest/scripts/ndjson-to-tsv.py --metadata-columns {params.metadata_columns} --metadata {output.metadata}
        """

rule copy_ingest_files:
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