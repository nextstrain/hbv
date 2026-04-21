"""
This part of the workflow handles fetching sequences and metadata from NCBI.

REQUIRED INPUTS:

    None

OUTPUTS:

    ndjson = "data/raw/ncbi/ncbi_records.ndjson"
    ndjson_genbank = "data/raw/entrez/genbank_records.ndjson" # for metadata only
    active_ndjson = "data/active/ncbi_records.ndjson"

Development mode can additionally create:

    data/dev/ncbi_records.dev_sample.ndjson

There are two different approaches for fetching data from NCBI. The "Fetching from Entrez" workflow was adapted to "Fetch from NCBI" using the Mumps repo (https://github.com/nextstrain/mumps/blob/main/ingest/rules/fetch_from_ncbi.smk) as a template.
Fetching from Entrez is still included to provide *self-described* HBV (sub)genotype metadata to be compared with Nextclade assignments at later steps.

Edit the workflow config to provide the correct parameter.

Workflow:
1. Fetch with NCBI Datasets (https://www.ncbi.nlm.nih.gov/datasets/)
    - requires `ncbi_taxon_id` config
    - Directly returns NDJSON without custom parsing
    - Fastest option for large datasets (e.g. SARS-CoV-2)
    - Only returns metadata fields that are available through NCBI Datasets
    - Only works for viral genomes

2. Fetch from Entrez (https://www.ncbi.nlm.nih.gov/books/NBK25501/)
    - requires `entrez_search_term` config (now `entrez_query`)
    - Returns all available data via a GenBank file
    - Requires a custom script to parse the necessary fields from the GenBank file

"""

# Fetch from NCBI Datasets.

rule fetch_ncbi_dataset_package:
    params:
        ncbi_taxon_id=config["ncbi_taxon_id"],
    output:
        dataset_package=temp("data/raw/ncbi/ncbi_dataset.zip"),
    # Allow retries in case of network errors
    retries: 5
    log:
        "logs/fetch_ncbi_dataset_package.txt"
    shell:
        r"""
        datasets download virus genome taxon {params.ncbi_taxon_id:q} \
            --no-progressbar \
            --filename {output.dataset_package:q} \
            > {log:q} 2>&1
        """

rule dump_ncbi_dataset_report:
    """
    This rule is not part of the default workflow.
    It is intended to be used as a specific target for users to be able
    to inspect and explore the full raw metadata from NCBI Datasets.
    """
    input:
        dataset_package="data/raw/ncbi/ncbi_dataset.zip",
    output:
        ncbi_dataset_tsv="data/raw/ncbi/ncbi_dataset_report_raw.tsv",
    log:
        "logs/dump_ncbi_dataset_report.txt"
    shell:
        r"""
        dataformat tsv virus-genome \
            --package {input.dataset_package:q} \
            > {output.ncbi_dataset_tsv:q} 2> {log:q}
        """

rule extract_ncbi_dataset_sequences:
    input:
        dataset_package="data/raw/ncbi/ncbi_dataset.zip",
    output:
        ncbi_dataset_sequences=temp("data/raw/ncbi/ncbi_dataset_sequences.fasta"),
    log:
        "logs/extract_ncbi_dataset_sequences.txt"
    shell:
        r"""
        unzip -jp {input.dataset_package} \
            ncbi_dataset/data/genomic.fna \
            > {output.ncbi_dataset_sequences:q} 2> {log:q}
        """

rule format_ncbi_dataset_report:
    input:
        dataset_package="data/raw/ncbi/ncbi_dataset.zip",
    output:
        ncbi_dataset_tsv=temp("data/raw/ncbi/ncbi_dataset_report.tsv"),
    params:
        ncbi_datasets_fields=",".join(config["ncbi_datasets_fields"]),
    log:
        "logs/format_ncbi_dataset_report.txt"
    shell:
        r"""
        (
            dataformat tsv virus-genome \
                --package {input.dataset_package:q} \
                --fields {params.ncbi_datasets_fields:q} \
                --elide-header \
                | csvtk fix-quotes -Ht \
                | csvtk add-header -t -n {params.ncbi_datasets_fields:q} \
                | csvtk rename -t -f accession -n accession_version \
                | csvtk -t mutate -f accession_version -n accession -p "^(.+?)\." --at 1
        ) > {output.ncbi_dataset_tsv:q} 2> {log:q}
        """

# Technically you can bypass this step and directly provide FASTA and TSV files
# as input files for the curate pipeline.
# We do the formatting here to have a uniform NDJSON file format for the raw
# data that we host on data.nextstrain.org
rule format_ncbi_datasets_ndjson:
    input:
        ncbi_dataset_sequences="data/raw/ncbi/ncbi_dataset_sequences.fasta",
        #ncbi_dataset_tsv="data/ncbi_dataset_report_with_strain.tsv",
        ncbi_dataset_tsv="data/raw/ncbi/ncbi_dataset_report.tsv",

    output:
        ndjson="data/raw/ncbi/ncbi_records.ndjson",
    log:
        "logs/format_ncbi_datasets_ndjson.txt",
    shell:
        r"""
        augur curate passthru \
            --metadata {input.ncbi_dataset_tsv} \
            --fasta {input.ncbi_dataset_sequences} \
            --seq-id-column accession_version \
            --seq-field sequence \
            --unmatched-reporting warn \
            --duplicate-reporting warn \
            > {output.ndjson:q} 2> {log:q}
        """

rule dev_subsample_ncbi_records:
    """
    Development-only helper target.
    This rule is only pulled into the DAG when `config["dev"]` is true and
    curate consumes the dev sample.
    """
    input:
        ndjson="data/raw/ncbi/ncbi_records.ndjson"
    output:
        temp("data/dev/ncbi_records.dev_sample.ndjson")
    params:
        n=config["dev_n_ingest"],
        ref=config["reference_accession"],
    shell:
        r"""
        (
          grep '"accession"[[:space:]]*:[[:space:]]*"{params.ref}"' {input.ndjson} || true
          head -n {params.n} {input.ndjson}
        ) | awk '!seen[$0]++' > {output}
        """


rule select_active_ncbi_records:
    """
    Select the NDJSON consumed by downstream curation.
    In development mode this is the subsampled dev NDJSON; otherwise it is the
    full raw NCBI NDJSON.
    """
    input:
        ndjson=(
            "data/dev/ncbi_records.dev_sample.ndjson"
            if config["dev"]
            else "data/raw/ncbi/ncbi_records.ndjson"
        )
    output:
        ndjson="data/active/ncbi_records.ndjson"
    shell:
        """
        cp {input.ndjson:q} {output.ndjson:q}
        """


rule write_active_ncbi_accessions:
    """
    Record the active NCBI accessions used for development-mode GenBank fetches.
    """
    input:
        ndjson="data/active/ncbi_records.ndjson"
    output:
        temp("data/dev/active_ncbi_accessions.txt")
    shell:
        r"""
        jq -r '.accession // empty' {input.ndjson:q} \
            | awk 'NF && !seen[$0]++' > {output:q}
        """

# Fetch from Entrez.

# overrides.smk
rule fetch_genbank:
    """Fetch GenBank records either from the full Entrez query or the active dev accession set."""
    input:
        accessions=(
            "data/dev/active_ncbi_accessions.txt"
            if config["dev"]
            else []
        )
    params:
        term=config["entrez_query"],
        accessions_arg=(
            "--accessions data/dev/active_ncbi_accessions.txt"
            if config["dev"]
            else ""
        ),
    output:
        genbank="data/raw/entrez/genbank_records.gb"
    shell:
        r"""
        python scripts/fetch-genbank.py \
            --term {params.term:q} \
            --output {output.genbank:q} \
            {params.accessions_arg}
        """

rule add_extra_genomes:
    """
    This step shouldn't be necessary but the NCBI reference genome, NC_003977,
    is not returned via the ENTREZ query. Even if we change the reference it's good
    to add this genome.
    """
    input:
        entrez = "data/raw/entrez/genbank_records.gb",
        ref = config['reference_genbank'],
    output:
        genbank = temp("data/raw/entrez/genbank_records_with_reference.gb"),
    shell:
        """
        cat {input.ref:q} {input.entrez:q} > {output.genbank:q}
        """

rule parse_genbank:
    """Parse fetched GenBank records into NDJSON for metadata-only curation."""
    input:
        genbank = "data/raw/entrez/genbank_records_with_reference.gb",
    output:
        ndjson = "data/raw/entrez/genbank_records.ndjson"
    shell:
        """
        scripts/parse-genbank.py --input {input.genbank} --output {output.ndjson}
        """
