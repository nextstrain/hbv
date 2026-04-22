"""
Rules for fetching the raw HBV record sets used by ingest.

Required external inputs:
- no upstream workflow outputs; this module starts from configured NCBI sources
- local defaults such as the reference GenBank file used by `add_extra_genomes`

Key outputs:
- `data/raw/ncbi/ncbi_records.ndjson`
- `data/raw/entrez/genbank_records.ndjson`
- `data/active/ncbi_records.ndjson`
- `data/active/active_entrez_accessions.txt` when active-set subsetting is used

The workflow uses both NCBI Datasets and Entrez. NCBI Datasets provides the
main sequence set used downstream, while Entrez provides GenBank records used
for additional HBV genotype and subgenotype annotations.

Workflow:
1. Fetch with NCBI Datasets (https://www.ncbi.nlm.nih.gov/datasets/)
    - requires `ncbi_taxon_id` config
    - directly returns NDJSON without custom parsing
    - fastest option for large datasets
    - only returns metadata fields available through NCBI Datasets
    - only works for viral genomes

2. Fetch from Entrez (https://www.ncbi.nlm.nih.gov/books/NBK25501/)
    - requires `entrez_query` config
    - returns full GenBank records
    - requires custom parsing to extract the fields used downstream
"""

# Fetch from NCBI Datasets.

rule fetch_ncbi_dataset_package:
    """Download the raw NCBI Datasets package used to derive the ingest FASTA and metadata report."""
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
    """Extract the genomic FASTA from the downloaded NCBI Datasets archive."""
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
    """Select and normalize the NCBI Datasets metadata fields used to build raw ingest NDJSON."""
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
    """Combine the formatted NCBI metadata table and FASTA into raw NDJSON records."""
    input:
        ncbi_dataset_sequences="data/raw/ncbi/ncbi_dataset_sequences.fasta",
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

rule select_active_ncbi_records:
    """
    Select the NDJSON consumed by downstream curation.
    In development mode, or when complete-genome Entrez filtering is enabled,
    this is the subset of raw NCBI records whose accessions were first selected
    from the Entrez GenBank sample. Otherwise it is the full raw NCBI NDJSON.
    """
    input:
        ndjson="data/raw/ncbi/ncbi_records.ndjson",
        accessions=(
            "data/active/active_entrez_accessions.txt"
            if SUBSET_ACTIVE_NCBI
            else []
        )
    output:
        ndjson="data/active/ncbi_records.ndjson"
    params:
        subset_active=str(SUBSET_ACTIVE_NCBI).lower(),
        verbose=str(config.get("verbose", False)).lower(),
    shell:
        r"""
        if [ "{params.subset_active}" = "true" ] || [ "{params.subset_active}" = "True" ] || [ "{params.subset_active}" = "1" ]; then
          HBV_VERBOSE={params.verbose} python scripts/subset_ndjson_by_accessions.py \
            --input {input.ndjson:q} \
            --accessions {input.accessions:q} \
            --output {output.ndjson:q}
        else
          cp {input.ndjson:q} {output.ndjson:q}
        fi
        """


rule write_active_entrez_accessions:
    """
    Record the fetched Entrez accessions used to define the active NCBI subset.
    """
    input:
        genbank="data/raw/entrez/genbank_records.gb"
    output:
        temp("data/active/active_entrez_accessions.txt")
    shell:
        r"""
        grep '^ACCESSION' {input.genbank:q} \
            | awk '{{print $2}}' \
            | awk 'NF && !seen[$0]++' > {output:q}
        """

# Fetch from Entrez.

rule fetch_genbank:
    """Fetch GenBank records either from the full Entrez query or a small development sample."""
    params:
        term=config["entrez_query"],
        verbose=str(config.get("verbose", False)).lower(),
        dev_args=(
            f"--limit {config['dev_n_ingest']} "
            f"--complete-genomes"
            if DEV_MODE
            else ""
        ),
        full_query_args=(
            "--complete-genomes"
            if config.get("entrez_complete_genomes_only", False) and not DEV_MODE
            else ""
        ),
    output:
        genbank="data/raw/entrez/genbank_records.gb"
    shell:
        r"""
        HBV_VERBOSE={params.verbose} python scripts/fetch-genbank.py \
            --term {params.term:q} \
            --output {output.genbank:q} \
            {params.full_query_args} \
            {params.dev_args}
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
    params:
        verbose=str(config.get("verbose", False)).lower(),
    shell:
        """
        HBV_VERBOSE={params.verbose} scripts/parse-genbank.py --input {input.genbank} --output {output.ndjson}
        """
