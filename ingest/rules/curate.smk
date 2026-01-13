"""
This part of the workflow handles the curation of data from NCBI

REQUIRED INPUTS:

    ndjson      = data/ncbi.ndjson

OUTPUTS:

    sequences = "data/circularised.fasta",
    metadata = "data/circularised.tsv",


"""




def format_field_map(field_map: dict[str, str]) -> str:
    """
    Format dict to `"key1"="value1" "key2"="value2"...` for use in shell commands.
    """
    return " ".join([f'"{key}"="{value}"' for key, value in field_map.items()])







# This curate pipeline is based on existing pipelines for pathogen repos using NCBI data.
# You may want to add and/or remove steps from the pipeline for custom metadata
# curation for your pathogen. Note that the curate pipeline is streaming NDJSON
# records between scripts, so any custom scripts added to the pipeline should expect
# the input as NDJSON records from stdin and output NDJSON records to stdout.
# The final step of the pipeline should convert the NDJSON records to two
# separate files: a metadata TSV and a sequences FASTA.
rule curate:
    input:
        sequences_ndjson="data/ncbi.subset3.ndjson",
    output:
        metadata = "data/curated-metadata.tsv",
        sequences = "data/curated-sequences.fasta",

    log:
        "logs/curate.txt",
    benchmark:
        "benchmarks/curate.txt"
    params:
        field_map=format_field_map(config["curate"]["field_map"]),
        strain_regex=config["curate"]["strain_regex"],
        strain_backup_fields=config["curate"]["strain_backup_fields"],
        date_fields=config["curate"]["date_fields"],
        expected_date_formats=config["curate"]["expected_date_formats"],
        genbank_location_field=config["curate"]["genbank_location_field"],
        articles=config["curate"]["titlecase"]["articles"],
        abbreviations=config["curate"]["titlecase"]["abbreviations"],
        titlecase_fields=config["curate"]["titlecase"]["fields"],
        authors_field=config["curate"]["authors_field"],
        authors_default_value=config["curate"]["authors_default_value"],
        abbr_authors_field=config["curate"]["abbr_authors_field"],
        id_field=config["curate"]["output_id_field"],
        sequence_field=config["curate"]["output_sequence_field"],
    shell:
        r"""
        exec &> >(tee {log:q})

        cat {input.sequences_ndjson} \
            | augur curate rename \
                --field-map {params.field_map} \
            | augur curate normalize-strings \
            | augur curate transform-strain-name \
                --strain-regex {params.strain_regex} \
                --backup-fields {params.strain_backup_fields} \
            | augur curate format-dates \
                --date-fields {params.date_fields} \
                --expected-date-formats {params.expected_date_formats} \
            | augur curate parse-genbank-location \
                --location-field {params.genbank_location_field} \
            | augur curate titlecase \
                --titlecase-fields {params.titlecase_fields} \
                --articles {params.articles} \
                --abbreviations {params.abbreviations} \
            | augur curate abbreviate-authors \
                --authors-field {params.authors_field} \
                --default-value {params.authors_default_value} \
                --abbr-authors-field {params.abbr_authors_field} \
            | scripts/fix_country_field.py \
            | vendored/apply-geolocation-rules --geolocation-rules defaults/geoLocationRules.tsv \
            | scripts/add-year.py \
            | jq -c '.name = (.name // .accession // .accession_version // "")' \
            | augur curate passthru \
                --output-metadata {output.metadata} \
                --output-fasta {output.sequences} \
                --output-id-field {params.id_field} \
                --output-seq-field {params.sequence_field}
        """

rule add_metadata_columns:
    """Add columns to metadata
    Notable columns:
    - [NEW] url: URL linking to the NCBI GenBank record ('https://www.ncbi.nlm.nih.gov/nuccore/*').
    """
    input:
        metadata = "data/curated-metadata.tsv"
    output:
        metadata = temp("data/all_metadata_added.tsv")
    log:
        "logs/add_metadata_columns.txt"
    params:
        accession=config['curate']['genbank_accession']
    shell:
        r"""
        exec &> >(tee {log:q})

        csvtk mutate2 -t \
          -n url \
          -e '"https://www.ncbi.nlm.nih.gov/nuccore/" + ${params.accession}' \
          {input.metadata} \
        > {output.metadata}
        """

rule subset_metadata:
    input:
        metadata="data/all_metadata_added.tsv",
    output:
        subset_metadata="data/subset_metadata.tsv",
    log:
        "logs/subset_metadata.txt"
    params:
        metadata_fields=",".join(config["curate"]["metadata_columns"]),
    shell:
        r"""
        exec &> >(tee {log:q})

        csvtk cut -t -f {params.metadata_fields} \
            {input.metadata} > {output.subset_metadata}
        """


rule recircularise:
    input:
        metadata = "data/subset_metadata.tsv",
        sequences = "data/curated-sequences.fasta",
    output:
        sequences = "data/circularised.fasta",
        metadata = "data/circularised.tsv",
    params:
        reference = config['reference_accession']
    shell:
        """
        scripts/re-circularise.py \
            --seqs-in {input.sequences} --meta-in {input.metadata} \
            --seqs-out {output.sequences} --meta-out {output.metadata} \
            --reference {params.reference}
        """


    
rule align_unrotated:
    """
    Align all genomes before rotation for parsing by our notebook.
    This rule must be called explicitly, it is not part of the DAG to produce the outputs of `rule all`
    """
    input:
        sequences = "data/curated-sequences.fasta",
    output:
        alignment = "data/curated-sequences.aligned.fasta",
    params:
        dataset = config['nextclade_dataset'],
    threads: 4
    shell:
        """
        nextclade run \
            -j {threads} --silent --replace-unknown \
            --input-dataset {params.dataset} \
            --output-fasta {output.alignment} \
            {input.sequences}
        """