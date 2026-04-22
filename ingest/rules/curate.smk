"""
Rules for curating fetched NCBI and GenBank records into the ingest handoff.

Required external inputs:
- `data/active/ncbi_records.ndjson`
- `data/raw/entrez/genbank_records.ndjson`
- local defaults such as `geoLocationRules.tsv` and the subgenotype map

Key outputs:
- `data/circularised/circularised_metadata.tsv`
- `data/circularised/circularised_sequences.fasta`
- `results/sequences.fasta`
"""

def format_field_map(field_map: dict[str, str]) -> str:
    """
    Format dict to `"key1"="value1" "key2"="value2"...` for use in shell commands.
    """
    return " ".join([f'"{key}"="{value}"' for key, value in field_map.items()])

rule curate_genbank_metadata:
    """
    Normalize curated GenBank metadata and apply the HBV subgenotype correction table.
    """
    input:
        ndjson = "data/raw/entrez/genbank_records.ndjson",
        geolocations = "defaults/geoLocationRules.tsv",
        mapping = "defaults/subgenotype_correction.tsv",

    output:
        metadata = "data/curated/genbank/curated_metadata.tsv",
        sequences = temp("data/curated/genbank/curated_sequences.fasta"),
        metadata_pre_mapping = temp(
            "data/curated/genbank/curated_metadata_pre_subgenotype_mapping.tsv"
        ),
    shell:
        """
        cat {input.ndjson} \
            | scripts/fix_country_field.py \
            | vendored/apply-geolocation-rules --geolocation-rules {input.geolocations} \
            | scripts/add-year.py \
            | augur curate passthru \
                --output-seq-field sequence --output-id-field accession \
                --output-metadata {output.metadata_pre_mapping} --output-fasta {output.sequences}

        python "scripts/subgenotype_mapping.py" \
            --metadata-in {output.metadata_pre_mapping} \
            --mapping {input.mapping} \
            --metadata-out {output.metadata}

        """

# This curate pipeline is based on existing pipelines for pathogen repos using NCBI data.
# You may want to add and/or remove steps from the pipeline for custom metadata
# curation for your pathogen. Note that the curate pipeline is streaming NDJSON
# records between scripts, so any custom scripts added to the pipeline should expect
# the input as NDJSON records from stdin and output NDJSON records to stdout.
# The final step of the pipeline should convert the NDJSON records to two
# separate files: a metadata TSV and a sequences FASTA.
rule curate_ncbi:
    """Run the main streaming NCBI curation pipeline and write curated metadata plus sequences."""
    input:
        sequences_ndjson="data/active/ncbi_records.ndjson",
        geolocations = "defaults/geoLocationRules.tsv"

    output:
        metadata = "data/curated/ncbi/curated_metadata.tsv",
        sequences = "data/curated/ncbi/curated_sequences.fasta",

    log:
        "logs/curate.txt",
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
        (
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
                | vendored/apply-geolocation-rules --geolocation-rules {input.geolocations} \
                | scripts/add-year.py \
                | jq -c '.name = (.name // .accession // .accession_version // "")' \
                | augur curate passthru \
                    --output-metadata {output.metadata} \
                    --output-fasta {output.sequences} \
                    --output-id-field {params.id_field} \
                    --output-seq-field {params.sequence_field}
        ) > {log:q} 2>&1
        """

rule add_metadata_columns:
    """Add columns to metadata
    Notable columns:
    - [NEW] url: URL linking to the NCBI GenBank record ('https://www.ncbi.nlm.nih.gov/nuccore/*').
    - "genotype_genbank", "subgenotype_genbank": (self-assigned) strain annotation from Genbank records
    """
    input:
        metadata="data/curated/ncbi/curated_metadata.tsv",
        metadata_genbank="data/curated/genbank/curated_metadata.tsv",
    output:
        metadata=temp("data/merged/metadata_with_genbank_annotations.tsv"),
    log:
        "logs/add_metadata_columns.txt"
    params:
        accession_col="accession",
        verbose=str(config.get("verbose", False)).lower(),
    shell:
        r"""
        HBV_VERBOSE={params.verbose} python scripts/add_genbank_metadata.py \
          --metadata {input.metadata} \
          --metadata-genbank {input.metadata_genbank} \
          --accession-col {params.accession_col} \
          --out {output.metadata} \
          > {log:q} 2>&1
    """

rule subset_metadata:
    """Keep only the metadata columns needed by the recircularisation step."""
    input:
        metadata="data/merged/metadata_with_genbank_annotations.tsv",
    output:
        subset_metadata=temp("data/merged/metadata_for_circularisation.tsv"),
    log:
        "logs/subset_metadata.txt"
    params:
        metadata_fields=",".join(config["curate"]["metadata_columns"]),
    shell:
        r"""
        csvtk cut -t -f {params.metadata_fields} \
            {input.metadata} > {output.subset_metadata} 2> {log:q}
        """

rule recircularise:
    """Rotate sequences to a consistent origin and annotate the corresponding metadata."""
    input:
        metadata = "data/merged/metadata_for_circularisation.tsv",
        sequences = "data/curated/ncbi/curated_sequences.fasta",
        reference_genbank = config["reference_genbank"],
    output:
        sequences = "data/circularised/circularised_sequences.fasta",
        metadata = "data/circularised/circularised_metadata.tsv",
    params:
        reference = config['reference_accession'],
        verbose=str(config.get("verbose", False)).lower(),
    shell:
        """
        HBV_VERBOSE={params.verbose} scripts/re-circularise.py \
            --seqs-in {input.sequences} --meta-in {input.metadata} \
            --seqs-out {output.sequences} --meta-out {output.metadata} \
            --reference {params.reference} \
            --reference-genbank {input.reference_genbank}
        """

rule copy_ingest_sequences:
    """Copy the circularised sequences to the final ingest results directory."""
    input:
        sequences = "data/circularised/circularised_sequences.fasta",
    output:
        sequences = "results/sequences.fasta",
    shell:
        """
        cp {input.sequences} {output.sequences}
        """

rule align_unrotated:
    """
    Align all genomes before rotation for parsing by our notebook.
    This rule must be called explicitly, it is not part of the DAG to produce the outputs of `rule all`
    """
    input:
        sequences = "data/curated/ncbi/curated_sequences.fasta",
    output:
        alignment = "data/qc/unrotated_aligned_sequences.fasta",
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
