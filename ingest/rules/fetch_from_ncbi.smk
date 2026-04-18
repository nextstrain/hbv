"""
This part of the workflow handles fetching sequences and metadata from NCBI.

REQUIRED INPUTS:

    None

OUTPUTS:

    ndjson = data/ncbi.ndjson
    ndjson_genbank = "data/genbank.ndjson" # for metadata only



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



###########################################################################
####################### 1. Fetch from NCBI Datasets #######################
###########################################################################



rule fetch_ncbi_dataset_package:
    params:
        ncbi_taxon_id=config["ncbi_taxon_id"],
    output:
        dataset_package=temp("data/ncbi_dataset.zip"),
    # Allow retries in case of network errors
    retries: 5
    log:
        "logs/fetch_ncbi_dataset_package.txt"
    benchmark:
        "benchmarks/fetch_ncbi_dataset_package.txt"
    shell:
        r"""
        exec &> >(tee {log:q})

        datasets download virus genome taxon {params.ncbi_taxon_id:q} \
            --no-progressbar \
            --filename {output.dataset_package}
        """

# Note: This rule is not part of the default workflow!
# It is intended to be used as a specific target for users to be able
# to inspect and explore the full raw metadata from NCBI Datasets.
rule dump_ncbi_dataset_report:
    input:
        dataset_package="data/ncbi_dataset.zip",
    output:
        ncbi_dataset_tsv="data/ncbi_dataset_report_raw.tsv",
    log:
        "logs/dump_ncbi_dataset_report.txt"
    shell:
        r"""
        exec &> >(tee {log:q})

        dataformat tsv virus-genome \
            --package {input.dataset_package} > {output.ncbi_dataset_tsv}
        """


rule extract_ncbi_dataset_sequences:
    input:
        dataset_package="data/ncbi_dataset.zip", 
    output:
        ncbi_dataset_sequences=temp("data/ncbi_dataset_sequences.fasta"),
    log:
        "logs/extract_ncbi_dataset_sequences.txt"
    benchmark:
        "benchmarks/extract_ncbi_dataset_sequences.txt"
    shell:
        r"""
        exec &> >(tee {log:q})

        unzip -jp {input.dataset_package} \
            ncbi_dataset/data/genomic.fna > {output.ncbi_dataset_sequences}
        """

rule format_ncbi_dataset_report:
    input:
        dataset_package="data/ncbi_dataset.zip",
    output:
        ncbi_dataset_tsv=temp("data/ncbi_dataset_report.tsv"),
    params:
        ncbi_datasets_fields=",".join(config["ncbi_datasets_fields"]),
    log:
        "logs/format_ncbi_dataset_report.txt"
    benchmark:
        "benchmarks/format_ncbi_dataset_report.txt"
    shell:
        r"""
        exec &> >(tee {log:q})

        dataformat tsv virus-genome \
            --package {input.dataset_package} \
            --fields {params.ncbi_datasets_fields:q} \
            --elide-header \
            | csvtk fix-quotes -Ht \
            | csvtk add-header -t -n {params.ncbi_datasets_fields:q} \
            | csvtk rename -t -f accession -n accession_version \
            | csvtk -t mutate -f accession_version -n accession -p "^(.+?)\." --at 1 \
            > {output.ncbi_dataset_tsv}
        """


# Technically you can bypass this step and directly provide FASTA and TSV files
# as input files for the curate pipeline.
# We do the formatting here to have a uniform NDJSON file format for the raw
# data that we host on data.nextstrain.org
rule format_ncbi_datasets_ndjson:
    input:
        ncbi_dataset_sequences="data/ncbi_dataset_sequences.fasta",
        #ncbi_dataset_tsv="data/ncbi_dataset_report_with_strain.tsv",
        ncbi_dataset_tsv="data/ncbi_dataset_report.tsv",

    output:
        ndjson="data/ncbi.ndjson",
    log:
        "logs/format_ncbi_datasets_ndjson.txt",
    benchmark:
        "benchmarks/format_ncbi_datasets_ndjson.txt"
    shell:
        r"""
        exec &> >(tee {log:q})

        augur curate passthru \
            --metadata {input.ncbi_dataset_tsv} \
            --fasta {input.ncbi_dataset_sequences} \
            --seq-id-column accession_version \
            --seq-field sequence \
            --unmatched-reporting warn \
            --duplicate-reporting warn \
            2> {log} > {output.ndjson}
        """


rule ncbi_active:
    input:
        ndjson="data/ncbi.ndjson"
    output:
        ACTIVE_NDJSON
    params:
        dev=config.get("dev", False),
        n=config.get("dev_n_ingest", 100),
        ref=config["reference_accession"],
    shell:
        r"""
        # mkdir -p data

        if [ "{params.dev}" = "true" ] || [ "{params.dev}" = "True" ] || [ "{params.dev}" = "1" ]; then
          (
            grep '"accession"[[:space:]]*:[[:space:]]*"{params.ref}"' {input.ndjson} || true
            head -n {params.n} {input.ndjson}
          ) | awk '!seen[$0]++' > {output}
        else
          cp {input.ndjson} {output}
        fi
        """




###########################################################################
########################## 2. Fetch from Entrez ###########################
###########################################################################


# overrides.smk
rule fetch_genbank:
    params:
        term=config["entrez_query"]
    output:
        genbank="data/entrez/genbank.gb"
    shell:
        r"""
        python - {params.term:q} {output.genbank:q} <<'PY'
import json, sys, time, random
from http.client import IncompleteRead
from urllib.error import HTTPError, URLError
from Bio import SeqIO, Entrez

Entrez.email = "hello@nextstrain.org"
BATCH_SIZE = 1000

def get_esearch_history(term):
    handle = Entrez.esearch(
        db="nucleotide",
        term=term,
        retmode="json",
        usehistory="y",
        retmax=0,
    )
    esearch_result = json.loads(handle.read())["esearchresult"]
    print(f"Search term {{term!r}} returned {{esearch_result['count']}} IDs.")
    return {{
        "count": int(esearch_result["count"]),
        "query_key": esearch_result["querykey"],
        "web_env": esearch_result["webenv"],
    }}

def fetch_batch(query_key, web_env, start, tries=8):
    for attempt in range(tries):
        try:
            handle = Entrez.efetch(
                db="nucleotide",
                query_key=query_key,
                webenv=web_env,
                retstart=start,
                retmax=BATCH_SIZE,
                rettype="gb",
                retmode="text",
            )
            return handle.read()
        except (IncompleteRead, HTTPError, URLError, OSError):
            time.sleep(min(60, (2 ** attempt) + random.random()))
    raise RuntimeError(f"efetch failed after {{tries}} tries at retstart={{start}}")

def main(term, out_path):
    h = get_esearch_history(term)
    count, query_key, web_env = h["count"], h["query_key"], h["web_env"]

    print(f"Fetching GenBank records in batches of n={{BATCH_SIZE}}")
    with open(out_path, "w") as output_handle:
        written = 0
        for start in range(0, count, BATCH_SIZE):
            records = fetch_batch(query_key, web_env, start)
            output_handle.write(records)
            output_handle.flush()
            written += records.count("\nLOCUS")
            print(f"[batch] total_written={{written}}")
            time.sleep(0.4)

if __name__ == "__main__":
    term = sys.argv[1]
    out_path = sys.argv[2]
    main(term, out_path)
PY
        """


rule add_extra_genomes:
    """
    This step shouldn't be necessary but the NCBI reference genome, NC_003977,
    is not returned via the ENTREZ query. Even if we change the reference it's good
    to add this genome.
    """
    input:
        entrez = "data/entrez/genbank.gb",
        ref = config['reference_genbank'],
    output:
        genbank = "data/entrez/genbank.with-reference.gb",
    shell:
        """
        cat {input.ref:q} {input.entrez:q} > {output.genbank:q}
        """


rule parse_genbank:
    input:
        genbank = "data/entrez/genbank.with-reference.gb",
    output:
        ndjson = "data/genbank.ndjson"
    shell:
        """
        scripts/parse-genbank.py --input {input.genbank} --output {output.ndjson}
        """

