# nextstrain.org/hbv/ingest

This is the ingest pipeline for Hepatitis B (HBV) virus sequences

> NOTE: This ingest pipeline is in development and the inferred metadata (especially, but not limited to, "clade_nextclade") should not be used for scientific results.

> NOTE: Ingest dev mode is incompatible with other directories. First run ingest properly before rerunning phylogenetics or exploratory.

## Software requirements

Follow the [standard installation instructions](https://docs.nextstrain.org/en/latest/install.html) for Nextstrain's suite of software tools.

## Usage

> NOTE: These command examples assume you are within the `ingest` directory.

```sh
snakemake --cores 4
```

> NOTE: The initial data download step can take a very long time, in some runs up to about an hour. For a much faster Entrez fetch, set `entrez_complete_genomes_only: true`; this currently restricts the GenBank branch to about 11.7k complete-genome records. Development mode already enforces this filter.

This produces a number of intermediate files in `data/` as well as three files in `results/` for downstream analysis:

- `results/metadata.tsv`
- `results/sequences.fasta`
- `results/aligned.fasta`

## Steps involved

#### GenBank data as inputs

GenBank sequences and metadata are fetched via a NCBI Entrez query.
As of April 21, 2026, the default `Hepatitis B virus[Organism]` query returns about 138k nucleotide records, while the optional `complete genome` filter reduces this to about 11.7k records.

> NOTE: Genotype and subgenotype annotations are extracted manually from free-text GenBank annotation. The current parsing script is likely still too harsh for some common note formats and should be reviewed again before relying on these fields.

Raw NCBI- and Entrez-derived inputs are organized under:

- `data/raw/ncbi/`
- `data/raw/entrez/`

The NDJSON actively consumed by curation is kept separately under:

- `data/active/`

Curated and post-curation outputs are organized under:

- `data/curated/ncbi/`
- `data/curated/genbank/`
- `data/merged/`
- `data/circularised/`
- `data/nextclade/`
- `data/qc/`

#### Genomes rotated to use a consistent origin

There is a jupyter notebook exploring the process behind this - see `../legacy/notebooks/alignment-qc.ipynb` (using old reference!)

#### Accuracy of Nextclade inference

Nextclade v3 is used to align all genomes and assign genotype based on a guide tree we have created.
The local dataset used for this step is configured via `nextclade_dataset` in `defaults/config.yaml`.
The translated CDS FASTAs in `data/nextclade/translations/` are also used downstream by the phylogenetics workflow.

Preliminary stats can be seen in `ingest/data/qc/metadata_summary.txt` after an ingest build has completed.

### Development mode

Set `dev: true` in `defaults/config.yaml` to generate a small development-only
cohort controlled by `dev_n_ingest`.

Downstream curation always reads `data/active/ncbi_records.ndjson`:

- with `dev: true`, ingest fetches the first `dev_n_ingest` Entrez/GenBank
  records from the query result and always restricts those queries to
  `complete genome` records
- with `dev: true`, the active NDJSON is then filtered to those accessions
- with `dev: false`, the active NDJSON is copied from the full raw NCBI NDJSON
  unless `entrez_complete_genomes_only: true`, in which case it is filtered to
  the fetched Entrez accession set

Outside development mode, the same Entrez complete-genome restriction can be
enabled independently with `entrez_complete_genomes_only: true`.

## Configuration

Configuration parameters are in `defaults/config.yaml`. These may be overridden by using Snakemake's `--configfile` or `--config` options.

### Environment Variables

None currently required

## `ingest/vendored`

This repository uses [`git subrepo`](https://github.com/ingydotnet/git-subrepo) to manage copies of ingest scripts in [ingest/vendored](./vendored), from [nextstrain/ingest](https://github.com/nextstrain/ingest).

See [vendored/README.md](vendored/README.md#vendoring) for instructions on how to update
the vendored scripts.
