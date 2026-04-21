# nextstrain.org/hbv/ingest

This is the ingest pipeline for Hepatitis B (HBV) virus sequences

> NOTE: This ingest pipeline is in development and the inferred metadata (especially, but not limited to, "clade_nextclade") should not be used for scientific results.

## Software requirements

Follow the [standard installation instructions](https://docs.nextstrain.org/en/latest/install.html) for Nextstrain's suite of software tools.

## Usage

> NOTE: These command examples assume you are within the `ingest` directory.

```sh
snakemake --cores 4
```

> NOTE: The initial data download step can take a very long time, in some runs up to about an hour.

This produces a number of intermediate files in `data/` as well as three files in `results/` for downstream analysis:

- `results/metadata.tsv`
- `results/sequences.fasta`
- `results/aligned.fasta`

## Steps involved

#### GenBank data as inputs

GenBank sequences and metadata are fetched via a NCBI Entrez query.
As of mid 2024 there are around ~11.5k genomes and the full GenBank file is ~150Mb.

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

There is a jupyter notebook exploring the process behind this - see `../notebooks/alignment-qc.ipynb`

#### Accuracy of Nextclade inference

Nextclade v3 is used to align all genomes and assign genotype based on a guide tree we have created.

Preliminary stats can be seen in `ingest/data/qc/metadata_summary.txt` after an ingest build has completed.

### Development mode

Set `dev: true` in `defaults/config.yaml` to generate a small development-only
sample at `data/dev/ncbi_records.dev_sample.ndjson`.

Downstream curation always reads `data/active/ncbi_records.ndjson`:

- with `dev: true`, the active NDJSON is built from the development sample
- with `dev: true`, the Entrez/GenBank fetch is restricted to accessions present
  in the active NCBI development sample
- with `dev: false`, the active NDJSON is built from the full raw NCBI NDJSON in
  `data/raw/ncbi/ncbi_records.ndjson`

## Configuration

Configuration parameters are in `defaults/config.yaml`. These may be overridden by using Snakemake's `--configfile` or `--config` options.

### Environment Variables

None currently required

## `ingest/vendored`

This repository uses [`git subrepo`](https://github.com/ingydotnet/git-subrepo) to manage copies of ingest scripts in [ingest/vendored](./vendored), from [nextstrain/ingest](https://github.com/nextstrain/ingest).

See [vendored/README.md](vendored/README.md#vendoring) for instructions on how to update
the vendored scripts.
