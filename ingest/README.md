# nextstrain.org/hbv/ingest

This is the ingest pipeline for Hepatitis B (HBV) virus sequences

## Software requirements

Follow the [standard installation instructions](https://docs.nextstrain.org/en/latest/install.html) for Nextstrain's suite of software tools.

## Usage

> NOTE: These command examples assume you are within the `ingest` directory.

```sh
snakemake --cores 4
```

This produces a number of intermediate files in `data/` as well as three files in `results/` for downstream analysis:

- `results/metadata.tsv`
- `results/sequences.fasta`
- `results/aligned.fasta`

## Configuration

Configuration parameters are in `defaults/config.yaml`. These may be overridden by using Snakemake's `--configfile` or `--config` options.

### Environment Variables

None currently required

## Input data

### GenBank data

GenBank sequences and metadata are fetched via a NCBI Entrez query.
As of mid 2024 there are around ~11.5k genomes and the full GenBank file is ~150Mb.


## `ingest/vendored`

This repository uses [`git subrepo`](https://github.com/ingydotnet/git-subrepo) to manage copies of ingest scripts in [ingest/vendored](./vendored), from [nextstrain/ingest](https://github.com/nextstrain/ingest).

See [vendored/README.md](vendored/README.md#vendoring) for instructions on how to update
the vendored scripts.
