# Nextstrain HBV (Hepatitis-B) builds

> This repo is currently experimental and the results may not be scientifically accurate

Based on work of James Hadfield, (and in turn on Katie Kistler's work in [blab/adaptive-evolution](https://github.com/blab/adaptive-evolution))

This repository contains two workflows for the analysis of HBV virus data:

- [`ingest/`](./ingest) - Download data from GenBank, curate metadata, rotate genomes, align and infer genotypes using Nextclade
- [`phylogenetics/`](./phylogenetics) - Filter sequences, construct phylogeny and export for visualization. Generate Nextclade dataset

Additionally available are:

- [`exploratory/`](./exploratory) - Exploration of HBV phylogenetic changes between trees of different genomic segments - to be removed
- [`legacy/`](./legacy) - Old functionality within the phylogenetic workflow by James Hadfield not yet fully incoporated in the current workflow

Each folder contains a `README.md` with more information.

For standard use, create the environment from [`nextstrain.yml`](./nextstrain.yml), run `ingest/` first, and then `phylogenetics/`:

```bash
mamba env create -f nextstrain.yml
conda activate hbv
cd ingest && snakemake --cores 4
cd ../phylogenetics && snakemake --cores 4 -pf
```

> Note that after running ingest in dev mode, the subsampling in phylogenetics will likely  not recover enough sequences! Run ingest in normal mode before proceeding!

> The initial data download in `ingest/` can take a very long time, in some runs up to about an hour.

This produces the main handoff files in `ingest/results/` (`metadata.tsv`, `sequences.fasta`, `aligned.fasta`) and then, by default, the Nextstrain/Auspice outputs in `phylogenetics/auspice_datasets/`. 

Running `phylogenetics/` with `defaults/nextclade/config_nextclade.yaml` instead builds a Nextclade dataset in [`nextclade_datasets/`](./nextclade_datasets). 
The `exploratory/` and `legacy/` directories are supplementary reference/development material rather than the main workflow path.


For more information about HBV, recombination, and workflow options, see the [`report`](https://polybox.ethz.ch/index.php/s/xeAyC7tJqgqndpw) associated with the April 2026 version of this repository.
