# Nextstrain HBV (Hepatitis-B) builds

> This repo is currently experimental and the results may not be scientifically accurate

Based on work of James Hadfield, and in turn on Katie Kistler's work in [blab/adaptive-evolution](https://github.com/blab/adaptive-evolution)

This repository contains two workflows for the analysis of HBV virus data:

- [`ingest/`](./ingest) - Download data from GenBank, curate metadata, rotate genomes, align and infer genotypes using Nextclade
- [`phylogenetics/`](./phylogenetics) - Filter sequences, construct phylogeny and export for visualization. Generate Nextclade dataset

Additionally available are:

- [`exploratory/`](./exploratory) - Exploration of HBV phylogenetic changes between trees of different genomic segments - to be removed
- [`legacy/`](./legacy) - Old functionality within the phylogenetic workflow by James Hadfield not yet fully incoporated in the current workflow

Each folder contains a `README.md` with more information.

For more information about HBV, recombination, and workflow options, see the [`report`](https://polybox.ethz.ch/index.php/s/xeAyC7tJqgqndpw) associated with the April 2026 version of this repository.
