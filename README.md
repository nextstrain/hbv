# Nextstrain HBV (Hepatitis-B) builds

> This repo is currently experimental and the results are not scientifically accurate
(most significantly, recombination has not been accounted for).
Latest available trees: [genotype A](https://nextstrain.org/staging/hbv/A), [genotype B](https://nextstrain.org/staging/hbv/B), [genotype C](https://nextstrain.org/staging/hbv/C), [genotype D](https://nextstrain.org/staging/hbv/D), [all 11k genomes](https://nextstrain.org/staging/hbv/all)


Based on Katie Kistler's work in [blab/adaptive-evolution](https://github.com/blab/adaptive-evolution)


This repository contains two workflows for the analysis of HBV virus data:

- [`ingest/`](./ingest) - Download data from GenBank, curate metadata, rotate genomes, align and infer genotypes using Nextclade
- [`phylogenetic/`](./phylogenetic) - Filter sequences, construct phylogeny and export for visualization

Each folder contains a README.md with more information.