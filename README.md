# Nextstrain Hepatitis B builds

Currently a WIP

Based on Katie Kistler's work in [blab/adaptive-evolution](https://github.com/blab/adaptive-evolution)

## Ingest

You can run the ingest part of the snakemake pipeline via:

```
snakemake --cores 4 -npf ingest/results/aligned.fasta ingest/results/sequences.fasta ingest/results/metadata.tsv
```

Ingest will fetch genomes from NCBI's Entrez API and write to `./ingest/data/genbank.gb`.
As of mid 2023 there are around ~11k genomes and the full GenBank file is ~150Mb.

## Phylo

For a small ~800-tip dev tree: `snakemake --cores 4 -pf auspice/hbv_dev.json`
For the entire human-HBV tree, with 11k tips (takes ~15min on a 4-core M1 machine): `snakemake --cores 4 -pf auspice/hbv_all.json`

(If you haven't run the ingest section of the pipeline then any phylo build will run it.)



## Nextclade


Currently an all-human HBV dataset exists, with reference JN182318

Create a new version of the dataset by copying a datestamped directory in `nextclade_datasets/references/JN182318/version`

Update the tree and the (example) sequences:

```bash
snakemake --cores 2 auspice/hbv_nextclade-tree.json
cp auspice/hbv_nextclade-tree.json nextclade_datasets/references/JN182318/versions/YYYY-MM-DD/tree.json

snakemake --cores 1 results/nextclade-sequences/filtered.fasta
cp results/nextclade-sequences/filtered.fasta nextclade_datasets/references/JN182318/versions/YYYY-MM-DD/sequences.fasta
```

Many of the files - especially `qc.json` - need to be optimised for HBV.

Example incantation to call genomes against the above dataset:
```bash
nextclade run \
  --input-dataset nextclade_datasets/references/JN182318/versions/2023-06-26 \
  --output-all results/nextclade-test \
  --output-basename test \
  results/nextclade-sequences/filtered.fasta
```
