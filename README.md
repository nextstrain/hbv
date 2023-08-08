# Nextstrain Hepatitis B builds

Currently a WIP

Based on Katie Kistler's work in [blab/adaptive-evolution](https://github.com/blab/adaptive-evolution)

This pipeline consists of a single Snakemake workflow with various targets.
The general approach is:
* The data ingest part of the pipeline:
  * Data is downloaded from NCBI via Entrez
  * Genomes are re-cut so that the origin matches the reference origin
  * Nextclade is used to both align genomes and infer the HBV genotype
* A general phylogenetic pipeline with various targets:
  * `all` builds a tree of the entire ~11k dataset
  * `dev` builds a ~500 tip tree
  * `nextclade-tree` is used to update the nextclade dataset (see below)


## Ingest

You can run the ingest part of the snakemake pipeline via:

```
snakemake --cores 4 -npf ingest/results/aligned.fasta ingest/results/sequences.fasta ingest/results/metadata.tsv
```

Ingest will fetch genomes from NCBI's Entrez API and write to `./ingest/data/genbank.gb`.
As of mid 2023 there are around ~11k genomes and the full GenBank file is ~150Mb.

### `ingest/vendored`

This repository uses [`git subrepo`](https://github.com/ingydotnet/git-subrepo) to manage copies of ingest scripts in `ingest/vendored`, from [nextstrain/ingest](https://github.com/nextstrain/ingest). To pull new changes from the central ingest repository, run:

```sh
git subrepo pull ingest/vendored
```

Changes should not be pushed using `git subrepo push`.

1. For pathogen-specific changes, make them in this repository via a pull request.
2. For pathogen-agnostic changes, make them on [nextstrain/ingest](https://github.com/nextstrain/ingest) via pull request there, then use `git subrepo pull` to add those changes to this repository.

## Phylo

* For a small ~500-tip dev tree: `snakemake --cores 4 -pf auspice/hbv_dev.json`

* For the entire human-HBV tree, with 11k tips (takes ~15min on a 4-core M1 machine): `snakemake --cores 4 -pf auspice/hbv_all.json`

* For genotype builds, each of ~500 tips: `snakemake --cores 4 -pf auspice/hbv_A.json`. Currently genotypes A-D are supported.

(If you haven't run the ingest section of the pipeline then any phylo build will run it.)

## Accuracy of Nextclade inference

This is reported within `./ingest/results/metadata.summary.txt`.
An example summary using genomes fetched on 2023-06-27 and `nextclade_datasets/references/JN182318/versions/2023-06-26`:

```
Overall Nextclade QC status
Failed align   *                      3.7% (n=425/11358)
bad            **                    11.9% (n=1346/11358)
mediocre       **                    10.7% (n=1212/11358)
good           ***************       73.7% (n=8375/11358)
```

```
Accuracy of inferred genotypes (nextclade) vs metadata-assigned genotypes (GenBank)
A              ********************  99.3% (n=900/906)
B              ********************  98.6% (n=1235/1252)
C              ********************  97.7% (n=2895/2962)
D              *******************   95.2% (n=1701/1786)
E              ******************** 100.0% (n=358/358)
F              ******************** 100.0% (n=292/292)
G              ****************      79.1% (n=34/43)
H              ******************** 100.0% (n=25/25)
I              ********************  99.2% (n=131/132)
recombinant                           0.0% (n=0/372)
```

## Updating Nextclade dataset


Currently a HBV dataset exists using reference JN182318 and a ~2000 tip tree attempting to cover observed human HBV diversity and genotypes.

Create a new version of the dataset by copying a previous datestamped directory in `nextclade_datasets/references/JN182318/version/`

Update the tree and potentially the example sequences:

```bash
snakemake --cores 2 auspice/hbv_nextclade-tree.json
cp auspice/hbv_nextclade-tree.json nextclade_datasets/references/JN182318/versions/YYYY-MM-DD/tree.json

snakemake --cores 1 results/nextclade-sequences/filtered.fasta
cp results/nextclade-sequences/filtered.fasta nextclade_datasets/references/JN182318/versions/YYYY-MM-DD/sequences.fasta
```

> NOTE: Many of the files - especially `qc.json` - still need to be optimised for HBV.

## TODO

* Sanitise metadata
* Check occurrences of a clade of (genbank-labelled) genotype G within the I clade
* Optimise nextclade parameters
* The mutation counts are incredible & need examining
