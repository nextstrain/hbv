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

This repository uses [`git subrepo`](https://github.com/ingydotnet/git-subrepo) to manage copies of ingest scripts in [`ingest/vendored`](./ingest/vendored), from [nextstrain/ingest](https://github.com/nextstrain/ingest). To pull new changes from the central ingest repository, run:

See [ingest/vendored/README.md](./ingest/vendored/README.md#vendoring) for instructions on how to update
the vendored scripts.

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

# Updating the Nextclade dataset

Nextclade datasets exist for reference `NC_003977` and `JN182318`, however the latter is deprecated and will be removed shortly.
The dataset includes a ~2000 tip tree attempting to cover observed human HBV diversity and genotypes, as well as a small set of example sequences which are useful for trialling the web interface.

> NOTE: Many of the files - especially `qc.json` - still need to be optimised for HBV.

### Updating the tree

Generate a new tree, and manually check the clades have been correctly annotated - this often requires updating `config/clades-genotypes.tsv` and re-generating the Auspice json:

```bash
snakemake --cores 2 auspice/hbv_nextclade-tree.json
```

When you are happy with the tree, create a new datestamped version of the dataset by copying a previous directory in `nextclade_datasets`
and replacing the `tree.json` with the newly generated tree (`auspice/hbv_nextclade-tree.json`).

Finally update `nextclade_dataset` in `config/config.yaml` to point to the new dataset


### Updating example sequences

A small set of example sequences can be generated via

```bash
snakemake --cores 1 results/nextclade-sequences/filtered.fasta
```

And then a new dataset created as described above.

### Changing the Nextclade reference


* Create a new dataset by copying files from a previous one and changing the path to reflect the new reference
* Remove the tree.json, change the reference related files
* update `nextclade_dataset` in `config/config.yaml` to point to the new dataset
* run through a "nextclade-tree" build, e.g. `snakemake -npf auspice/hbv_nextclade-tree.json`
* update the clade TSV & rerun until you're happy (all clade definitions will need updating if you changed the reference)
* copy `auspice/hbv_nextclade-tree.json` to be the `tree.json` in your dataset

# TODO

* Sanitise metadata
* Check occurrences of a clade of (genbank-labelled) genotype G within the I clade
* Optimise nextclade parameters
* The mutation counts are incredible & need examining
