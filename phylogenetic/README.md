# nextstrain.org/hbv

This is the experimental Nextstrain phylogenetic workflow behind the (as yet unreleased) HBV datasets.

## Usage

> NOTE: These command examples assume you are within the `phylogenetic` directory.


```
snakemake --cores 4 -pf {target}
```

Where target is one (or more) of the following auspice datasets:
- `auspice/hbv_dev.json` for a small ~500-tip dev tree
- `auspice/hbv_{A,B,C,D}.json` for genotype builds, each of ~500 tips. Currently only genotypes A-D are supported.
- `auspice/hbv_all.json` the entire human-HBV tree, with 11k tips (takes ~15min on a 4-core M1 machine)


## Configuration

_Work in progress_

### Input data

The phylogenetics workflow expects a number of files to exist which are produced by the ingest workflow.
Please see `../ingest/README.md` for how to generate these files.

##  Updating the Nextclade dataset

> The following is out of date and remains here as a reminder for me to update!

Nextclade datasets exist for reference `NC_003977`. The dataset includes a ~2000 tip tree attempting to cover observed human HBV diversity and genotypes, as well as a small set of example sequences which are useful for trialling the web interface.

> NOTE: Many of the files - especially `qc.json` - still need to be optimised for HBV.

#### Updating the tree

Generate a new tree, and manually check the clades have been correctly annotated - this often requires updating `defaults/clades-genotypes.tsv` and re-generating the Auspice json:

```bash
snakemake --cores 2 auspice/hbv_nextclade-tree.json
```

When you are happy with the tree, create a new datestamped version of the dataset by copying a previous directory in `../nextclade_datasets`
and replacing the `tree.json` with the newly generated tree (`auspice/hbv_nextclade-tree.json`).

Finally update `nextclade_dataset` in `config/config.yaml` to point to the new dataset


#### Updating example sequences

A small set of example sequences can be generated via

```bash
snakemake --cores 1 results/nextclade-sequences/filtered.fasta
```

And then a new dataset created as described above.

#### Changing the Nextclade reference


* Create a new dataset by copying files from a previous one and changing the path to reflect the new reference
* Remove the tree.json, change the reference related files
* update `nextclade_dataset` in `config/config.yaml` to point to the new dataset
* run through a "nextclade-tree" build, e.g. `snakemake -npf auspice/hbv_nextclade-tree.json`
* update the clade TSV & rerun until you're happy (all clade definitions will need updating if you changed the reference)
* copy `auspice/hbv_nextclade-tree.json` to be the `tree.json` in your dataset
