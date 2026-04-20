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

Nextclade datasets exist for reference `NC_003977`. The dataset includes a ~2000 tip tree attempting to cover observed human HBV diversity and genotypes, as well as a small set of example sequences which are useful for trialling the web interface.

> NOTE: Many of the files - especially `qc.json` - still need to be optimised for HBV.

#### Updating the tree


```bash
snakemake --cores 4 configfile defaults/config_nextclade.yaml
```

This creates a new version of the dataset in `../nextclade_datasets`

Specify a new date-stamped version name in version in the config file and finally update `nextclade_dataset` in the ingest `../ingest/default/config.yaml` to point to the new dataset


#### Updating example sequences

A small set of example sequences is automatically generated in the nextclade workflow and can be configurated in the nextclade config file

And then a new dataset created as described above.