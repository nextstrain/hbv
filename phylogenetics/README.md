# nextstrain.org/hbv

This is the experimental Nextstrain phylogenetic workflow behind the (as yet unreleased) HBV datasets.

## Usage

> NOTE: These command examples assume you are within the `phylogenetic` directory.

```
snakemake --cores 4 -pf
```

The following build is automatically generated (change the mask between P and S genes in config file):

- `auspice_datasets/P_masked/main-clades.json` the entire human-HBV tree, with appr. 2k tips using stitched genotypes

The following alternative builds are additionally generated in dev mode (set in config file):

- `auspice/P_masked/full-tree.jso` the entire human-HBV tree, without the stitching procedure
- `auspice/hbv_{A,B,C,D,E,F,G,H,I}.json` for single genotype builds. Note that some of these are very small and one should consider disabling filtering by subgenotype annotation availability (via config)

## Configuration

_Work in progress_

### Input data

The phylogenetics workflow expects a number of files to exist which are produced by the ingest workflow.
Please see `../ingest/README.md` for how to generate these files.

##  Updating the Nextclade dataset

Nextclade datasets exist for reference `NC_003977`. The dataset includes a ~2000 tip tree attempting to cover observed human HBV diversity and genotypes, as well as a small set of example sequences which are useful for trialling the web interface.

> NOTE: Note that alignment parameters and QC metrics are set to those suggested for highly diverse viruses and not adapted for HBV specifically!

#### Updating the tree

```bash
snakemake --cores 4 --configfile defaults/nextclade/config_nextclade.yaml
```

This creates a new version of the dataset in `../nextclade_datasets`

Specify a new date-stamped version name in version in the config file and finally update `nextclade_dataset` in the ingest `../ingest/default/config.yaml` to point to the new dataset

#### Updating example sequences

A small set of example sequences is automatically generated in the nextclade workflow and can be configurated in the nextclade config file.
Example sequences are sampled evenly across genotypes, recombinants, and qc status but can be alternatively sampled at random (config option).
