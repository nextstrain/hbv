# Nextclade datasets

- [`references/`](./references/NC_003977/versions/) -  Temporary folder for storing nextclade datasets for currently reference genome NC_003977
- `dataset.json` - top-level Nextclade dataset descriptor defining the dataset name and default reference

Available are:
- [`2023-08-22/`](./references/NC_003977/versions/2023-08-22) Previous version with recombination-unaware phylogeny - by James Hadfiled

- The phylogenetics workflow will automatically create a new dataset here for a build defined by the config file if called with

> NOTE: This command assumes you are within the `phylogenetics` directory.

```bash
snakemake --cores 4 --configfile defaults/nextclade/config_nextclade.yaml
```

To view the generated dataset, navigate to the dataset version e.g. [`2026_04_jonas/stitched_P_masked/`](./references/NC_003977/versions/2026_04_jonas/stitched_P_masked) and exectute:

> NOTE: This command assumes you are within the proper dataset directory within [`references/`](./references/NC_003977/versions/) .

```bash
auspice view --datasetDir test_output
```

Then select nextclade.auspice and filter by "Node Type → New" to see example sequences.
