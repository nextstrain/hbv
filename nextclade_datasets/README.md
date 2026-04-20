# Nextclade datasets

- [`references/`](./references/NC_003977/versions/) -  Temporary folder for storing nextclade datasets for currently reference genome NC_003977


Available are: 
- [`2023-08-22/`](./references/NC_003977/versions/2023-08-22) Previous version with recombination-unaware phylogeny - by James Hadfiled


> NOTE: These command examples assume you are within the `phylogenetics` directory.

- The phylogenetics workflow will create a new dataset here for a build defined by the config file if called with 

```bash
snakemake --cores 4 --configfile defaults/config_nextclade.yaml
```

To view the generated dataset, navigate to the dataset version e.g. [`2026_04_jonas/stitched_P_masked/`](./references/NC_003977/versions/2026_04_jonas/stitched_P_masked) and exectute:


```bash
auspice view --datasetDir test_output
```

Then select nextclade.auspice and filter by "Node Type → New" to see example sequences. 
