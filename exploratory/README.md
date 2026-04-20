# Exploration of HBV phylogenetic changes between trees of different genomic segments 



> This folder is highly experimental, i.e. not part of the nextstrain procedure, and should be removed before publishing 

Generate naive (recombination-unaware) phylogenies of different regions in of HBV with:

```bash
snakemake --cores 4 
```
Regions can be specified in the config file, currently the HBV genes and breakpoints identified by an early study are identified.

The trees can be compared by seperately targeting the rules in  [`./rules/tree_comparisons.smk`](./rules/tree_comparisons.smk) for different metrics 


> NOTE THAT THESE ARE NOT FINAL SCIENTIFIC RESULTS! The breakpoint lifting procedure between references needs to be reviewed again and the current results cannot be trusted! Running the rules with current tip sampling numbers takes very long so change these in the config file beforehand. Please note that the results of the workflow and calculated metrics are available for reference (/auspice_darasets, /results respectively) so that the repo does not have to be run!