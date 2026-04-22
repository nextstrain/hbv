# Exploration of HBV phylogenetic changes between trees of different genomic segments

> This folder is highly experimental, i.e. not part of the nextstrain procedure, and should be removed before publishing
>
> Ingest dev mode is incompatible with other directories. First run ingest properly before rerunning phylogenetics or exploratory.

Generate naive (recombination-unaware) phylogenies of different regions of HBV with:

```bash
snakemake --cores 4
```
Regions can be specified in the config file; currently, the HBV genes and breakpoints identified by an early study are identified.

Use genes by default:

```bash
snakemake --cores 1 comparison_targets --config dev_mode=true
```

Use breakpoint-derived segments instead:

```bash
snakemake --cores 1 comparison_targets --configfile defaults/config.segments.yaml --config dev_mode=true --rerun-incomplete
```

> The full exploratory trees are large and there is no automatic subsampling in normal mode, so a full run can take a long time. Use `dev_mode=true` for faster test runs.

To run only the comparison workflow once the trees exist, use the dedicated aggregate target:

```bash
snakemake --cores 4 comparison_targets
```

For faster testing, dev mode keeps 300 sequences for tree building and uses 10 sampled tip pairs/shared tips in the comparison metrics.

The comparison outputs can also be targeted individually. For the current breakpoint-based segment setup, use:

```bash
snakemake --cores 4 results/segments/pw_distances/corr_matrix/pw_tip_distance_correlation_matrix.patristic.png
snakemake --cores 4 results/segments/rf/RF.pdf
snakemake --cores 4 results/segments/treeknit/runs/segment1_segment2/results_summary.txt
snakemake --cores 4 results/segments/treeknit/compare_trees_treeknit.csv
snakemake --cores 4 results/segments/treeknit/plots/.done
```

If the config is switched to genes instead of segments, replace `results/segments/...` with `results/genes/...`.

> NOTE THAT THESE ARE NOT FINAL SCIENTIFIC RESULTS! The breakpoint-lifting procedure between references needs to be reviewed again, and the current results cannot be trusted! Running the rules with the current tip sampling numbers takes very long so change these in the config file beforehand.
