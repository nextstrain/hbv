# Legacy code within the phylogenetic workflow 

> Contains previous work by James Hadfield not yet fully incoporated in the current workflow


## What is still useful here

Most of the core tree-building logic now exists elsewhere in the repo, so this `legacy/` folder is mainly useful as a reference for older HBV-specific workflow ideas which have not been fully ported into the current main procedure.

Useful files to look at:

- `phylogenetic/defaults/clades-genotypes.tsv`: manually curated genotype-defining mutations for use with `augur clades`: For manually redefining clade annotation for future nextstrain datasets
- `phylogenetic/defaults/color_ordering.tsv` and `phylogenetic/defaults/color_schemes.tsv`: the old Auspice coloring setup, where category order and palette assignment were generated automatically instead of only being described in a static config JSON.
- `phylogenetic/defaults/lat-longs.tsv`: explicit lat/long mapping for Augur export.
- `phylogenetic/scripts/attach_root_mutations.py`: adds all root-vs-reference nucleotide and amino-acid differences back onto the root before export.

## Notebook overview

There are two notebook-style resources:

- `notebook/`: an Observable Framework website intended for GitHub Pages. It is a polished presentation layer rather than an analysis scratchpad. In particular, `notebook/docs/rotated.md` shows why rotating circular HBV genomes improves apparent alignment coverage, and `notebook/docs/index.md` is the landing page for that mini-site.
- `notebooks/alignment-qc.ipynb`: a Jupyter notebook for local exploratory QC (using old reference!). It is a more direct Python/matplotlib analysis notebook which reads metadata and alignments, counts reference-matching sites, and plots alignment quality by genotype/host.

## Legacy Snakemake functionality not yet present in the main workflow


- Generated Auspice color tables: the `colors` rule plus `assign-colors.py` created a `colors.tsv` dynamically from metadata, category ordering, and palette files. This is more sophisticated than the current static Auspice config setup and is the main legacy reference to check if you want to reintroduce manual categorical color control.
- `augur clades` from a manual TSV: the legacy pipeline inferred genotype labels from `clades-genotypes.tsv`. In the current workflow this has effectively been replaced by the custom `assign_clades.py` tree/metadata-based approach, so the old TSV-based clade definition system is still a distinct reference.
- Forced include logic before filtering: `include_file` explicitly keeps the chosen root, configured outgroups, and the reference accession in broad builds. The current workflow has filtering and downsampling logic, but not this exact "always include these sequences" file-based step.
- Root mutation attachment before export: the legacy `attach_root_mutations` rule ensured that all differences between the alignment reference and the inferred root were attached to the root node. This exact step is not part of the current main workflow.
- Richer build-specific Auspice config generation: the legacy `auspice_config` rule generated JSON with different default colorings, branch labels, filters, and extra QC/circularisation colorings depending on the build. The current workflow already has Auspice config JSONs, but they are simpler and mostly static.
- Explicit `lat-longs.tsv` use during `augur export`: the legacy workflow passed both colors and geographic coordinates directly into export. The current workflow exports with static config plus metadata-driven colorings, but not this old explicit color-table/lat-long-table combination.



For additional context on how the workflow evolved, it is also worth checking the earlier 2023 version of the repo.
