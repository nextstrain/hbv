# Nextstrain Hepatitis B builds

Currently a WIP

Based on Katie Kistler's work in [blab/adaptive-evolution](https://github.com/blab/adaptive-evolution)

## Ingest

1. Obtain `./ingest/genbank_sequences.gb`
    - example search: [(Hepatitis B virus) AND (complete genome)](https://www.ncbi.nlm.nih.gov/nuccore/?term=(Hepatitis+B+virus)+AND+(complete+genome)), which returns ~11k samples
      - Could also explore `"Hepatitis B virus"[porgn:__txid10407]`
    - download as "Complete Record" / "File" / "GenBank (full)", a ~100Mb file.

2. `./ingest/curate_genbank_to_fasta.ipynb`
    - produces `./ingest/hepatitisB_all.fasta`


## Phylo

1. `snakemake --cores 4 -pf auspice/hepatitisB_all.json`


