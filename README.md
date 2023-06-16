# Nextstrain Hepatitis B builds

Currently a WIP

Based on Katie Kistler's work in [blab/adaptive-evolution](https://github.com/blab/adaptive-evolution)

## Obtain data

1. Obtain `./ingest/genbank_sequences.gb`
    - example search: [(Hepatitis B virus) AND (complete genome)](https://www.ncbi.nlm.nih.gov/nuccore/?term=(Hepatitis+B+virus)+AND+(complete+genome)), which returns ~11k samples
      - Could also explore `"Hepatitis B virus"[porgn:__txid10407]`
    - download as "Complete Record" / "File" / "GenBank (full)", a ~100Mb file.


## Ingest

You can run the ingest part of the snakemake pipeline via:

```
snakemake --cores 4 -npf ingest/results/aligned.fasta ingest/results/sequences.fasta ingest/results/metadata.tsv
```


## Phylo

For a small ~800-tip dev tree: `snakemake --cores 4 -pf auspice/hbv_dev.json`
For the entire human-HBV tree, with 11k tips (takes ~15min on a 4-core M1 machine): `snakemake --cores 4 -pf auspice/hbv_all.json`

(If you haven't run the ingest section of the pipeline then any phylo build will run it.)

