from os.path import join

rule generate_example_sequences:
    input:
        metadata="../ingest/results/metadata.tsv",
        sequences="../ingest/results/sequences.fasta",
        script="scripts/select_examples_nextclade.py"
    output:
        sequences="results/nextclade/example_sequences/sequences.fasta",
        selected="results/nextclade/example_sequences/selected_examples.tsv"
    params:
        n_total=config["example_sequences"]["n"],
        recomb_frac=config["example_sequences"]["recombinant_fraction"],
        sampling=config["example_sequences"]["sampling"]
    shell:
        r"""
        python {input.script} \
          --metadata {input.metadata} \
          --output {output.selected} \
          --n-total {params.n_total} \
          --recomb-frac {params.recomb_frac}\
          --sampling {params.sampling}

        tail -n +2 {output.selected} | cut -f1 | seqkit grep \
          --pattern-file /dev/stdin \
          {input.sequences} \
          > {output.sequences}
        """

rule assemble_dataset:
    input:
        tree=join(RESULTS, config["nextclade_tree_build"]),
        sequences = "results/nextclade/example_sequences/sequences.fasta",
        reference = "defaults/nextclade/reference.fasta",
        annotation= config["reference"]["gff"],
        pathogen = "defaults/nextclade/pathogen.json"
    output:
        tree=       DATASET_DIR + "tree.json",
        annotation= DATASET_DIR + "genome_annotation.gff3",
        readme=     DATASET_DIR + "README.md",
        changelog=  DATASET_DIR + "CHANGELOG.md",
        reference=  DATASET_DIR + "reference.fasta",
        sequences=  DATASET_DIR + "sequences.fasta",
        pathogen =  DATASET_DIR + "pathogen.json"
    shell:
        """
        cp {input.tree} {output.tree}
        cp {input.annotation} {output.annotation}
        cp {input.reference} {output.reference}
        cp {input.sequences} {output.sequences}
        cp {input.pathogen} {output.pathogen}
        printf "# Example dataset for HepB virus\n\nDataset for Hepatitis B Virus. Work in progress. \n\nNote that alignment parameters are set to those suggested for highly diverse viruses and not adapted for HBV specifically." > {output.readme}
        printf "## Unreleased\n\nInitial release.\n" > {output.changelog}
        """

rule test_dataset:
    input:
        sequences=  DATASET_DIR + "sequences.fasta",
        tree=       DATASET_DIR + "tree.json",
        annotation= DATASET_DIR + "genome_annotation.gff3",
        readme=     DATASET_DIR + "README.md",
        changelog=  DATASET_DIR + "CHANGELOG.md",
        reference=  DATASET_DIR + "reference.fasta",
        pathogen =  DATASET_DIR + "pathogen.json"
    output:
        outdir=directory(DATASET_DIR + "test_output"),
    params:
        dataset_dir=DATASET_DIR,
    shell:
        """
        nextclade3 run \
            {input.sequences} \
            --input-dataset {params.dataset_dir} \
            --output-all {output.outdir}
        """

#______________________________________________________________________________________________________________________________________________________________________________________________
#______________________________________________________________________________________________________________________________________________________________________________________________

# TODO: UPLOAD DATASETS TO NEXTSTRAIN.ORG

# Does not run by default as part of rule all
rule deploy_to_nextstrain_staging:
    input:
        rules.all.input
    shell:
        """
        nextstrain deploy s3://nextstrain-staging {input}
        """

rule download:
   "Downloading ingested sequences and metadata from data.nextstrain.org"
   output:
       sequences="nextclade/data/sequences.fasta.zst",
       metadata="nextclade/data/metadata.tsv.zst",
       alignment="nextclade/data/alignment.fasta.zst",
   params:
       metadata_url="https://data.nextstrain.org/files/workflows/hbv/metadata.tsv.zst",         # Those do not exist yet
       sequences_url="https://data.nextstrain.org/files/workflows/hbv/sequences.fasta.zst",
       alignment_url="https://data.nextstrain.org/files/workflows/hbv/alignment.fasta.zst",
   shell:
       """
       curl -fsSL {params.sequences_url:q} --output {output.sequences}
       curl -fsSL {params.metadata_url:q} --output {output.metadata}
       curl -fsSL {params.alignment_url:q} --output {output.alignment}
       """
