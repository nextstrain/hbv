rule fetch_genbank:
    params:
        term = 'Hepatitis+B+virus[All+Fields]complete+genome[All+Fields]'
    output:
        genbank = "ingest/data/genbank.gb",
    retries: 1  # Requires snakemake 7.7.0 or later
    shell:
        """
        ingest/scripts/fetch-genbank.py --term {params.term:q} --output {output.genbank}
        """

rule parse_genbank:
    input:
        genbank = "ingest/data/genbank.gb",
    output:
        sequences = "ingest/results/genbank.fasta",
        metadata = "ingest/results/genbank.tsv",
    shell:
        """
        ingest/scripts/parse-genbank.py --genbank {input.genbank} \
            --output-sequences {output.sequences} --output-metadata {output.metadata}
        """

rule re_circularise:
    input:
        metadata = "ingest/results/genbank.tsv",
        sequences = "ingest/results/genbank.fasta",
    output:
        sequences = "ingest/results/genbank.recircular.fasta",
        metadata = "ingest/results/genbank.recircular.tsv",
    params:
        reference = "JN182318" # should be the same as the one used for the big alignment
    shell:
        """
        ingest/scripts/re-circularise.py \
            --seqs-in {input.sequences} --meta-in {input.metadata} \
            --seqs-out {output.sequences} --meta-out {output.metadata} \
            --reference {params.reference}
        """

rule align_everything:
    input:
        sequences = "ingest/results/genbank.recircular.fasta",
        metadata = "ingest/results/genbank.recircular.tsv",
        reference = "ingest/config/JN182318.fasta"
    output:
        log = temp('ingest/results/nextalign-errors.txt'),
        alignment = "ingest/results/aligned.fasta",
        metadata = "ingest/results/metadata.tsv",
    params:
        match_rate = "0.2" # default is 0.3
    threads: 4
    shell:
        """
        nextalign run \
            --output-fasta {output.alignment} --output-errors {output.log} \
            --min-match-rate {params.match_rate} -r {input.reference} \
            --replace-unknown -j {threads} {input.sequences};
        ingest/scripts/append-nextalign-log.py \
            --meta-in {input.metadata} --meta-out {output.metadata} --nextalign {output.log}
        """

rule provision:
    input:
        sequences = "ingest/results/genbank.recircular.fasta",
    output:
        sequences = "ingest/results/sequences.fasta",
    shell:
        """
        cp {input.sequences} {output.sequences}
        """