
rule length_filter:
    input:
        sequences = "../ingest/results/sequences.fasta",
        metadata  = "../ingest/results/metadata.tsv",
        alignment = "../ingest/results/aligned.fasta",
    output:
        sequences = "data/filtered/sequences.len_filtered.fasta",
        alignment = "data/filtered/alignment.len_filtered.fasta",
        metadata  = "data/filtered/metadata.len_filtered.tsv",
    params:
        min_length = config["min_length"],
        max_length =config["max_length"],
        length_filtering=config["length_filtering"],
    shell:
        r"""
        mkdir -p data
        kept_ids="$(mktemp)"

        if [ "{params.length_filtering}" = "true" ] || [ "{params.length_filtering}" = "True" ] || [ "{params.length_filtering}" = "1" ]; then
            seqkit seq -g -m={params.min_length} -M={params.max_length} {input.sequences} > {output.sequences}

            seqkit seq -n {output.sequences} > "$kept_ids"

            seqkit grep -f "$kept_ids" {input.alignment} > {output.alignment}

            awk 'NR==FNR {{a[$1]=1; next}} FNR==1 || ($1 in a)' "$kept_ids" {input.metadata} > {output.metadata}
        else
            cp {input.sequences} {output.sequences}
            cp {input.alignment} {output.alignment}
            cp {input.metadata} {output.metadata}
        fi

        rm -f "$kept_ids"
        """


## TODO - there are a number of nextclade QC status' we can filter on here.
## Currently the settings in the nextclade dataset need to be looked at as 
## around 40% of all sequences (including the entirety of some genotypes)
## have QC=bad mainly due to frameshifts and stop codons.
# TODO include dev and nextclade build back here 

def args_for(mode, key):
    # dev subsampling
    dev = str(config.get("dev", False)).lower() in ("1", "true", "yes")
    dev_n = int(config.get("dev_n", 100))
    max_n = dev_n if dev else 4000

    if mode == "basic":
        # no genotype query, just (optional) balanced subsampling
        return f"--group-by year --subsample-max-sequences {max_n}"

    if mode == "stitched" or mode=="single-clade":
        # per-genotype query + balanced subsampling
        build = key
        if build == "C":
            query = """--query "(clade_nextclade=='C') | (clade_nextclade=='C_re')\""""
        else:
            query = f"""--query "clade_nextclade=='{build}'\""""
        return f"{query} --group-by year --subsample-max-sequences {max_n}"

    #elif wildcards.build == "dev":
    #    return "--group-by genotype_genbank --subsample-max-sequences 500"
    #elif wildcards.build == "nextclade-tree":
    #    return "--group-by genotype_genbank --subsample-max-sequences 2000"
    #elif wildcards.build == "nextclade-sequences":
    #    return "--group-by genotype_genbank --subsample-max-sequences 25"

    raise Exception("Unknown build parameter")


# TODO include more here
rule include_file:
    output:
        file="results/{mode}/{key}/include.txt",
    params:
        ref=lambda wc: config["reference"]["id"],
    shell:
        r"""
        mkdir -p results/{wildcards.mode}/{wildcards.key}
        printf "%s\n" "{params.ref}" > {output.file}
        """


rule filter_by_clade:
    input:
        alignment="data/filtered/alignment.len_filtered.fasta",
        metadata="data/filtered/metadata.len_filtered.tsv",
        include = "results/{mode}/{key}/include.txt",   
        #exclude = "defaults/exclude.txt",
    output:
        alignment="results/{mode}/{key}/filtered.fasta",
        metadata="results/{mode}/{key}/filtered.tsv",
    params:
        args=lambda wc: args_for(wc.mode, wc.key),
    wildcard_constraints:
        mode="basic|stitched|single-clade",
        key="all|" + "|".join(ALL_GTS),
    shell:
        r"""
        mkdir -p results/{wildcards.mode}/{wildcards.key}
        augur filter \
          --sequences {input.alignment} --metadata {input.metadata} \
          --metadata-id-columns accession \
          --include {input.include} \
          {params.args} \
          --output-sequences {output.alignment} \
          --output-metadata {output.metadata} 
        """

#____________________________________________________________________________________________________________________________________________________________________________________________


rule specify_genomic_regions_genes:
    input:
        ref_gb=config["reference"]["genbank"],
        script="scripts/specify_genomic_regions_genes.py",
    output:
        regions="defaults/genomic_regions_genes.txt",
    shell:
        r"""
        python {input.script} {input.ref_gb} {output.regions} 
        """

#TODO make this work for CDS specifically and not for genes! then keep same naming conventions as in nextclade

# TODO have this work for other genes as well in different builds (e.g. S)
rule mask_gene:
    input:
        regions="defaults/genomic_regions_genes.txt",
        alignment="results/{mode}/{key}/filtered.fasta",
    output:
        alignment="results/{mode}/{key}/masked.{gene}.fasta",
    params:
        gene_mask=config["gene_mask"],
        ref=config["reference"]["id"]
    wildcard_constraints:
        gene="|".join(config["genes"]),
    shell:
        r"""
        mkdir -p data/masked

        read -r start end < <(
        awk -v g="{params.gene_mask}" '$0 !~ /^#/ && $1==g {{print $2, $3; exit}}' "{input.regions}"
        )

        L=$(seqkit grep -n -p "^{params.ref}$" "{input.alignment}" | seqkit seq -s | tr -d '.-' | wc -c | tr -d ' ')
        maskfile="$(mktemp)"

        if [ "$start" -le "$end" ]; then
        # mask [1,start-1]
        if [ "$start" -gt 1 ]; then
            seq 1 $((start-1)) >> "$maskfile"
        fi
        # mask [end+1,L]
        if [ "$end" -lt "$L" ]; then
            seq $((end+1)) "$L" >> "$maskfile"
        fi
        else
        # gene wraps -> mask [end+1, start-1]
        seq $((end+1)) $((start-1)) >> "$maskfile"
        fi

        augur mask \
        --sequences "{input.alignment}" \
        --mask "$maskfile" \
        --output "{output.alignment}"

        rm -f "$maskfile"
        """


