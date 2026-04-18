
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
        length_filtering = config["length_filtering"]["use_as_filter"],
        min_length       = config["length_filtering"]["min_length"],
        max_length       = config["length_filtering"]["max_length"],
    shell:
        r"""
        # mkdir -p data
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
# TODO include nextclade build back here 
def define_filters(mode, key):
    
    # Subsample based on dev mode and build     
    if str(config.get("dev", False)).lower() in ("1", "true", "yes"):
        max_n = int(config .get("dev_n_stitched_parts", 100)) if mode == "stitched" else int(config .get("dev_n_totaltree", 500))
    else: 
        max_n = int(config .get("n_stitched_parts", 800)      if mode == "stitched" else int(config .get("n_totaltree", 3000))


    query_exprs = []

    if mode in ("stitched", "single-clade"):
        if key == "C":
            query_exprs.append('(clade_nextclade=="C") | (clade_nextclade=="C_re")')
        else:
            query_exprs.append(f'clade_nextclade=="{key}"')

        if mode=="stitched" and config.get("subgenotype_filtering_fulltree", False) and key in ["A", "B", "C", "D", "F"]: #"I"
            query_exprs.append('subgenotype_genbank.notnull() & (subgenotype_genbank != "")')
        if mode=="single-clade" and config.get("subgenotype_filtering_singleclade", False) and key in ["A", "B", "C", "D", "F"]: #"I"
            query_exprs.append('subgenotype_genbank.notnull() & (subgenotype_genbank != "")')

    elif mode != "basic":
        raise Exception("Unknown build parameter")


    if config.get("filter_genbank_vs_nextclade", False):
        query_exprs.append('genotype_genbank.notnull() & (genotype_genbank != "") & clade_nextclade.notnull() & (clade_nextclade != "")')
        if mode in ("stitched", "single-clade"):
            if key == "C":
                query_exprs.append('genotype_genbank == "C"')
            else:
                query_exprs.append('genotype_genbank == clade_nextclade')


    if config.get("augur_custom_filter"):
        query_exprs.append(config["augur_custom_filter"])    

    args = []
    if query_exprs:
        combined = " & ".join(f"({q})" for q in query_exprs)
        args.append(f"--query '{combined}'")
    
    args.append(f"--group-by year region --subsample-max-sequences {max_n}") 

    return " ".join(args)

    #elif wildcards.build == "nextclade-tree":
    #    return "--group-by genotype_genbank --subsample-max-sequences 2000"
    #elif wildcards.build == "nextclade-sequences":
    #    return "--group-by genotype_genbank --subsample-max-sequences 25"




# # TODO right now this is empty! 
# rule include_file:
#     output:
#         file="results/{mode}/{key}/include.txt",
#     params:
#         #ref=lambda wc: config["reference"]["id"],
#     shell:
#         r"""
#         mkdir -p results/{wildcards.mode}/{wildcards.key}
#         touch {output.file}
#         """
#         #printf "%s\n" "{params.ref}" > {output.file}



rule filter_by_clade:
    input:
        alignment="data/filtered/alignment.len_filtered.fasta",
        metadata="data/filtered/metadata.len_filtered.tsv",
        #include = "results/{mode}/{key}/include.txt",   
        #exclude = "defaults/exclude.txt",
    output:
        alignment="results/{mode}/{key}/filtered.fasta",
        metadata="results/{mode}/{key}/filtered.tsv",
    params:
        args=lambda wc: define_filters(wc.mode, wc.key), # accept wc , in function get mode key
    wildcard_constraints:
        mode="basic|stitched|single-clade",
        key="all|" + "|".join(ALL_GTS),
    shell:
        r"""
        # mkdir -p results/{wildcards.mode}/{wildcards.key}
        augur filter \
          --sequences {input.alignment} --metadata {input.metadata} \
          --metadata-id-columns accession \
          {params.args} \
          --output-sequences {output.alignment} \
          --output-metadata {output.metadata} 
        """
        #--include {input.include} \

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

rule mask_gene:
    input:
        regions="defaults/genomic_regions_genes.txt",
        alignment="results/{mode}/{key}/filtered.fasta",
    output:
        alignment="results/{mode}/{key}/{gene}_masked/{gene}_masked_aln.fasta",
    params:
        ref=config["reference"]["id"]
    wildcard_constraints:
        gene="|".join(config["gene_mask"]),
    shell:
        r"""
        # mkdir -p data/masked

        read -r start end < <(
        awk -v g="{wildcards.gene}" '$0 !~ /^#/ && $1==g {{print $2, $3; exit}}' "{input.regions}"
        )

        L=$(seqkit grep -n -p "^{params.ref}$" "{input.alignment}" \
            | seqkit seq -s \
            | head -n1 \
            | tr -d '\n' \
            | wc -c \
            | tr -d ' ')

        if [ -z "$L" ]; then
        L=$(seqkit seq -s "{input.alignment}" \
            | head -n1 \
            | tr -d '\n' \
            | wc -c \
            | tr -d ' ')
        fi

        if [ -z "$L" ]; then
        echo "mask_gene: could not determine alignment length" >&2
        exit 1
        fi
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


