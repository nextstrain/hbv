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



def define_filters(mode, key):
    
    # Subsample based on dev mode and build     
    if str(config.get("dev", False)).lower() in ("1", "true", "yes"):
        max_n = int(config.get("dev_n_stitched_parts", 100)) if mode == "stitched" else int(config.get("dev_n_totaltree", 500))
    else:
        max_n = int(config.get("n_stitched_parts", 800)) if mode == "stitched" else int(config.get("n_totaltree", 3000))

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


    if config.get("sampling", {}).get("augur_custom_filter"):
        query_exprs.append(config["sampling"]["augur_custom_filter"])    

    args = []
    if query_exprs:
        combined = " & ".join(f"({q})" for q in query_exprs)
        args.append(f"--query '{combined}'")
    
    args.append(
    f"--group-by {' '.join(config['sampling']['group_by'])} "
    f"--subsample-max-sequences {max_n}"
)

    return " ".join(args)

    #elif wildcards.build == "nextclade-tree":
    #    return "--group-by genotype_genbank --subsample-max-sequences 2000"
    #elif wildcards.build == "nextclade-sequences":
    #    return "--group-by genotype_genbank --subsample-max-sequences 25"


def get_filter_args(wc):
    return define_filters(wc.mode, wc.key)


PREVIOUS_EXCLUDE_FILE = f"defaults/{config['workflow']}/exclude.txt"


def get_previous_exclude_file(wildcards):
    if config.get("filter_previously_excluded", False):
        return PREVIOUS_EXCLUDE_FILE
    return []


rule combine_previous_excludes:
    output:
        exclude=PREVIOUS_EXCLUDE_FILE,
    params:
        source_dir=RESULTS,
    shell:
        r"""
        if [ -d "{params.source_dir}" ]; then
          find "{params.source_dir}" -type f -name '*_exclude.txt' -exec cat {{}} + \
            | tr -s '[:space:]' '\n' \
            | awk 'NF' \
            | sort -u \
            > "{output.exclude}"
        else
          : > "{output.exclude}"
        fi

        n=$(wc -l < "{output.exclude}" | tr -d ' ')
        echo "$n unique accessions written to {output.exclude}"
        """



rule filter_by_clade:
    input:
        alignment="data/filtered/alignment.len_filtered.fasta",
        metadata="data/filtered/metadata.len_filtered.tsv",
        exclude = get_previous_exclude_file,
    output:
        alignment=RESULTS + "/{mode}/{key}/filtered.fasta",
        metadata=RESULTS + "/{mode}/{key}/filtered.tsv",
    params:
        args=get_filter_args,
        filter_previously_excluded=str(config.get("filter_previously_excluded", False)).lower(),
    wildcard_constraints:
        mode="basic|stitched|single-clade",
        key="all|" + "|".join(ALL_GTS),
    shell:
        r"""
        exclude_arg=""
        if [ "{params.filter_previously_excluded}" = "true" ]; then
          if [ -s "{input.exclude}" ]; then
            exclude_arg="--exclude {input.exclude}"
          else
            echo "No previous exclude accessions found in {input.exclude}; continuing without --exclude"
          fi
        fi

        augur filter \
          --sequences {input.alignment} --metadata {input.metadata} \
          --metadata-id-columns accession \
          $exclude_arg \
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

rule write_gene_mask:
    input:
        regions="defaults/genomic_regions_genes.txt",
        ref_gb=config["reference"]["genbank"],
        script="scripts/write_gene_mask.py",
    output:
        mask=temp(RESULTS + "/masks/{gene}_mask.txt"),
    wildcard_constraints:
        gene="|".join(config["gene_mask"]),
    shell:
        r"""
        python {input.script} \
          --regions "{input.regions}" \
          --reference-genbank "{input.ref_gb}" \
          --gene "{wildcards.gene}" \
          --output "{output.mask}"
        """


rule mask_gene:
    input:
        alignment= RESULTS + "/{mode}/{key}/filtered.fasta",
        mask=RESULTS + "/masks/{gene}_mask.txt",
    output:
        alignment=RESULTS + "/{mode}/{key}/{gene}_masked/{gene}_masked_aln.fasta",
    wildcard_constraints:
        gene="|".join(config["gene_mask"]),
    shell:
        r"""
        augur mask \
          --sequences "{input.alignment}" \
          --mask "{input.mask}" \
          --output "{output.alignment}"
        """
