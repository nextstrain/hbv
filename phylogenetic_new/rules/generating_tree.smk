
# Reference sequence is to be included in the alignment but excluded from tree building to avoid it appearing in wrong individual clade-trees.
rule augur_tree:
    input:
        alignment = "results/{mode}/{key}/{gene}_masked/{gene}_masked_aln.fasta",
    output:
        tree      = "results/{mode}/{key}/{gene}_masked/{gene}_masked.tree.nwk",
    threads: 4
    shell:
        """
        augur tree \
          --alignment {input.alignment} \
          --method fasttree \
          --output {output.tree}
        """


# TO DO: Pruned nodes are written to exclude files that could be used in future runs directly. This is not implemented yet however. 
rule prune_tree:
    input:
        tree_nwk="results/{mode}/{key}/{gene}_masked/{gene}_masked.tree.nwk",
        metadata="results/{mode}/{key}/filtered.tsv",
        script="scripts/prune_trees.py"
    output:
        out_tree="results/{mode}/{key}/{gene}_masked/{gene}_masked.pruned.tree.nwk",
        metadata="results/{mode}/{key}/{gene}_masked/{gene}_metadata.pruned.tsv",
        exclude= "results/{mode}/{key}/{gene}_masked/{gene}_exclude.txt",
    params:
        long_branch_threshold = config["long_branch_threshold"],
        tip_branch_threshold  = config["tip_branch_threshold"],
        metadata_col="subgenotype_genbank", # for clade purity filtering,

        purity_args=(
            f"--purity_metadata_col subgenotype_genbank "
            f"--minimal_monophyletic_purity {config['clade_purity']['minimal_monophyletic_purity']} "
            f"--maximal_monophyletic_fraction {config['clade_purity']['clade_fraction']}"
            if config["clade_purity"]["use_as_filter"] else ""
        ),

        minclade_args=(
            f"--minclade_metadata_cols genotype_genbank subgenotype_genbank "
            f"--prune_min_counts {config['clade_settings']['min_count_per_clade']} {config['clade_settings']['min_count_per_subclade']}"
            if config["clade_settings"]["min_count_mode"] == "use_for_pruning" else ""
        ),
    shell:
        r"""
        python {input.script} \
          --tree {input.tree_nwk} \
          --metadata_in {input.metadata} \
          --metadata_out {output.metadata} \
          --cutoff_allbranches {params.long_branch_threshold} \
          --cutoff_tips {params.tip_branch_threshold} \
          --out_tree {output.out_tree} \
          --exclude {output.exclude} \
          {params.purity_args} \
          {params.minclade_args}     
        
        n=$(tr ' ' '\n' < {output.exclude} | sed '/^$/d' | sort -u | wc -l)
        echo "$n unique accessions in {output.exclude}"
        """       



rule augur_refine:
    input:
        tree      = "results/{mode}/{key}/{gene}_masked/{gene}_masked.pruned.tree.nwk",
        alignment = "results/{mode}/{key}/{gene}_masked/{gene}_masked_aln.fasta"
    output:
        tree      = "results/{mode}/{key}/{gene}_masked/{gene}_masked_refined.tree.nwk",
    threads: 4
    shell:
        """
        augur refine \
          --tree {input.tree} \
          --alignment {input.alignment} \
          --root mid_point \
          --output-tree {output.tree}
        """


# NOTE THAT THIS MAPPING INCLUDES PRE-REGIONS IN C AND S 
rule alias_translations_for_augur:
    input:
        dir="../ingest/data/nextclade",
        pol="../ingest/data/nextclade/cds_pol.fasta",
    params: # for genes with muultiple transcripts 
        s=config["gene_products_for_ancestral"]["S"], # "envL", "envM" or "envS"
        c=config["gene_products_for_ancestral"]["C"], # "pre-capsid" or "capsid"
    output:
        p="../ingest/data/nextclade/cds_P.fasta", # maybe make those temp files...
        s="../ingest/data/nextclade/cds_S.fasta",
        c="../ingest/data/nextclade/cds_C.fasta"
    shell:
        r"""
        cp {input.pol} {output.p}
        cp {input.dir}/cds_{params.s}.fasta {output.s}
        cp {input.dir}/cds_{params.c}.fasta {output.c}
        """


# TODO take ref sequence as root here for reconstruction as well?
rule ancestral:
    input:
        tree=     "results/{mode}/{key}/{gene}_masked/{gene}_masked_refined.tree.nwk",
        alignment="results/{mode}/{key}/filtered.fasta", # Using non-masked alignment for ancestral reconstruction
        annotation= config["reference"]["gff"],
        translations=expand("../ingest/data/nextclade/cds_{g}.fasta", g=config["genes"]),
        root = "../nextclade_datasets/references/NC_003977/versions/2023-08-22/reference.fasta",   # Mutations are relative to the reference sequence
 
    output:
        node_data = "results/{mode}/{key}/{gene}_masked/ancestral/{gene}.json",
        sequences = "results/{mode}/{key}/{gene}_masked/ancestral/{gene}.fasta",
    params: 
        genes=" ".join(config["ancestral_genes"]),
        translation_pattern="../ingest/data/nextclade/cds_%GENE.fasta",

    threads: 4
    shell:
        r"""
        augur ancestral \
          --tree {input.tree} \
          --alignment {input.alignment} \
          --annotation {input.annotation} \
          --genes {params.genes} \
          --translations {params.translation_pattern} \
          --output-node-data {output.node_data} \
          --output-sequences {output.sequences} \
          --root-sequence {input.root}

        """



