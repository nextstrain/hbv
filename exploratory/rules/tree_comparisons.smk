
"""
This part of the workflow to compare trees between different regions of HBV 

REQUIRED INPUTS:

...

OUTPUTS:

...


"""

TAG = make_tag(REGION_NAMES)
OUTDIR_CORR = f"results/correlation_analysis_{TAG}"


rule compare_pairwise_distances:
    input:
        trees = expand(f"{OUTDIR}/{{name}}/{{name}}.tree.nwk", name=REGION_NAMES),
    output:
        matrix = f"{OUTDIR_CORR}/corr_matrix/pw_tip_distance_correlation_matrix.patristic.png",
        grid   = f"{OUTDIR_CORR}/corr_grid/grid_patristic_scatter.pdf",
    params:
        n_pairs = config.get("n_subsamples_pairwise_distance_comparison", 100),
        outdir  = OUTDIR_CORR,
    shell:
        r"""
        set -euo pipefail
        mkdir -p {params.outdir}/corr_matrix {params.outdir}/corr_grid
        python scripts/pairwise_tip_distance.py {input.trees} {params.n_pairs} {params.outdir}
        """


rule compare_trees_RF_pair:
    input:
        t1 = f"{OUTDIR}" + "/{seg1}/{seg1}.tree.nwk",
        t2 = f"{OUTDIR}" + "/{seg2}/{seg2}.tree.nwk",
    output:
        "results/rf/{seg1}_{seg2}.tsv"
    shell:
        r"""
        set -euo pipefail
        mkdir -p results/rf
        python scripts/compare_region_trees.py {input.t1} {input.t2} {output} 
        """


rule compare_trees_RF_all:
    input:
        expand("results/rf/{seg1}_{seg2}.tsv", zip,
               seg1=[a for a,b in PAIRS],
               seg2=[b for a,b in PAIRS])
    output:
        f"results/compare_trees_RF_{TAG}.tsv"

    shell:
        r"""
        set -euo pipefail
        mkdir -p results
        head -n 1 {input[0]} > {output}
        for f in {input}; do
            tail -n +2 "$f" >> {output}
        done
        rm -rf results/rf
        """

rule display_RF:
    input:
        table=f"results/compare_trees_RF_{TAG}.tsv"
    output:
        grid= f"results/RF_{TAG}.pdf"
    shell:
        r"""
        python scripts/visualise_rf.py {input.table} {output.grid}
        """



rule tree_knit_pair:
    input:
        t1 = "data/regions/{seg1}/{seg1}.tree.nwk",
        t2 = "data/regions/{seg2}/{seg2}.tree.nwk",
    output:
        outdir = directory("results/treeknit/{seg1}_{seg2}")
    params:
        tmp1 = "results/treeknit/{seg1}_{seg2}.out1.nwk",
        tmp2 = "results/treeknit/{seg1}_{seg2}.out2.nwk",
        prune = "scripts/treeknit/prune_to_shared.py",
        summarize = "scripts/treeknit/treeknit_summarize.py",
        n_subset = config.get("treeknit", {}).get("n_subset", "None"),
        seed = config.get("treeknit", {}).get("seed", "None"),
    shell:
        r"""
        set -euo pipefail
        mkdir -p {output.outdir}

        python {params.prune} {input.t1} {input.t2} {params.tmp1} {params.tmp2} {params.n_subset} {params.seed}

        JULIA_PROJECT="{workflow.basedir}"
        julia --project="$JULIA_PROJECT" \
        -e 'using TreeKnit; TreeKnit.treeknit(ARGS[1], ARGS[2]; outdir=ARGS[3], no_likelihood=true)' \
        {params.tmp1} {params.tmp2} {output.outdir}

        python {params.summarize} {input.t1} {input.t2} {output.outdir} {params.n_subset}
        rm -f {params.tmp1} {params.tmp2}
        """



rule treeknit_all:
    input:
        individual=expand("results/treeknit/{seg1}_{seg2}/results_summary.txt",zip,
            seg1=[a for a,b in PAIRS],
            seg2=[b for a,b in PAIRS],
        )
    output:
        summary=f"results/compare_trees_treeknit_{TAG}.csv"

    shell: 
        r"""
        set -euo pipefail
        rm -f {output.summary}
        for f in {input.individual}; do
            python scripts/treeknit_summary_to_csv.py "$f" "{output.summary}"
        done        
        """



 rule treeknit_visualisation:
    input:
        summary="results/compare_trees_treeknit_{TAG}.csv"
    output:
        done="results/treeknit_plots/.done"
    shell:
        r"""
        python scripts/treeknit/treeknit_visualisation.py {input.summary}
        touch {output.done}
        """