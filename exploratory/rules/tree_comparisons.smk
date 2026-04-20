"""
This part of the workflow to compare trees between different regions of HBV

REQUIRED INPUTS:
-  trees: expand(f"data/regions/{{name}}/{{name}}.tree.nwk", name=REGION_NAMES),

OUTPUTS:
-   matrix = f"{PW_RESULTS}/corr_matrix/pw_tip_distance_correlation_matrix.patristic.png",
-   grid = f"{PW_RESULTS}/corr_grid/grid_patristic_scatter.pdf",
-   f"{RF_RESULTS}/compare_trees_RF.tsv"
-   grid= f"{RF_RESULTS}/RF.pdf"
-   done=f"{TREEKNIT_RESULTS}/plots/.done"
"""

TAG = make_tag(REGION_NAMES)
PW_RESULTS = f"{RESULTS}/pw_distances"
RF_RESULTS = f"{RESULTS}/rf"
TREEKNIT_RESULTS = f"{RESULTS}/treeknit"

rule compare_pairwise_distances:
    input:
        trees = expand(f"{OUTDIR}/{{name}}/{{name}}.tree.nwk", name=REGION_NAMES),
    output:
        matrix = f"{PW_RESULTS}/corr_matrix/pw_tip_distance_correlation_matrix.patristic.png",
        grid   = f"{PW_RESULTS}/corr_grid/grid_patristic_scatter.pdf",
    params:
        n_pairs = config.get("n_subsamples_pairwise_distance_comparison", 100),
        outdir  = PW_RESULTS,
    shell:
        r"""
        set -euo pipefail
        python scripts/pairwise_tip_distance.py {input.trees} {params.n_pairs} {params.outdir}
        """

rule compare_trees_RF_pair:
    input:
        t1 = f"{OUTDIR}" + "/{seg1}/{seg1}.tree.nwk",
        t2 = f"{OUTDIR}" + "/{seg2}/{seg2}.tree.nwk",
    output:
        f"{RF_RESULTS}/pairwise/{{seg1}}_{{seg2}}.tsv"
    shell:
        r"""
        set -euo pipefail
        python scripts/compare_region_trees.py {input.t1} {input.t2} {output}
        """

rule compare_trees_RF_all:
    input:
        expand(f"{RF_RESULTS}/pairwise/{{seg1}}_{{seg2}}.tsv", zip,
               seg1=[a for a,b in PAIRS],
               seg2=[b for a,b in PAIRS])
    output:
        f"{RF_RESULTS}/compare_trees_RF.tsv"

    shell:
        r"""
        set -euo pipefail
        head -n 1 {input[0]} > {output}
        for f in {input}; do
            tail -n +2 "$f" >> {output}
        done
        """

rule display_RF:
    input:
        table=f"{RF_RESULTS}/compare_trees_RF.tsv"
    output:
        grid= f"{RF_RESULTS}/RF.pdf"
    shell:
        r"""
        python scripts/visualise_rf.py {input.table} {output.grid}
        """

rule tree_knit_pair:
    input:
        t1 = "data/regions/{seg1}/{seg1}.tree.nwk",
        t2 = "data/regions/{seg2}/{seg2}.tree.nwk",
    output:
        outdir = directory(f"{TREEKNIT_RESULTS}/runs/{{seg1}}_{{seg2}}")
    params:
        tmp1 = f"{TREEKNIT_RESULTS}/runs/{{seg1}}_{{seg2}}/input1.nwk",
        tmp2 = f"{TREEKNIT_RESULTS}/runs/{{seg1}}_{{seg2}}/input2.nwk",
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
        individual=expand(f"{TREEKNIT_RESULTS}/runs/{{seg1}}_{{seg2}}/results_summary.txt",zip,
            seg1=[a for a,b in PAIRS],
            seg2=[b for a,b in PAIRS],
        )
    output:
        summary=f"{TREEKNIT_RESULTS}/compare_trees_treeknit.csv"

    shell:
        r"""
        set -euo pipefail
        rm -f {output.summary}
        for f in {input.individual}; do
            python scripts/treeknit/treeknit_summary_to_csv.py "$f" "{output.summary}"
        done
        """

rule treeknit_visualisation:
    input:
        summary=f"{TREEKNIT_RESULTS}/compare_trees_treeknit.csv"
    output:
        done=f"{TREEKNIT_RESULTS}/plots/.done"
    params:
        outdir=f"{TREEKNIT_RESULTS}/plots",
    shell:
        r"""
        python scripts/treeknit/treeknit_visualisation.py {input.summary} --outdir {params.outdir}
        touch {output.done}
        """
