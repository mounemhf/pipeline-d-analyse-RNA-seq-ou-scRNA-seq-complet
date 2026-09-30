# Cell-type annotation (CellTypist) and per-cluster marker genes.

MAJORITY_VOTING_FLAG = "--majority-voting" if config["annotation"]["majority_voting"] else ""


rule annotate:
    """Annotate cell types with a CellTypist immune reference model."""
    input:
        h5ad=f"{RESULTS}/processed/pbmc3k_clustered.h5ad",
    output:
        h5ad=f"{RESULTS}/processed/pbmc3k_annotated.h5ad",
        by_cluster=f"{TABLES}/celltype_by_cluster.csv",
        cluster_labels=f"{TABLES}/cluster_annotation.csv",
    params:
        figures_dir=FIGURES,
        tables_dir=TABLES,
        model=config["annotation"]["celltypist_model"],
        majority_voting=MAJORITY_VOTING_FLAG,
        dpi=config["figures"]["dpi"],
    conda:
        "../envs/annotate.yml"
    log:
        f"{LOG_DIR}/annotate.log"
    threads:
        config["threads"]
    shell:
        """
        python workflow/scripts/annotate.py \
            --input {input.h5ad} \
            --output {output.h5ad} \
            --tables-dir {params.tables_dir} \
            --figures-dir {params.figures_dir} \
            --model {params.model} \
            {params.majority_voting} \
            --dpi {params.dpi} \
            > {log} 2>&1
        """


rule markers:
    """Identify marker genes per Leiden cluster."""
    input:
        h5ad=f"{RESULTS}/processed/pbmc3k_annotated.h5ad",
    output:
        all=f"{TABLES}/marker_genes_all.csv",
        top=f"{TABLES}/marker_genes_top{config['markers']['n_genes']}.csv",
    params:
        figures_dir=FIGURES,
        tables_dir=TABLES,
        n_genes=config["markers"]["n_genes"],
        method=config["markers"]["method"],
        dpi=config["figures"]["dpi"],
    conda:
        "../envs/markers.yml"
    log:
        f"{LOG_DIR}/markers.log"
    threads:
        config["threads"]
    shell:
        """
        python workflow/scripts/markers.py \
            --input {input.h5ad} \
            --tables-dir {params.tables_dir} \
            --figures-dir {params.figures_dir} \
            --n-genes {params.n_genes} \
            --method {params.method} \
            --dpi {params.dpi} \
            > {log} 2>&1
        """
