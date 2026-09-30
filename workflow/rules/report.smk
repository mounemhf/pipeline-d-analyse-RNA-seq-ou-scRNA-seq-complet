# Final Markdown report assembly.

REPORT_DIR = config["paths"]["report"]


rule report:
    """Assemble the final report with figures and biological interpretation."""
    input:
        preprocess_summary=f"{TABLES}/preprocess_summary.json",
        cluster_summary=f"{TABLES}/cluster_summary.json",
        cluster_annotation=f"{TABLES}/cluster_annotation.csv",
        celltype_by_cluster=f"{TABLES}/celltype_by_cluster.csv",
        markers_top=f"{TABLES}/marker_genes_top{config['markers']['n_genes']}.csv",
    output:
        f"{REPORT_DIR}/report.md",
    params:
        figures_dir=FIGURES,
        dataset="10x Genomics PBMC 3k (v3 chemistry)",
        model=config["annotation"]["celltypist_model"],
    conda:
        "../envs/report.yml"
    log:
        f"{LOG_DIR}/report.log"
    shell:
        """
        python workflow/scripts/generate_report.py \
            --preprocess-summary {input.preprocess_summary} \
            --cluster-summary {input.cluster_summary} \
            --cluster-annotation {input.cluster_annotation} \
            --celltype-by-cluster {input.celltype_by_cluster} \
            --markers-top {input.markers_top} \
            --figures-dir {params.figures_dir} \
            --output {output} \
            --dataset "{params.dataset}" \
            --model {params.model} \
            > {log} 2>&1
        """
