# Doublet detection (Scrublet), Leiden clustering and UMAP embedding.


rule cluster:
    """Detect doublets, then cluster and embed cells."""
    input:
        h5ad=f"{RESULTS}/processed/pbmc3k_preprocessed.h5ad",
    output:
        h5ad=f"{RESULTS}/processed/pbmc3k_clustered.h5ad",
        summary=f"{TABLES}/cluster_summary.json",
    params:
        figures_dir=FIGURES,
        expected_doublet_rate=config["doublets"]["expected_doublet_rate"],
        n_neighbors=config["preprocessing"]["n_neighbors"],
        n_pcs=config["preprocessing"]["n_pcs"],
        leiden_resolution=config["clustering"]["leiden_resolution"],
        dpi=config["figures"]["dpi"],
        seed=config["seed"],
    conda:
        "../envs/cluster.yml"
    log:
        f"{LOG_DIR}/cluster.log"
    threads:
        config["threads"]
    shell:
        """
        python workflow/scripts/cluster.py \
            --input {input.h5ad} \
            --output {output.h5ad} \
            --summary {output.summary} \
            --figures-dir {params.figures_dir} \
            --expected-doublet-rate {params.expected_doublet_rate} \
            --n-neighbors {params.n_neighbors} \
            --n-pcs {params.n_pcs} \
            --leiden-resolution {params.leiden_resolution} \
            --dpi {params.dpi} \
            --seed {params.seed} \
            > {log} 2>&1
        """
