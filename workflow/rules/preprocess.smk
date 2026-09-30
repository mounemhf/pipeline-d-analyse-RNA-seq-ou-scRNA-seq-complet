# Cell-level preprocessing of the filtered gene-barcode matrix:
# QC filtering, normalisation, highly variable genes, PCA.

RESULTS = config["paths"]["results"]
FIGURES = config["paths"]["figures"]
TABLES = config["paths"]["tables"]


rule preprocess:
    """Filter, normalise and embed the count matrix (scanpy)."""
    input:
        matrix=f"{MATRIX_DIR}/matrix.mtx",
        barcodes=f"{MATRIX_DIR}/barcodes.tsv",
        genes=f"{MATRIX_DIR}/genes.tsv",
    output:
        h5ad=f"{RESULTS}/processed/pbmc3k_preprocessed.h5ad",
        summary=f"{TABLES}/preprocess_summary.json",
    params:
        figures_dir=FIGURES,
        min_genes=config["qc"]["min_genes_per_cell"],
        min_cells=config["qc"]["min_cells_per_gene"],
        max_counts=config["qc"]["max_counts"],
        max_pct_mito=config["qc"]["max_pct_mito"],
        mito_prefix=config["qc"]["mito_prefix"],
        target_sum=config["preprocessing"]["target_sum"],
        n_top_hvg=config["preprocessing"]["n_top_hvg"],
        n_pcs=config["preprocessing"]["n_pcs"],
        dpi=config["figures"]["dpi"],
        seed=config["seed"],
    conda:
        "../envs/preprocess.yml"
    log:
        f"{LOG_DIR}/preprocess.log"
    threads:
        config["threads"]
    shell:
        """
        python workflow/scripts/preprocess.py \
            --input {MATRIX_DIR} \
            --output {output.h5ad} \
            --summary {output.summary} \
            --figures-dir {params.figures_dir} \
            --min-genes {params.min_genes} \
            --min-cells {params.min_cells} \
            --max-counts {params.max_counts} \
            --max-pct-mito {params.max_pct_mito} \
            --mito-prefix {params.mito_prefix} \
            --target-sum {params.target_sum} \
            --n-top-hvg {params.n_top_hvg} \
            --n-pcs {params.n_pcs} \
            --dpi {params.dpi} \
            --seed {params.seed} \
            > {log} 2>&1
        """
