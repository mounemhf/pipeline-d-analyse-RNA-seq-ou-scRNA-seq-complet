#!/usr/bin/env python
"""Preprocess the 10x PBMC 3k filtered gene-barcode matrix.

Steps: QC metric computation, cell/gene filtering, normalisation,
log1p transformation, highly variable gene selection, scaling and PCA.
Saves the processed AnnData object, QC figures and a JSON summary.
"""

import argparse
import json
import logging
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import scanpy as sc

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
log = logging.getLogger(__name__)

QC_VARS = ["n_genes_by_counts", "total_counts", "pct_counts_mt"]


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--input", required=True,
                   help="10x filtered matrix directory (matrix.mtx, barcodes.tsv, genes.tsv)")
    p.add_argument("--output", required=True, help="Output .h5ad file")
    p.add_argument("--summary", required=True, help="Output JSON summary file")
    p.add_argument("--figures-dir", required=True, help="Directory for QC figures")
    p.add_argument("--min-genes", type=int, default=200)
    p.add_argument("--min-cells", type=int, default=3)
    p.add_argument("--max-counts", type=float, default=15000)
    p.add_argument("--max-pct-mito", type=float, default=5.0)
    p.add_argument("--mito-prefix", default="MT-")
    p.add_argument("--target-sum", type=float, default=1e4)
    p.add_argument("--n-top-hvg", type=int, default=2000)
    p.add_argument("--n-pcs", type=int, default=40)
    p.add_argument("--dpi", type=int, default=300)
    p.add_argument("--seed", type=int, default=42)
    return p.parse_args()


def save_current_figure(path: Path, dpi: int) -> None:
    plt.savefig(path, dpi=dpi, bbox_inches="tight")
    plt.close()
    log.info("Wrote figure %s", path)


def main():
    args = parse_args()
    figdir = Path(args.figures_dir)
    figdir.mkdir(parents=True, exist_ok=True)
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    Path(args.summary).parent.mkdir(parents=True, exist_ok=True)
    sc.settings.set_figure_params(dpi=80, dpi_save=args.dpi, frameon=False)

    log.info("Loading 10x matrix from %s", args.input)
    adata = sc.read_10x_mtx(args.input, var_names="gene_symbols", cache=False)
    adata.var_names_make_unique()
    n_cells_raw, n_genes_raw = adata.shape
    log.info("Raw data: %d cells x %d genes", n_cells_raw, n_genes_raw)

    # Per-cell QC metrics, including mitochondrial fraction
    adata.var["mt"] = adata.var_names.str.startswith(args.mito_prefix)
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None,
                               log1p=False, inplace=True)

    sc.pl.violin(adata, QC_VARS, jitter=0.4, multi_panel=True, show=False)
    save_current_figure(figdir / "qc_violins_prefilter.png", args.dpi)

    # Cell and gene filtering
    sc.pp.filter_cells(adata, min_genes=args.min_genes)
    sc.pp.filter_genes(adata, min_cells=args.min_cells)
    adata = adata[adata.obs["pct_counts_mt"] < args.max_pct_mito].copy()
    adata = adata[adata.obs["total_counts"] < args.max_counts].copy()
    n_cells_filt, n_genes_filt = adata.shape
    log.info("After filtering: %d cells x %d genes", n_cells_filt, n_genes_filt)

    sc.pl.violin(adata, QC_VARS, jitter=0.4, multi_panel=True, show=False)
    save_current_figure(figdir / "qc_violins_postfilter.png", args.dpi)

    # Keep raw counts for doublet detection (Scrublet) downstream
    adata.layers["counts"] = adata.X.copy()

    # Normalisation; keep a log-normalised copy in .raw for marker analysis
    sc.pp.normalize_total(adata, target_sum=args.target_sum)
    sc.pp.log1p(adata)
    adata.raw = adata

    # Highly variable genes, then PCA on HVGs only
    sc.pp.highly_variable_genes(adata, n_top_genes=args.n_top_hvg, flavor="seurat")
    sc.pl.highly_variable_genes(adata, show=False)
    save_current_figure(figdir / "highly_variable_genes.png", args.dpi)

    sc.pp.scale(adata, max_value=10)
    sc.pp.pca(adata, n_comps=args.n_pcs, mask_var="highly_variable",
              svd_solver="arpack", random_state=args.seed)
    sc.pl.pca_variance_ratio(adata, n_pcs=args.n_pcs, log=True, show=False)
    save_current_figure(figdir / "pca_variance_ratio.png", args.dpi)

    adata.write_h5ad(args.output)
    log.info("Wrote %s", args.output)

    summary = {
        "n_cells_raw": int(n_cells_raw),
        "n_genes_raw": int(n_genes_raw),
        "n_cells_filtered": int(n_cells_filt),
        "n_genes_filtered": int(n_genes_filt),
        "median_genes_per_cell": float(adata.obs["n_genes_by_counts"].median()),
        "median_counts_per_cell": float(adata.obs["total_counts"].median()),
        "params": vars(args),
    }
    Path(args.summary).write_text(json.dumps(summary, indent=2) + "\n")
    log.info("Wrote %s", args.summary)


if __name__ == "__main__":
    main()
