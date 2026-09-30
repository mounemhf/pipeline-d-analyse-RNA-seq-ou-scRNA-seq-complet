#!/usr/bin/env python
"""Automated cell-type annotation with CellTypist.

Annotates cells with a reference immune model on the log-normalised
expression stored in .raw, optionally smoothed by majority voting over
the Leiden clusters. Exports cluster-level annotation tables and a UMAP
coloured by predicted cell type.
"""

import argparse
import logging
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc
import celltypist

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
log = logging.getLogger(__name__)


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--input", required=True, help="Clustered .h5ad file")
    p.add_argument("--output", required=True, help="Output annotated .h5ad file")
    p.add_argument("--tables-dir", required=True, help="Directory for output tables")
    p.add_argument("--figures-dir", required=True, help="Directory for figures")
    p.add_argument("--model", default="Immune_All_High.pkl",
                   help="CellTypist model name (auto-downloaded on first use)")
    p.add_argument("--majority-voting", action="store_true",
                   help="Smooth predictions by majority voting over Leiden clusters")
    p.add_argument("--dpi", type=int, default=300)
    return p.parse_args()


def main():
    args = parse_args()
    figdir = Path(args.figures_dir)
    tabdir = Path(args.tables_dir)
    figdir.mkdir(parents=True, exist_ok=True)
    tabdir.mkdir(parents=True, exist_ok=True)
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    sc.settings.set_figure_params(dpi=80, dpi_save=args.dpi, frameon=False)

    log.info("Loading %s", args.input)
    adata = sc.read_h5ad(args.input)

    # CellTypist expects log1p-normalised expression (10k per cell), kept in .raw
    log.info("Annotating with CellTypist model %s", args.model)
    adata_log = adata.raw.to_adata()
    predictions = celltypist.annotate(
        adata_log,
        model=args.model,
        majority_voting=args.majority_voting,
        over_clustering="leiden" if args.majority_voting else None,
    )
    # to_adata() does not take the target AnnData; copy the prediction
    # columns explicitly. conf_score is not stored in predicted_labels:
    # like celltypist's own to_adata(), derive it from probability_matrix.
    label_col = "majority_voting" if args.majority_voting else "predicted_labels"
    insert_cols = ["predicted_labels"]
    if args.majority_voting:
        insert_cols.append("majority_voting")
    adata.obs[insert_cols] = predictions.predicted_labels[insert_cols]
    adata.obs["conf_score"] = predictions.probability_matrix.max(axis=1).values

    adata.obs["cell_type"] = adata.obs[label_col]
    log.info("Predicted cell types: %s",
             adata.obs["cell_type"].value_counts().to_dict())

    # Cluster x cell-type contingency table and per-cluster consensus label
    crosstab = pd.crosstab(adata.obs["leiden"], adata.obs["cell_type"])
    crosstab.to_csv(tabdir / "celltype_by_cluster.csv")
    cluster_annotation = pd.DataFrame({
        "leiden": crosstab.index,
        "cell_type": crosstab.idxmax(axis=1).values,
        "n_cells": adata.obs["leiden"].value_counts().sort_index().values,
    })
    cluster_annotation.to_csv(tabdir / "cluster_annotation.csv", index=False)
    log.info("Wrote annotation tables to %s", tabdir)

    sc.pl.umap(adata, color="cell_type", show=False)
    plt.savefig(figdir / "umap_celltypes.png", dpi=args.dpi, bbox_inches="tight")
    plt.close()

    adata.write_h5ad(args.output)
    log.info("Wrote %s", args.output)


if __name__ == "__main__":
    main()
