#!/usr/bin/env python
"""Doublet detection, Leiden clustering and UMAP embedding.

Runs Scrublet on the raw counts layer, removes predicted doublets,
then builds the neighbourhood graph on the PCA space computed during
preprocessing, clusters with Leiden and embeds with UMAP.
"""

import argparse
import json
import logging
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import scanpy as sc
import scrublet as scr

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
log = logging.getLogger(__name__)


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--input", required=True, help="Preprocessed .h5ad file")
    p.add_argument("--output", required=True, help="Output clustered .h5ad file")
    p.add_argument("--summary", required=True, help="Output JSON summary file")
    p.add_argument("--figures-dir", required=True, help="Directory for figures")
    p.add_argument("--expected-doublet-rate", type=float, default=0.06)
    p.add_argument("--n-neighbors", type=int, default=10)
    p.add_argument("--n-pcs", type=int, default=40)
    p.add_argument("--leiden-resolution", type=float, default=0.5)
    p.add_argument("--dpi", type=int, default=300)
    p.add_argument("--seed", type=int, default=42)
    return p.parse_args()


def main():
    args = parse_args()
    figdir = Path(args.figures_dir)
    figdir.mkdir(parents=True, exist_ok=True)
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    Path(args.summary).parent.mkdir(parents=True, exist_ok=True)
    sc.settings.set_figure_params(dpi=80, dpi_save=args.dpi, frameon=False)

    log.info("Loading %s", args.input)
    adata = sc.read_h5ad(args.input)

    # --- Doublet detection on raw counts ---
    log.info("Running Scrublet (expected doublet rate %.3f)", args.expected_doublet_rate)
    scrub = scr.Scrublet(
        adata.layers["counts"],
        expected_doublet_rate=args.expected_doublet_rate,
        random_state=args.seed,
    )
    doublet_scores, predicted_doublets = scrub.scrub_doublets(verbose=False)
    if predicted_doublets is None:
        # Automatic thresholding failed; fall back to the expected-rate quantile
        cutoff = np.quantile(doublet_scores, 1 - args.expected_doublet_rate)
        predicted_doublets = doublet_scores > cutoff
        log.warning("Scrublet threshold failed; using score quantile %.4f", cutoff)

    adata.obs["doublet_score"] = doublet_scores
    adata.obs["predicted_doublet"] = predicted_doublets

    fig, _ = scrub.plot_histogram()
    fig.savefig(figdir / "doublet_histogram.png", dpi=args.dpi, bbox_inches="tight")
    plt.close(fig)

    n_doublets = int(predicted_doublets.sum())
    log.info("Detected %d / %d doublets (%.1f%%)",
             n_doublets, adata.n_obs, 100 * n_doublets / adata.n_obs)
    adata = adata[~adata.obs["predicted_doublet"]].copy()

    # --- Neighbourhood graph, Leiden clustering, UMAP ---
    log.info("Computing neighbourhood graph (%d neighbors, %d PCs)",
             args.n_neighbors, args.n_pcs)
    sc.pp.neighbors(adata, n_neighbors=args.n_neighbors, n_pcs=args.n_pcs,
                    random_state=args.seed)
    sc.tl.leiden(adata, resolution=args.leiden_resolution,
                 random_state=args.seed, flavor="leidenalg", key_added="leiden")
    sc.tl.umap(adata, random_state=args.seed)

    n_clusters = adata.obs["leiden"].nunique()
    log.info("Leiden clustering: %d clusters", n_clusters)

    sc.pl.umap(adata, color="leiden", legend_loc="on data", show=False)
    plt.savefig(figdir / "umap_clusters.png", dpi=args.dpi, bbox_inches="tight")
    plt.close()

    adata.write_h5ad(args.output)
    log.info("Wrote %s", args.output)

    summary = {
        "n_doublets_detected": n_doublets,
        "n_cells_after_doublet_removal": int(adata.n_obs),
        "n_clusters": int(n_clusters),
        "cluster_sizes": adata.obs["leiden"].value_counts().sort_index().to_dict(),
        "params": vars(args),
    }
    Path(args.summary).write_text(json.dumps(summary, indent=2) + "\n")
    log.info("Wrote %s", args.summary)


if __name__ == "__main__":
    main()
