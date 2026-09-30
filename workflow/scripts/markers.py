#!/usr/bin/env python
"""Per-cluster marker gene identification.

Differential expression with scanpy's rank_genes_groups on the
log-normalised expression (.raw). Exports the full statistics table
(names, scores, log fold-changes, p-values, adjusted p-values) and a
top-N table per cluster, plus a dotplot of the top markers.
"""

import argparse
import logging
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import scanpy as sc

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
log = logging.getLogger(__name__)


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--input", required=True, help="Annotated .h5ad file")
    p.add_argument("--tables-dir", required=True, help="Directory for output tables")
    p.add_argument("--figures-dir", required=True, help="Directory for figures")
    p.add_argument("--groupby", default="leiden", help=".obs column to compare")
    p.add_argument("--method", default="wilcoxon",
                   choices=["wilcoxon", "t-test", "t-test_overestim_var", "logreg"])
    p.add_argument("--n-genes", type=int, default=25,
                   help="Top markers per cluster to export")
    p.add_argument("--dotplot-n-genes", type=int, default=5,
                   help="Markers per cluster shown on the dotplot")
    p.add_argument("--dpi", type=int, default=300)
    return p.parse_args()


def main():
    args = parse_args()
    figdir = Path(args.figures_dir)
    tabdir = Path(args.tables_dir)
    figdir.mkdir(parents=True, exist_ok=True)
    tabdir.mkdir(parents=True, exist_ok=True)
    sc.settings.set_figure_params(dpi=80, dpi_save=args.dpi, frameon=False)

    log.info("Loading %s", args.input)
    adata = sc.read_h5ad(args.input)

    log.info("rank_genes_groups by '%s' (%s)", args.groupby, args.method)
    sc.tl.rank_genes_groups(adata, groupby=args.groupby, method=args.method,
                            use_raw=True)

    # Full table: names, scores, logfoldchanges, pvals, pvals_adj per group
    df = sc.get.rank_genes_groups_df(adata, group=None)
    df.to_csv(tabdir / "marker_genes_all.csv", index=False)

    top = df.groupby("group").head(args.n_genes)
    top_path = tabdir / f"marker_genes_top{args.n_genes}.csv"
    top.to_csv(top_path, index=False)
    log.info("Wrote %s (%d rows) and full table (%d rows)",
             top_path, len(top), len(df))

    sc.pl.rank_genes_groups_dotplot(adata, n_genes=args.dotplot_n_genes,
                                    standard_scale="var", show=False)
    plt.savefig(figdir / "marker_dotplot.png", dpi=args.dpi, bbox_inches="tight")
    plt.close()


if __name__ == "__main__":
    main()
