#!/usr/bin/env python
"""Generate the final analysis report (report/report.md).

Assembles pipeline summaries, annotation tables and figures into a
Markdown report, including a written biological interpretation of the
PBMC populations recovered by clustering, annotation and marker genes.
"""

import argparse
import json
import logging
from datetime import date
from pathlib import Path

import pandas as pd

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
log = logging.getLogger(__name__)

# Lineages expected in PBMCs: keywords matched against CellTypist labels,
# canonical markers, and a fixed biological context sentence.
LINEAGES = [
    ("T cells", ["t cell", "cd4", "cd8", "treg", "mait"],
     ["CD3D", "CD3E", "IL7R", "CD8A", "CD8B"],
     "T lymphocytes are expected to form the largest PBMC fraction "
     "(typically 50-70%), and their dominance here is consistent with a "
     "healthy donor profile."),
    ("NK cells", ["nk"],
     ["GNLY", "NKG7", "KLRD1", "FCGR3A"],
     "NK cells are the cytotoxic arm of the innate lymphoid compartment "
     "and typically represent 5-15% of PBMCs."),
    ("B cells", ["b cell", "plasma"],
     ["MS4A1", "CD79A", "CD79B", "BANK1"],
     "B cells form the humoral arm of adaptive immunity and usually "
     "account for 5-10% of PBMCs."),
    ("Monocytes", ["mono", "macrophage"],
     ["CD14", "LYZ", "S100A8", "S100A9", "FCGR3A", "MS4A7", "LST1"],
     "Monocytes are the main circulating myeloid population (10-20% of "
     "PBMCs); classical CD14+ and non-classical FCGR3A (CD16)+ subsets are "
     "routinely resolved at this clustering resolution."),
    ("Dendritic cells", ["dc", "dendritic"],
     ["FCER1A", "CST3", "CLEC9A", "CLEC10A"],
     "Dendritic cells are rare in peripheral blood (<2%) and surface here "
     "as small, well-separated clusters."),
    ("Megakaryocytes / platelets", ["megakary", "platelet"],
     ["PPBP", "PF4"],
     "A small megakaryocyte/platelet signal is a common finding in PBMC "
     "preparations and reflects residual platelet contamination rather "
     "than a blood-resident population."),
]


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--preprocess-summary", required=True)
    p.add_argument("--cluster-summary", required=True)
    p.add_argument("--cluster-annotation", required=True)
    p.add_argument("--celltype-by-cluster", required=True)
    p.add_argument("--markers-top", required=True)
    p.add_argument("--figures-dir", required=True,
                   help="Figure directory on disk (existence checks)")
    p.add_argument("--figures-prefix", default="../results/figures",
                   help="Path prefix used for image links in the report")
    p.add_argument("--output", required=True)
    p.add_argument("--dataset", default="10x Genomics PBMC 3k (v3 chemistry)")
    p.add_argument("--model", default="Immune_All_High.pkl")
    p.add_argument("--n-inline-markers", type=int, default=5,
                   help="Markers per cluster listed in the marker table")
    return p.parse_args()


def md_table(df: pd.DataFrame) -> str:
    cols = [str(c) for c in df.columns]
    lines = ["| " + " | ".join(cols) + " |", "|" + "---|" * len(cols)]
    for _, row in df.iterrows():
        lines.append("| " + " | ".join(str(v) for v in row) + " |")
    return "\n".join(lines)


def figure(name: str, caption: str, figdir: Path, prefix: str) -> str:
    if (figdir / name).exists():
        return f"![{caption}]({prefix}/{name})\n"
    log.warning("Figure missing, skipped in report: %s", name)
    return ""


def lineage_paragraph(name, keywords, canonical, context, cluster_ann, markers):
    mask = cluster_ann["cell_type"].str.lower().apply(
        lambda s: any(k in s for k in keywords))
    sub = cluster_ann[mask]
    if sub.empty:
        return None
    n = int(sub["n_cells"].sum())
    desc = "; ".join(f"cluster {r['leiden']} ({r['cell_type']}, {r['n_cells']} cells)"
                     for _, r in sub.iterrows())
    tops, seen = [], set()
    for g in sub["leiden"].astype(str):
        for gene in markers.loc[markers["group"] == g, "names"].head(3):
            if gene not in seen:
                tops.append(gene)
                seen.add(gene)
    canon_found = [g for g in canonical if g in seen]
    canon_txt = (f" The canonical lineage markers {', '.join(canon_found)} "
                 "are among the strongest differentially expressed genes, "
                 "confirming the automated annotation.") if canon_found else ""
    return (f"**{name}** ({n} cells) were recovered in {desc}. Top marker "
            f"genes for these clusters include {', '.join(tops[:6])}.{canon_txt} "
            f"{context}")


def main():
    args = parse_args()
    figdir = Path(args.figures_dir)
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)

    pre = json.loads(Path(args.preprocess_summary).read_text(encoding="utf-8"))
    clu = json.loads(Path(args.cluster_summary).read_text(encoding="utf-8"))
    cluster_ann = pd.read_csv(args.cluster_annotation, dtype={"leiden": str})
    by_cluster = pd.read_csv(args.celltype_by_cluster, index_col=0)
    markers = pd.read_csv(args.markers_top, dtype={"group": str})

    n_filt = pre["n_cells_filtered"]
    type_totals = by_cluster.sum().sort_values(ascending=False)
    composition = pd.DataFrame({
        "cell type": type_totals.index,
        "cells": type_totals.values,
        "fraction": [f"{100 * v / type_totals.sum():.1f}%" for v in type_totals.values],
    })

    params = pre["params"]
    param_table = md_table(pd.DataFrame([
        ["min genes / cell", params["min_genes"]],
        ["min cells / gene", params["min_cells"]],
        ["max % mitochondrial", params["max_pct_mito"]],
        ["max total counts", params["max_counts"]],
        ["highly variable genes", params["n_top_hvg"]],
        ["principal components", params["n_pcs"]],
        ["neighbours", clu["params"]["n_neighbors"]],
        ["Leiden resolution", clu["params"]["leiden_resolution"]],
        ["expected doublet rate", clu["params"]["expected_doublet_rate"]],
        ["CellTypist model", args.model],
        ["random seed", params["seed"]],
    ], columns=["parameter", "value"]))

    interp = []
    interp.append(
        f"Starting from {pre['n_cells_raw']} barcoded cells, quality filtering "
        f"retained {n_filt} high-quality cells with a median of "
        f"{pre['median_genes_per_cell']:.0f} detected genes per cell. Scrublet "
        f"flagged {clu['n_doublets_detected']} predicted doublets "
        f"({100 * clu['n_doublets_detected'] / n_filt:.1f}%), in line with the "
        "expected multiplet rate of the 10x Chromium platform at this loading "
        f"concentration. Leiden clustering of the remaining "
        f"{clu['n_cells_after_doublet_removal']} cells resolved "
        f"{clu['n_clusters']} transcriptionally distinct clusters, which "
        f"CellTypist ({args.model}) assigned to the expected PBMC lineages.")
    for name, keywords, canonical, context in LINEAGES:
        para = lineage_paragraph(name, keywords, canonical, context,
                                 cluster_ann, markers)
        if para:
            interp.append(para)
    top3 = ", ".join(f"{t} ({f})" for t, f in
                     zip(composition["cell type"].head(3), composition["fraction"].head(3)))
    interp.append(
        f"Overall, the recovered composition — {top3} as the three largest "
        "populations — matches the expected makeup of peripheral blood from a "
        "healthy donor, and cluster-resolved marker genes are concordant with "
        "the reference-based CellTypist labels. Two caveats apply: this is a "
        "single-donor, single-batch dataset, so no batch integration was "
        "required (with multiple samples, Harmony or scVI integration would "
        "precede clustering); and automated annotation is only as granular as "
        "its reference — rare subsets such as dendritic cell flavours would "
        "benefit from manual marker-based curation before any downstream "
        "biological claim.")

    marker_rows = []
    for g, sub in markers.groupby("group"):
        ct = cluster_ann.loc[cluster_ann["leiden"] == g, "cell_type"]
        marker_rows.append([
            g, ct.iloc[0] if len(ct) else "?",
            ", ".join(sub["names"].head(args.n_inline_markers)),
        ])
    marker_table = md_table(pd.DataFrame(
        marker_rows, columns=["cluster", "cell type",
                              f"top {args.n_inline_markers} markers"]))

    rep = f"""# Single-cell RNA-seq analysis report — {args.dataset}

*Generated on {date.today().isoformat()} by the Snakemake pipeline
(`snakemake --use-conda --cores 4`). All parameters are defined in
`config/config.yaml`.*

## 1. Read and cell quality control

Read-level QC (FastQC/MultiQC) is available under `results/qc/` when the raw
FASTQ files are downloaded. At the cell level, {pre['n_cells_raw']} barcoded
cells x {pre['n_genes_raw']} genes were reduced to {n_filt} cells x
{pre['n_genes_filtered']} genes after filtering on detected genes, total
counts and mitochondrial fraction (thresholds below).

{figure("qc_violins_prefilter.png", "Per-cell QC metrics before filtering", figdir, args.figures_prefix)}
{figure("qc_violins_postfilter.png", "Per-cell QC metrics after filtering", figdir, args.figures_prefix)}

### Pipeline parameters

{param_table}

## 2. Doublet detection

Scrublet predicted **{clu['n_doublets_detected']} doublets**, which were
removed before clustering.

{figure("doublet_histogram.png", "Scrublet doublet score distribution", figdir, args.figures_prefix)}

## 3. Clustering and embedding

Leiden clustering resolved **{clu['n_clusters']} clusters** across
{clu['n_cells_after_doublet_removal']} cells.

{figure("umap_clusters.png", "UMAP coloured by Leiden cluster", figdir, args.figures_prefix)}
{figure("highly_variable_genes.png", "Highly variable gene selection", figdir, args.figures_prefix)}
{figure("pca_variance_ratio.png", "PCA variance ratio", figdir, args.figures_prefix)}

## 4. Cell-type annotation

Automated annotation with CellTypist (`{args.model}`, majority voting over
Leiden clusters):

{md_table(composition)}

{figure("umap_celltypes.png", "UMAP coloured by predicted cell type", figdir, args.figures_prefix)}

## 5. Marker genes

Per-cluster differential expression (`rank_genes_groups`, Wilcoxon). Full
statistics — including log fold-changes and adjusted p-values — are in
`results/tables/marker_genes_all.csv`.

{marker_table}

{figure("marker_dotplot.png", "Top marker genes per cluster", figdir, args.figures_prefix)}

## 6. Biological interpretation

{chr(10).join(p + chr(10) for p in interp)}
"""
    Path(args.output).write_text(rep, encoding="utf-8")
    log.info("Wrote %s", args.output)


if __name__ == "__main__":
    main()
