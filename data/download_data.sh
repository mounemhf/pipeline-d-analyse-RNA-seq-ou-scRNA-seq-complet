#!/usr/bin/env bash
# Download the 10x Genomics PBMC 3k dataset (v3 chemistry).
#
# Two artefacts are available:
#   - filtered gene-barcode matrix (~28 MB): entry point of the analysis
#   - raw FASTQ files (~4.5 GB): only used for FastQC/MultiQC read-level QC
#
# The pipeline itself consumes the filtered matrix; read alignment and UMI
# counting (Cell Ranger) are out of scope and documented in the README.
#
# Usage:
#   bash data/download_data.sh                  # matrix only (default)
#   bash data/download_data.sh --with-fastq     # matrix + FASTQ
#   bash data/download_data.sh --fastq-only     # FASTQ only
#   bash data/download_data.sh --matrix-only    # matrix only (explicit)
#   bash data/download_data.sh --force          # re-download even if present

set -euo pipefail

MATRIX_URL="https://cf.10xgenomics.com/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz"
FASTQ_URL="https://cf.10xgenomics.com/samples/cell/pbmc3k/pbmc3k_fastqs.tar"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RAW_DIR="${SCRIPT_DIR}/raw"
MATRIX_DIR="${RAW_DIR}/filtered_gene_bc_matrices/hg19"
FASTQ_DIR="${RAW_DIR}/pbmc3k_fastqs"

WITH_MATRIX=1
WITH_FASTQ=0
FORCE=0

usage() {
    sed -n '2,16p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'
    exit "${1:-0}"
}

for arg in "$@"; do
    case "${arg}" in
        --with-fastq)  WITH_FASTQ=1 ;;
        --fastq-only)  WITH_MATRIX=0; WITH_FASTQ=1 ;;
        --matrix-only) WITH_MATRIX=1; WITH_FASTQ=0 ;;
        --force)       FORCE=1 ;;
        -h|--help)     usage 0 ;;
        *) echo "Unknown option: ${arg}" >&2; usage 1 ;;
    esac
done

mkdir -p "${RAW_DIR}"

download() {
    # $1 = URL, $2 = destination file
    echo "[download] $1"
    curl -fSL --retry 3 --retry-delay 5 -o "$2" "$1"
}

fetch_matrix() {
    if [[ -f "${MATRIX_DIR}/matrix.mtx" && "${FORCE}" -eq 0 ]]; then
        echo "[skip] filtered matrix already present in ${MATRIX_DIR}"
        return
    fi
    local archive="${RAW_DIR}/pbmc3k_filtered_gene_bc_matrices.tar.gz"
    download "${MATRIX_URL}" "${archive}"
    tar -xzf "${archive}" -C "${RAW_DIR}"
    rm -f "${archive}"
    for f in matrix.mtx barcodes.tsv genes.tsv; do
        if [[ ! -f "${MATRIX_DIR}/${f}" ]]; then
            echo "[error] expected file missing after extraction: ${MATRIX_DIR}/${f}" >&2
            exit 1
        fi
    done
    echo "[ok] filtered matrix ready in ${MATRIX_DIR}"
}

fetch_fastq() {
    if compgen -G "${FASTQ_DIR}/*.fastq.gz" > /dev/null && [[ "${FORCE}" -eq 0 ]]; then
        echo "[skip] FASTQ files already present in ${FASTQ_DIR}"
        return
    fi
    echo "[warn] FASTQ archive is ~4.5 GB; this may take a while."
    local archive="${RAW_DIR}/pbmc3k_fastqs.tar"
    download "${FASTQ_URL}" "${archive}"
    tar -xf "${archive}" -C "${RAW_DIR}"
    rm -f "${archive}"
    if ! compgen -G "${FASTQ_DIR}/*.fastq.gz" > /dev/null; then
        echo "[error] no .fastq.gz files found after extraction in ${FASTQ_DIR}" >&2
        exit 1
    fi
    echo "[ok] FASTQ files ready in ${FASTQ_DIR}"
}

[[ "${WITH_MATRIX}" -eq 1 ]] && fetch_matrix
[[ "${WITH_FASTQ}" -eq 1 ]] && fetch_fastq

echo "[done] data available under ${RAW_DIR}"
