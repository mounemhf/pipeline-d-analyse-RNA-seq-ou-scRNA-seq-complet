# Read-level quality control: FastQC per FASTQ file, then MultiQC aggregation.
# Runs only when the raw FASTQ files have been downloaded
# (bash data/download_data.sh --with-fastq); otherwise the pipeline skips
# this module (see rule all in the Snakefile).

import glob
import os

FASTQ_DIR = config["paths"]["fastq_dir"]
QC_DIR = f"{config['paths']['results']}/qc"

FASTQ_FILES = (
    sorted(glob.glob(os.path.join(FASTQ_DIR, "*.fastq.gz")))
    if os.path.isdir(FASTQ_DIR)
    else []
)


rule fastqc:
    """Per-file read quality control with FastQC."""
    input:
        FASTQ_FILES,
    output:
        directory(f"{QC_DIR}/fastqc"),
    conda:
        "../envs/fastqc.yml"
    log:
        f"{LOG_DIR}/fastqc.log"
    threads:
        config["threads"]
    shell:
        "mkdir -p {output} && fastqc -t {threads} -o {output} {input} > {log} 2>&1"


rule multiqc:
    """Aggregate all FastQC reports into a single MultiQC report."""
    input:
        f"{QC_DIR}/fastqc",
    output:
        f"{QC_DIR}/multiqc_report.html",
    conda:
        "../envs/multiqc.yml"
    log:
        f"{LOG_DIR}/multiqc.log"
    shell:
        "multiqc --force -o {QC_DIR} -n multiqc_report.html {input} > {log} 2>&1"
