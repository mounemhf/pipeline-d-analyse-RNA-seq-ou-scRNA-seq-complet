# Data

Raw and processed data are **not** tracked by git (see `.gitignore`).

Download the 10x Genomics PBMC 3k dataset (v3 chemistry) with:

```bash
bash data/download_data.sh                # filtered matrix only (~28 MB)
bash data/download_data.sh --with-fastq   # also raw FASTQ (~17 GB, for FastQC/MultiQC)
```

The analysis starts from the official Cell Ranger filtered gene-barcode
matrix. The raw FASTQ files are used exclusively for read-level QC
(FastQC/MultiQC); alignment and UMI counting are out of scope.
