#!/bin/bash
# Bash script to convert base call files (BCLs) for each flow cell directory into FASTQ files. 
# Important Arguments
# --id: Required. A unique run ID string (e.g., Sample1_Fastq).
# --run: Required. The path of Illumina BCL run folder.
# --csv	Optional. Path to a simple CSV with lane, sample, and index columns, which describe the way to demultiplex the flow cell. The index 
# column should contain a 10x Genomics sample dual-index name (e.g., SI-TT-A12).
# https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/mkfastq

# Note: this requires high resources (ram and processing core),
# Run cellranger mkfastq for BCL to FASTQ conversion.

cellranger mkfastq --id=Sample1_Fastq \
                     --run=/path/to/BCL \
                     --csv=cellranger-bcl-to-fastq-Sample1.csv

