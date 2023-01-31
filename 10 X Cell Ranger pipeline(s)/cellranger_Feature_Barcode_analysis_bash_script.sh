#!/bin/bash

# This bash script pipeline outputs a unified feature-barcode matrix that contains gene expression counts alongside Feature Barcode counts for # each cell barcode. 
# To enable Feature Barcode analysis, cellranger count needs two new inputs:
# 1. Feature Reference CSV is passed to cellranger count with the --feature-ref flag and declares the set of Feature Barcode reagents in use in # the experiment. 
# 2. Libraries CSV is passed to cellranger count with the --libraries flag, and declares the FASTQ files and library type for each input 
# dataset. In a typical Feature Barcode analysis there will be two input libraries: one for the normal single-cell gene expression reads, and # one for the Feature Barcode reads.

# Other important Arguments
# --id: Required. A unique run ID string (e.g., Sample1_Fastq).
# --run: Required. The path of Illumina BCL run folder.
# --csv	Optional. Path to a simple CSV with lane, sample, and index columns, which describe the way to demultiplex the flow cell. The index 
# column should contain a 10x Genomics sample dual-index name (e.g., SI-TT-A12).
# https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/feature-bc-analysis


# Note: this requires high resources (ram and processing core).
# Run cellranger count to generate single cell feature counts.
# Required. Sample name as specified in the sample sheet supplied to cellranger mkfastq.

cellranger count --id=Sample1_Feature_Barcode_outPut \
                   --libraries=/Functional_Heterogeneity/Feature_Barcode/library.csv \
                   --transcriptome=/opt/refdata-cellranger-GRCh38-3.0.0 \
		   --feature-ref=/Functional_Heterogeneity/Feature_Barcode/feature_ref.csv 
