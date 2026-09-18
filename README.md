# TE5primeIsoforms

Analysis code and processed data for TE-associated 5′ isoforms and their potential coding products in cancer cell lines.

## Overview

This repository contains the analysis scripts used to identify and evaluate transposable element (TE)-associated 5′ transcript isoforms from public PacBio Iso-Seq data across cancer cell lines.

The workflow includes long-read transcript processing and quality control, TE annotation, first-exon event identification, transcription start site support assessment, matched short-read RNA-seq splice-junction support, manual structural curation, ORF and peptide analysis, HLA-I binding and presentation prediction, reference-sequence comparison, and reanalysis of public mass spectrometry data.

## Workflow

![Analysis workflow](figures/workflow.png)

The analysis starts from public PacBio Iso-Seq datasets and integrates transcript structure annotation, repeat annotation, transcription start site evidence, matched short-read RNA-seq, coding-sequence analysis, peptide prediction, and complementary mass spectrometry evidence.

## Repository structure

```text
TE5primeIsoforms/
├── scripts/
│   ├── 01_pbmm2.sh
│   ├── 02_collapse.sh
│   ├── 03_sqanti3.sh
│   ├── ...
│   └── 18_dedup_orf_vs_uniprot.py
│
├── figures/
│   └── workflow.png
│
└── README.md
