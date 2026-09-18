# Analysis of TE-associated 5′ isoforms in cancer cell lines

## Overview

This code supports the identification, structural curation, and coding-potential analysis of transposable element (TE)-associated 5′ transcript isoforms in cancer cell lines. The main workflow includes data preparation, full-length isoform generation and quality control, TE annotation and first-exon event identification, transcription start site support assessment, matched short-read RNA-seq support assessment, manual structural curation and event consolidation, open reading frame (ORF) and peptide analysis, human leukocyte antigen class I (HLA-I) binding and antigen processing/presentation prediction, normal-tissue and reference-protein background assessment, and reanalysis of public mass spectrometry data.

The code was organized according to the main analysis workflow used in the accompanying manuscript. Representative scripts, filtering rules, and result-generation steps are retained to document the analysis of TE-associated 5′ isoforms from long-read identification and structural validation to coding-potential and peptide-level evaluation.

Large raw sequencing files and public reference databases are not duplicated in this repository. Public dataset accessions, sample information, and major analysis outputs are provided in the manuscript and Supplementary Tables.

---

## 1. Data inputs

This step defines the main data and reference files used throughout the analysis. Inputs include:

- PacBio Iso-Seq long-read RNA-seq data
- Matched short-read RNA-seq data
- Human GRCh38/hg38 reference genome
- GENCODE v47 gene annotation
- FANTOM5 CAGE and refTSS annotations
- Poly(A) motif annotations used during transcript evaluation
- RepeatMasker repeat annotations/resources
- HLA-I allele list
- Public mass spectrometry data used for peptide-level evaluation

Raw sequencing data and large reference files are not redistributed here. Sample information is summarized in Supplementary Table S1, and transcript-level TE annotations and filtering information are provided in Supplementary Table S2.

---

## 2. Full-length isoform generation and quality control

This step generates and annotates full-length transcript isoforms from PacBio Iso-Seq data.

### Scripts

```text
01_pbmm2.sh
02_collapse.sh
03_sqanti3.sh
04_extract_filter_sqanti3_out.sh
```

These scripts cover long-read alignment, isoform collapsing, SQANTI3 structural annotation, and initial result filtering. The resulting transcript models contain genomic coordinates, exon structures, ORF annotations, and transcription start site support information used in subsequent analyses.

---

## 3. TE annotation and first-exon event identification

This step identifies TE fragments within transcript isoforms and determines their positions relative to transcript exon structures.

### Scripts

```text
05_run_repeatmasker_batch.sh
06_parseRM_step1_get_allinfo.py
07_parseRM_step2_merge_allinfo.py
```

Transcript sequences are annotated with RepeatMasker, and TE coordinates are integrated with transcript and genomic coordinates. The analysis records the exon containing each TE fragment, its boundary relationship with the exon, truncation status, and overlap extent.

TE overlaps are classified according to their positions within transcript structures, including first exons, internal exons, last exons, and intronic regions. First-exon-associated records are further evaluated according to whether the TE overlaps the 5′ portion of exon 1, covers most or nearly all of exon 1, or extends from exon 1 into the following exon.

The resulting annotations form the basis for manual structural curation and event-level consolidation.

---

## 4. Transcription start site support assessment

This step evaluates whether the 5′ ends of TE-associated first-exon isoforms are supported by transcription start site annotations.

### Related script

```text
03_sqanti3.sh
```

SQANTI3 output is used together with CAGE-related annotations and refTSS information to assess transcription start site support. FANTOM5 CAGE data can also be inspected manually for selected events.

FANTOM5: https://fantom.gsc.riken.jp/5/

Transcription start site support is used as one layer of evidence during manual structural review and is not used alone to define the final event catalog.

---

## 5. Matched short-read RNA-seq support

This step evaluates splice-junction support for the long-read isoforms using matched short-read RNA-seq data.

### Scripts

```text
09_get_candidate_junctions.py
10_get_star_junction_support.py
```

Splice junctions are extracted from transcript GTF structures and matched to STAR `SJ.out.tab` files. A junction is considered matched when chromosome, intron start, intron end, and strand are identical between the transcript-derived junction and the STAR output.

The analysis retains uniquely mapped split-read counts, multi-mapping split-read counts, and maximum overhang values as measures of junction support. The exon 1–exon 2 junction is used as the primary summary of 5′ structural support in the manuscript, while complete junction-level results are provided in Supplementary Tables S6 and S7.

Matched short-read transcript abundance was also estimated against a custom transcriptome reference combining GENCODE v47 transcripts with the long-read isoforms. Transcript abundance was summarized as TPM across the matched short-read RNA-seq datasets, with the corresponding results provided in Supplementary Table S8.

---

## 6. Manual structural curation and event consolidation

This step integrates the evidence generated above and produces the curated event-level catalog used for downstream analyses.

Transcript structures are reviewed using TE annotations, transcript and gene annotations, transcription start site support, long-read evidence, matched short-read splice-junction support, and other relevant structural information.

A single transcript may contain more than one TE fragment. TE fragments that define the first-exon event are distinguished from additional TE fragments present elsewhere in the same isoform. Sample-level observations are subsequently consolidated into nonredundant events for event-level analyses.

This step is based mainly on manual review and result integration and therefore does not have a separate standalone script. Curated sample-level records and nonredundant event-level results are provided in Supplementary Tables S4 and S5.

---

## 7. ORF and peptide analysis

This step evaluates predicted ORFs in the curated TE-associated 5′ isoforms and links coding models to transcript structure and TE annotations.

### Script

```text
11_get_ORF_origin.py
```

ORF annotations are derived from SQANTI3 outputs and integrated with transcript exon structures and TE annotations. The analysis determines the exon containing the ORF initiation site, whether the ORF spans the exon 1–exon 2 junction, and whether the predicted coding region overlaps a TE annotation.

ORFs encoding identical amino acid sequences can be collapsed into nonredundant ORFs. For exon 1-initiated coding models spanning the exon 1–exon 2 junction, junction-spanning peptide sequences are generated for downstream reference-sequence comparison and HLA-I analysis.

ORF-level results are summarized in Supplementary Table S9.

---

## 8. HLA-I binding and antigen processing/presentation prediction

This step evaluates selected peptide sequences for predicted HLA-I binding and MHC class I pathway presentation.

### Scripts

```text
08_get_netMHC_input.py
12_run_netmhcpan_pipeline.py
13_parse_netmhcpan_out.py
14_run_netctlpan_pipeline.py
15_parse_netctlpan_out.py
16_merge_epitope_predictions.py
```

NetMHCpan is used for HLA-I binding prediction, and NetCTLpan is used for complementary prediction of MHC class I pathway presentation. The manuscript analysis uses a predefined panel of 48 HLA-I alleles.

The scripts prepare peptide inputs, run the prediction tools, parse the corresponding outputs, and integrate peptide-level prediction results for downstream analysis.

---

## 9. Normal-tissue and reference-protein background assessment

This step evaluates whether the curated transcript, ORF, or peptide sequences overlap known normal-tissue transcript backgrounds or reference human protein sequences.

### Scripts

```text
17_gffcompare_with_gtexV9SupTab5_ids.py
18_dedup_orf_vs_uniprot.py
```

GTEx normal-tissue transcript information is used to annotate the normal-tissue transcriptional background of TE-associated isoforms. This information is used for interpretation and is not treated as a standalone criterion for defining the final event catalog.

Predicted ORF products and derived peptides are compared with reference human protein sequences, including GENCODE v47 protein-coding translations and reviewed human UniProtKB proteins. Exact sequence matching is used to identify peptides already present in reference proteins. UCSC BLAT (GRCh38/hg38) can be used for manual sequence-source review of selected results.

This step reduces ambiguity arising from known reference proteins, normal-tissue transcript backgrounds, or identical peptide sequences with multiple possible sources.

---

## 10. Public mass spectrometry data reanalysis

This step provides complementary peptide-level evidence for selected predicted coding regions.

Public mass spectrometry raw data are reanalyzed with MaxQuant v2.8.0.0 using a custom protein sequence database containing reviewed human UniProtKB proteins together with the predicted ORF products. Search settings are adjusted according to the experimental design of the corresponding source datasets.

Identified peptides are compared with the predicted ORF and peptide sequences. Representative MS2 spectra are shown in the manuscript and Supplementary Figures.

No standalone script is provided for this step.

---

## Software environment

```text
Linux 8.7, x86_64
Python v3.11.0
Iso-Seq v4.3.0
SQANTI3 v5.3.6
RepeatMasker v4.1.6
Salmon v1.10.2
netMHCpan v4.1
netCTLpan v1.1
MaxQuant v2.8.0.0
```

Additional external tools used in individual steps include pbmm2, STAR, gffcompare, and UCSC BLAT. Paths to software, reference files, and input/output directories should be adjusted for the local computing environment before running the scripts.

---

## Public reference resources

| Resource | URL | Use in this study |
| --- | --- | --- |
| GRCh38/hg38 | https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47/GRCh38.primary_assembly.genome.fa.gz | Long-read alignment, genomic coordinate assignment, and sequence review |
| GENCODE v47 | https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47/gencode.v47.primary_assembly.annotation.gtf.gz | Transcript structure annotation and reference sequence comparison |
| FANTOM5 CAGE | https://fantom.gsc.riken.jp/5/ | Transcription start site support |
| GTEx | https://www.gtexportal.org/home/downloads/adult-gtex/long_read_data | Normal-tissue transcript background assessment |
| UniProtKB | https://www.uniprot.org/uniprotkb?query=Human | Reviewed human protein sequence comparison |
| UCSC BLAT | https://genome.ucsc.edu/cgi-bin/hgBlat | Manual sequence-source review |

---

## Notes

- This repository contains the analysis scripts used to document the principal computational workflow of the study.
- Large public raw datasets and reference databases are not redistributed here.
- Some steps include manual review and therefore are not represented by a single standalone script.
- File paths and resource locations in the scripts may need to be modified before reuse in another computing environment.

