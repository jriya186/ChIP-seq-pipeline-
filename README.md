# RUNX1 ChIP-seq Analysis Pipeline

A complete ChIP-seq data analysis pipeline for studying **RUNX1 transcription factor binding** in breast cancer (MCF-7) cells, built with Nextflow (DSL2) and analyzed using a Python Jupyter notebook.

---

## Note on Results and Images

This repository contains the pipeline code and analysis notebook, but **result files, output images, and intermediate data are not included**. These were generated on a university HPC server (Boston University SCC / `projectnb/bf528`) to which access is no longer available, and the raw data files were too large to upload to GitHub. As a result:

- Images embedded in the Jupyter notebook (`analysis.ipynb`) will not render — their source files are missing.
- Output files referenced in the notebook (e.g., `results/annotated_peaks.txt`, BigWig files, motif results) are not present.
- The notebook can be read as a documented record of the analysis steps, findings, and interpretations, even without the rendered figures.

---

## Project Overview

This project processes and analyzes RUNX1 ChIP-seq data to characterize the genome-wide binding landscape of RUNX1 in MCF-7 breast cancer cells. The analysis replicates and compares results from a published study (GEO accession: GSE75070), covering everything from raw FASTQ processing through to pathway enrichment.

**Reference genome:** hg38  
**Cell line:** MCF-7 (human breast cancer)  
**Target protein:** RUNX1 (transcription factor)  
**Replicates:** 2 IP replicates + 2 Input controls

---

## Repository Structure

```
ChIP-seq-pipeline/
├── main.nf                  # Main Nextflow pipeline (DSL2)
├── modules/                 # Individual tool modules
│   ├── fastqc/
│   ├── trimmomatic/
│   ├── bowtie2_build/
│   ├── bowtie2_align/
│   ├── samtools_sort/
│   ├── samtools_idx/
│   ├── samtools_flagstat/
│   ├── deeptools_bamcoverage/
│   ├── multiqc/
│   ├── deeptools_multibwsummary/
│   ├── deeptools_plotcorrelation/
│   ├── macs3_callpeak/
│   ├── bedtools_intersect/
│   ├── bedtools_remove/
│   ├── homer_annotatepeaks/
│   ├── deeptools_computematrix/
│   ├── deeptools_plotprofile/
│   └── findmotif/
├── analysis.ipynb           # Downstream analysis and visualizations
└── README.md
```

---

## Pipeline Overview (`main.nf`)

The pipeline is written in **Nextflow DSL2** and organized into modular processes. It takes a CSV samplesheet as input and processes all samples through the following stages:

### 1. Quality Control & Preprocessing
- **FastQC** — raw read quality assessment
- **Trimmomatic** — adapter trimming and low-quality base removal
- **MultiQC** — aggregates FastQC, trimming logs, and alignment stats into one report

### 2. Alignment
- **Bowtie2** — genome index build (`BOWTIE2_BUILD`) and read alignment (`BOWTIE2_ALIGN`) to hg38
- **SAMtools** — BAM sorting, indexing, and alignment statistics (`flagstat`)

### 3. Signal Track Generation
- **deepTools `bamCoverage`** — generates BigWig coverage tracks from sorted BAMs

### 4. Sample Correlation
- **deepTools `multiBigwigSummary`** — summarizes signal across BigWig files
- **deepTools `plotCorrelation`** — Spearman correlation heatmap between samples

### 5. Peak Calling & Filtering
- **MACS3** — peak calling (note: due to package issues, pre-called narrowPeak files from course materials were used for downstream steps)
- **bedtools intersect** (`-f 0.5 -r`) — reproducible peak identification with 50% reciprocal overlap between replicates
- **bedtools intersect** (`-v`) — blacklist region removal using the hg38 ENCODE blacklist

### 6. Peak Annotation & Downstream Analysis
- **HOMER `annotatePeaks.pl`** — annotates peaks with genomic feature context (promoter-TSS, intron, intergenic, etc.) using hg38 FASTA and GTF
- **deepTools `computeMatrix`** + **`plotProfile`** — signal intensity profiles across gene bodies (±2000 bp padding) for IP replicates
- **HOMER `findMotifsGenome.pl`** — de novo and known motif enrichment on filtered peaks

---

## Downstream Analysis (`analysis.ipynb`)

The Jupyter notebook documents the post-pipeline analysis and interpretation. Sections include:

### Quality Control
Assessment of the MultiQC report across all four samples (IP_rep1, IP_rep2, Input_rep1, Input_rep2). Key findings:
- High base quality (Phred > 30) across all samples
- Effective adapter removal by Trimmomatic
- **Elevated duplication rates in IP samples** — IP_rep2 exceeded 74% duplicated reads, suggesting PCR overamplification or low library complexity. Downstream duplicate removal was recommended.

### Reproducible Peaks & Blacklist Filtering
- Spearman correlation was chosen for BigWig comparison as it is rank-based, robust to outliers, and appropriate for genomic data that often violates normality assumptions.
- Reproducible peaks were defined by 50% reciprocal overlap between replicates using `bedtools intersect -f 0.5 -r`.
- Blacklisted regions were removed using `bedtools intersect -v`, retaining only clean, high-confidence peaks.

### Signal Intensity Profiles *(images not available)*
deepTools `computeMatrix` + `plotProfile` were used to visualize read enrichment across gene bodies in IP_rep1 and IP_rep2. Both replicates showed strong signal enrichment at the **transcription start site (TSS)**, consistent with RUNX1's known role as a promoter-binding transcription factor.

### Motif Enrichment *(table not available)*
HOMER motif analysis on filtered peaks revealed:
- Strong enrichment for **RUNX family motifs** (RUNX1, RUNX2, RUNX-AML) — present in >30.5% of target sequences vs. background, confirming on-target ChIP enrichment.
- Enrichment of **Forkhead family factors** (FOXA1, FOXA2, FOXM1), suggesting cooperative regulatory activity.
- Additional enrichment of GRHL2 and FOSL2 motifs.

### Comparison with Published Study (GSE75070)
The notebook replicates key figures from the original publication:

- **Figure 2D/2E (Venn diagram):** Overlap between replicate peak sets was visualized. Discrepancies from the published counts were attributed to differences in peak-calling parameters, q-value thresholds, or overlap methods (IDR vs. coordinate overlap).
- **Figure 2F (RUNX1-bound DEGs):** RNA-seq data (MCF-7 shRUNX1 vs. shNS) was downloaded from GEO and intersected with ChIP-seq peaks to quantify RUNX1 binding at differentially expressed genes (±5 kb TSS and ±20 kb gene body). Fewer overlaps were found compared to the paper, likely due to differences in gene boundary definitions, annotation tools, or peak sets used.
- **Read count table (Supplementary Table S2A):** Raw and mapped read counts were compared across samples. INPUT_rep2 showed notably lower sequencing depth (~10.9M reads) and missing mapped read data in MultiQC output.

### Pathway Enrichment *(images not available)*
Promoter-TSS annotated peaks were extracted and gene names submitted to **Enrichr** for pathway analysis against KEGG and Reactome databases. Top enriched pathways included:
- **Reactome:** "Formation of WDR5-containing Histone-Modifying Complexes", "Gene Expression", "RNA Transport"
- **KEGG:** Cellular senescence, autophagy

These results are consistent with RUNX1's known role in chromatin remodeling and transcriptional regulation in breast cancer.

---

## Tools & Dependencies

| Tool | Version / Notes | Purpose |
|---|---|---|
| Nextflow | DSL2 | Pipeline orchestration |
| FastQC | — | Raw read QC |
| Trimmomatic | — | Adapter trimming |
| Bowtie2 | — | Read alignment to hg38 |
| SAMtools | — | BAM processing |
| deepTools | — | BigWig, matrix, correlation, profile |
| MACS3 | — | Peak calling |
| bedtools | — | Peak intersection & filtering |
| HOMER | — | Peak annotation & motif finding |
| MultiQC | — | QC report aggregation |
| Python | 3.13 | Downstream analysis |
| pandas | — | Data manipulation |
| matplotlib | — | Plotting |
| matplotlib-venn | — | Venn diagrams |
| requests | — | GEO data download |
| Singularity | — | Containerized tool execution |
| Conda | — | Environment management |

---

## How to Run

### Prerequisites
- Nextflow installed
- Singularity (or Docker) for containerized modules
- Access to hg38 reference genome FASTA and GTF

### Configuration
Create a `params` file or pass parameters at runtime:
```bash
nextflow run main.nf \
  --samplesheet samplesheet.csv \
  --genome /path/to/hg38.fa \
  --gtf /path/to/hg38.gtf \
  --adapter_fa /path/to/adapters.fa \
  --blacklist /path/to/hg38_blacklist.bed
```

### Samplesheet Format
The pipeline expects a CSV with `name` and `path` columns:
```csv
name,path
IP_rep1,/path/to/IP_rep1.fastq.gz
IP_rep2,/path/to/IP_rep2.fastq.gz
Input_rep1,/path/to/Input_rep1.fastq.gz
Input_rep2,/path/to/Input_rep2.fastq.gz
```

---

## Key Findings Summary

- RUNX1 binding in MCF-7 cells is strongly enriched at **gene promoters and TSSs**, consistent with its function as a transcription factor.
- Motif analysis confirms on-target ChIP enrichment, with RUNX family motifs as the top hit.
- Pathway enrichment links RUNX1-bound genes to **transcriptional regulation, chromatin remodeling, and RNA processing**, as well as cellular stress responses relevant to breast cancer biology.
- Results broadly recapitulate the published GSE75070 dataset, with quantitative differences attributable to methodological variations.

---

## Data Availability

Raw ChIP-seq data is publicly available on NCBI GEO: **[GSE75070](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE75070)**
