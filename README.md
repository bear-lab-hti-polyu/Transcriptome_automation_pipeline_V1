
#  Transcriptome Assembly Pipeline (v1.1)

**Date:** 2026-06-06  
**Support:** Multi-sample & Multi-platform (Illumina, Nanopore, PacBio)

---

##  Overview
This pipeline provides a robust, automated workflow for **transcriptome assembly** and **variant calling**. It is designed to handle single or multiple sequencing samples, guiding data from raw FASTQ files through quality control, alignment, and consensus generation to final protein prediction.

### Key Features
*   **Multi-Platform:** Optimized presets for Illumina (short-read), Nanopore, and PacBio (long-read).
*   **Sample Merging:** Automatically merges multiple input files into a single unified assembly.
*   **End-to-End:** Integrated variant calling (SNPs) and ORF prediction (TransDecoder).

---

##  Prerequisites & Dependencies
It is highly recommended to manage the following dependencies via **Conda** or **Mamba**. Ensure these tools are accessible in your `$PATH`:

| Category | Tools |
| :--- | :--- |
| **Quality Control** | `fastp`, `NanoFilt` |
| **Alignment** | `hisat2`, `minimap2`, `samtools` |
| **Assembly** | `stringtie` |
| **Variant Calling** | `bcftools` |
| **Sequence Processing** | `gffread`, `bedtools` |
| **ORF Prediction** | `TransDecoder` |

---

##  Usage

The script is flexible. You can process individual files or a comma-separated list for multi-sample merging.

### Basic Syntax
```bash
./transcriptome_pipeline_multi.sh -i <input_files> -r <reference_fasta> -o <output_prefix> [options]
```

### Examples
**1. Illumina (Short-read) - Single Sample:**
```bash
./transcriptome_pipeline_multi.sh -i sample1.fastq -r genome.fna -o output_name
```

**2. Nanopore (Long-read) - Multiple Samples (Skipping QC):**
```bash
./transcriptome_pipeline_multi.sh -i s1.fq,s2.fq,s3.fq -r genome.fna -o merged_output -m nanopore -s
```

---

##  Parameters

| Flag | Description | Default |
| :--- | :--- | :--- |
| `-i` | **Required:** Input FASTQ file(s) (comma-separated for multiple) | - |
| `-r` | **Required:** Reference genome FASTA file | - |
| `-o` | **Required:** Output prefix / Sample name | - |
| `-m` | **Mode:** `illumina`, `nanopore`, or `pacbio` | `illumina` |
| `-t` | Number of CPU threads | `32` |
| `-p` | Organism ploidy | `1` |
| `-q` | Minimum base quality score | `20` |
| `-d` | Minimum variant depth for filtering | `10` |
| `-s` | **Skip** quality filtering (use raw data directly) | Off |
| `-h` | Display help information | - |

---

##  Pipeline Workflow

1.  **Quality Filtering:** Raw data cleaning using `fastp` (Illumina) or `NanoFilt` (Long-read).
2.  **Mapping:** Splice-aware alignment to the reference via `hisat2` or `minimap2`.
3.  **Merging:** Automated merging of BAM files if multiple inputs are provided.
4.  **Assembly:** Transcript reconstruction using `StringTie`.
5.  **Variant Calling:** SNP/Indel identification using `bcftools`.
6.  **Consensus:** Generation of a sample-specific consensus genome.
7.  **Extraction:** Transcript sequence retrieval via `gffread` or `bedtools`.
8.  **Protein Prediction:** Identification of coding regions (CDS) and proteins using `TransDecoder`.

---

##  Directory Structure

Upon completion, the following directories will be created:

*   `02_clean/` : Filtered FASTQ files and QC reports.
*   `04_mapping/` : Alignment BAM files and mapping statistics (`flagstat`).
*   `05_stringtie/` : Final transcript annotations (GTF).
*   `06_variant/` : Filtered VCF files and variant stats.
*   `07_consensus/` : Sample-specific consensus genome (FASTA).
*   `08_output_transcript/` : Final extracted transcript sequences.
*   `09_proteins/` : Predicted protein sequences and CDS.

---

##  Notes
*   **Storage:** Ensure you have sufficient disk space. Intermediate files (BAMs) can be large; the pipeline typically requires **5–10x** the size of your raw input data.
*   **Memory:** For large genomes or high-depth long-read data, ensure your system has adequate RAM for `minimap2` and `stringtie` operations.
