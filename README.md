# shotgun-metagenomics-pipeline
Automated bash pipeline for end-to-end quality control, assembly, and profiling of shotgun metagenomic data.
# Automated Shotgun Metagenomic Analysis Pipeline

[![Language: Bash](https://img.shields.io/badge/Language-Bash-4EAA25.svg)](https://www.gnu.org/software/bash/)
[![Environment: Conda](https://img.shields.io/badge/Environment-Conda-green.svg)](https://docs.conda.io/)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)

An automated, checkpoint-aware Bash pipeline for end-to-end processing of shotgun metagenomic sequencing data. It orchestrates raw SRA retrieval, quality trimming, *de novo* assembly, contig-level coverage mapping, metagenomic binning, and quality assessment using isolated Conda environments.

---

## Key Features

- **Automated End-to-End Orchestration:** Chains SRA ingestion through to CheckM bin assessment without manual intervention.
- **Environment Isolation:** Uses dynamic `conda run` calls to manage tool-specific environments independently, preventing dependency and Python version conflicts.
- **Checkpoint-Aware & Resumable:** Built-in verification (`SKIP_EXISTING_STEPS`) skips already-computed stages, allowing straightforward recovery from interruptions.
- **Differential Coverage Binning:** Uses `bwa mem` and `jgi_summarize_bam_contig_depths` to calculate precise multi-sample contig depth profiles for MaxBin2.
- **Comprehensive Logging:** Real-time dual output (stdout and timestamped log file) with status checks (`[OK]` / `[ERROR]`) after each computational stage.

---

## Workflow Overview

| Step | Stage | Tool | Description |
| :---: | :--- | :--- | :--- |
| **1** | SRA Retrieval | **SRA Toolkit** (`fastq-dump`) | Fetches accessions and dumps paired-end FASTQ files (`--split-3 --gzip`). |
| **2** | Raw Read QC | **FastQC** | Evaluates base quality, adapter contamination, and GC content. |
| **3** | Preprocessing | **fastp** | Performs adapter clipping, poly-G/poly-X trimming, and low-quality filtering. |
| **4** | Post-Trim QC | **FastQC** | Validates read improvements post-filtering. |
| **5** | *De Novo* Assembly | **MEGAHIT** | Co-assembles paired-end reads into high-contiguity contigs. |
| **6** | Read Mapping | **BWA & Samtools** | Maps reads back to contigs to generate indexed, sorted BAM files. |
| **7** | Depth Profiling | **MetaBAT2** (`jgi_summarize_bam_contig_depths`) | Computes contig depth across all aligned BAM files. |
| **8** | MAG Binning | **MaxBin 2.0** | Groups assembled contigs into discrete genome bins using tetranucleotide frequencies and coverage profiles. |
| **9** | Bin Evaluation | **CheckM** (`lineage_wf`) | Evaluates completeness, contamination, coding density, and GC distributions. |

---

## Directory Structure

text
├── shotgun_metagenomics_pipeline.sh   # Main orchestration script
├── sampleID.txt                       # SRA accession list (one per line)
├── README.md                          # Documentation
├── LICENSE                            # MIT License
└── results/                           # Working directory (generated at runtime)
    ├── fastq_reads/                   # Raw FASTQ files
    ├── fastqc_initial_output/         # Raw read QC reports
    ├── fastp_output/                  # Trimmed reads and HTML/JSON summaries
    ├── fastqc_trimmed_output/         # Post-trimming QC reports
    ├── megahit_assembly/              # Final contigs (final.contigs.fa)
    ├── depth_calculation/             # Sorted BAM files and depth.txt
    ├── maxbin_output/                 # Assembled FASTA bins
    └── checkm_output/                 # QA summaries and PDF diagnostic plots

1. Prerequisites & Conda Environments
The script uses dedicated Conda environments for each module to prevent dependency collisions. You can create them using conda or mamba:

Bash
# SRA Toolkit
conda create -y -n sra_toolkit -c bioconda sra-tools

# Quality Control
conda create -y -n fastqc -c bioconda fastqc
conda create -y -n fastp -c bioconda fastp

# Assembly & Mapping
conda create -y -n megahit -c bioconda megahit
conda create -y -n mapping -c bioconda bwa samtools

# Binning & Depth Profiling
conda create -y -n maxbin -c bioconda maxbin2 metabat2

# Quality Assessment
conda create -y -n checkm -c bioconda -c conda-forge checkm-genome
2. Input Configuration
Create a file named sampleID.txt in your execution directory containing your target SRA accession numbers (one per line):

Plaintext
SRR12345678
SRR12345679
3. Pipeline Configuration
Open shotgun_metagenomics_pipeline.sh and configure the editable header section:

Bash
# Main working directory where all output will be stored
WORKING_DIR="/path/to/your/working/directory"

# File containing list of SRA accessions
SRA_ID_LIST="sampleID.txt"

# CPU threads allocated
THREADS=16

# Skip already completed steps (true/false)
SKIP_EXISTING_STEPS=true
4. Execution
Grant execution permissions and run the script:

Bash
chmod +x shotgun_metagenomics_pipeline.sh
./shotgun_metagenomics_pipeline.sh
Outputs
Upon successful completion, the pipeline produces:

megahit_assembly/final.contigs.fa: Assembled metagenome contigs.

maxbin_output/*.fasta: Recovered Metagenome-Assembled Genomes (MAGs).

checkm_output/checkm_qa_results.txt: Tabulated completeness and contamination metrics per MAG.

checkm_output/*.pdf: Quality distribution plots (GC content and coding density).

pipeline_log_*.log: Timestamped execution trace.

Author
Muhammad Bilal - https://github.com/Muhammad-Bilal-998
