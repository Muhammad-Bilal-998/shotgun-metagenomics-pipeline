#!/bin/bash

#==============================================================================#
#       S H O T G U N   M E T A G E N O M I C S   P I P E L I N E            #
#==============================================================================#
#
# This script automates a complete shotgun metagenomic workflow, from raw SRA
# reads through quality control, assembly, read mapping, binning, and quality
# assessment. The pipeline is modular, allowing users to skip completed steps
# and resume from any checkpoint.
#
# Author: Muhammad Bilal
# Last updated: 2026
#
# Usage:
#   bash shotgun_metagenomics_pipeline.sh
#
# Before running:
#   1. Edit the Editable Section below with your paths and sample IDs
#   2. Ensure your conda environments are set up correctly
#   3. Place your SRA accession list file in the working directory
#
#==============================================================================#

#------------------------------------------------------------------------------#
#                 EDITABLE SECTION: Customize for Your System                  #
#------------------------------------------------------------------------------#

# Main working directory where all output will be stored
WORKING_DIR="/path/to/your/working/directory"

# File containing a list of SRA accessions (one per line)
SRA_ID_LIST="sampleID.txt"

# Number of CPU threads to use for parallel processing
THREADS=16

# Set to "true" to skip steps if output already exists, "false" to rerun all steps
SKIP_EXISTING_STEPS=true

# Conda environment names (edit these to match your system)
ENV_SRA="sra_toolkit"          # fastq-dump
ENV_FASTQC="fastqc"            # FastQC
ENV_FASTP="fastp"              # fastp
ENV_MAPPING="mapping"          # bwa, samtools
ENV_MEGAHIT="megahit"          # megahit
ENV_MAXBIN="maxbin"            # MaxBin, jgi_summarize_bam_contig_depths
ENV_CHECKM="checkm"            # CheckM

#==============================================================================#
#                  DO NOT EDIT ANYTHING BELOW THIS LINE                        #
#==============================================================================#

# Color codes for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Create directories and log file
mkdir -p "$WORKING_DIR"
cd "$WORKING_DIR" || exit 1

LOG_FILE="pipeline_log_$(date +%Y%m%d_%H%M%S).log"
{
    echo "=========================================="
    echo "Shotgun Metagenomic Pipeline Started"
    echo "Start time: $(date)"
    echo "=========================================="
    echo "Working directory: $WORKING_DIR"
    echo "Sample ID list: $SRA_ID_LIST"
    echo "CPU threads: $THREADS"
    echo "Skip existing steps: $SKIP_EXISTING_STEPS"
    echo "=========================================="
} | tee "$LOG_FILE"

# Function to check if a command was successful
check_command() {
    if [ $? -ne 0 ]; then
        echo -e "${RED}[ERROR]${NC} $1 failed. Exiting." | tee -a "$LOG_FILE"
        exit 1
    else
        echo -e "${GREEN}[OK]${NC} $1 completed successfully." | tee -a "$LOG_FILE"
    fi
}

# Function to check if a conda environment exists
check_environment() {
    if ! conda env list | grep -q "^$1 "; then
        echo -e "${RED}[ERROR]${NC} Conda environment '$1' not found." | tee -a "$LOG_FILE"
        echo "Available environments:" | tee -a "$LOG_FILE"
        conda env list | tee -a "$LOG_FILE"
        exit 1
    fi
}

# Verify all required conda environments exist
echo "Checking conda environments..." | tee -a "$LOG_FILE"
for ENV in "$ENV_SRA" "$ENV_FASTQC" "$ENV_FASTP" "$ENV_MAPPING" "$ENV_MEGAHIT" "$ENV_MAXBIN" "$ENV_CHECKM"; do
    check_environment "$ENV"
done
echo -e "${GREEN}All environments found.${NC}" | tee -a "$LOG_FILE"

# Verify input files exist
if [ ! -f "$SRA_ID_LIST" ]; then
    echo -e "${RED}[ERROR]${NC} Sample ID list file '$SRA_ID_LIST' not found." | tee -a "$LOG_FILE"
    exit 1
fi

#------------------------------------------------------------------------------#
# Step 1: Download and convert SRA data to FASTQ                              #
#------------------------------------------------------------------------------#

echo "" | tee -a "$LOG_FILE"
echo "Step 1: Downloading and converting SRA reads to FASTQ format..." | tee -a "$LOG_FILE"

if [ "$SKIP_EXISTING_STEPS" = true ] && [ -d "fastq_reads" ] && [ "$(ls -A fastq_reads 2>/dev/null)" ]; then
    echo "  Skipping: fastq_reads directory already populated" | tee -a "$LOG_FILE"
else
    rm -rf fastq_reads
    mkdir -p fastq_reads
    while IFS= read -r SRA_ID; do
        [ -z "$SRA_ID" ] && continue
        echo "  Processing SRA ID: $SRA_ID" | tee -a "$LOG_FILE"
        conda run -n "$ENV_SRA" fastq-dump --split-3 --gzip --skip-technical --outdir fastq_reads "$SRA_ID" 2>&1 | tee -a "$LOG_FILE"
        check_command "fastq-dump for $SRA_ID"
    done < "$SRA_ID_LIST"
fi

#------------------------------------------------------------------------------#
# Step 2: Initial Quality Control with FastQC                                 #
#------------------------------------------------------------------------------#

echo "" | tee -a "$LOG_FILE"
echo "Step 2: Running initial quality control (FastQC)..." | tee -a "$LOG_FILE"

if [ "$SKIP_EXISTING_STEPS" = true ] && [ -d "fastqc_initial_output" ] && [ "$(ls -A fastqc_initial_output 2>/dev/null)" ]; then
    echo "  Skipping: FastQC reports already generated" | tee -a "$LOG_FILE"
else
    rm -rf fastqc_initial_output
    mkdir -p fastqc_initial_output
    conda run -n "$ENV_FASTQC" fastqc fastq_reads/*.fastq.gz -o fastqc_initial_output -q 2>&1 | tee -a "$LOG_FILE"
    check_command "FastQC initial QC"
fi

#------------------------------------------------------------------------------#
# Step 3: Quality Trimming with fastp                                         #
#------------------------------------------------------------------------------#

echo "" | tee -a "$LOG_FILE"
echo "Step 3: Trimming low-quality bases (fastp)..." | tee -a "$LOG_FILE"

if [ "$SKIP_EXISTING_STEPS" = true ] && [ -d "fastp_output" ] && [ "$(ls -A fastp_output 2>/dev/null)" ]; then
    echo "  Skipping: Trimmed reads already exist" | tee -a "$LOG_FILE"
else
    rm -rf fastp_output
    mkdir -p fastp_output
    for R1_FILE in fastq_reads/*_1.fastq.gz; do
        R2_FILE="${R1_FILE/_1.fastq.gz/_2.fastq.gz}"
        SAMPLE_NAME=$(basename "$R1_FILE" _1.fastq.gz)
        OUTPUT_R1="fastp_output/${SAMPLE_NAME}_trimmed_1.fastq.gz"
        OUTPUT_R2="fastp_output/${SAMPLE_NAME}_trimmed_2.fastq.gz"
        JSON_REPORT="fastp_output/${SAMPLE_NAME}_fastp.json"
        HTML_REPORT="fastp_output/${SAMPLE_NAME}_fastp.html"

        if [ -f "$R2_FILE" ]; then
            echo "  Trimming: $SAMPLE_NAME" | tee -a "$LOG_FILE"
            conda run -n "$ENV_FASTP" fastp -i "$R1_FILE" -o "$OUTPUT_R1" -I "$R2_FILE" -O "$OUTPUT_R2" \
                --thread "$THREADS" --json "$JSON_REPORT" --html "$HTML_REPORT" 2>&1 | tee -a "$LOG_FILE"
            check_command "fastp trimming for $SAMPLE_NAME"
        else
            echo "  Warning: No R2 file found for $SAMPLE_NAME, skipping." | tee -a "$LOG_FILE"
        fi
    done
fi

#------------------------------------------------------------------------------#
# Step 4: Post-Trimming QC with FastQC                                        #
#------------------------------------------------------------------------------#

echo "" | tee -a "$LOG_FILE"
echo "Step 4: Running post-trimming quality control (FastQC)..." | tee -a "$LOG_FILE"

if [ "$SKIP_EXISTING_STEPS" = true ] && [ -d "fastqc_trimmed_output" ] && [ "$(ls -A fastqc_trimmed_output 2>/dev/null)" ]; then
    echo "  Skipping: Post-trimming FastQC reports already generated" | tee -a "$LOG_FILE"
else
    rm -rf fastqc_trimmed_output
    mkdir -p fastqc_trimmed_output
    conda run -n "$ENV_FASTQC" fastqc fastp_output/*.fastq.gz -o fastqc_trimmed_output -q 2>&1 | tee -a "$LOG_FILE"
    check_command "FastQC post-trimming QC"
fi

#------------------------------------------------------------------------------#
# Step 5: Assembly with MEGAHIT                                               #
#------------------------------------------------------------------------------#

echo "" | tee -a "$LOG_FILE"
echo "Step 5: Assembling trimmed reads (MEGAHIT)..." | tee -a "$LOG_FILE"

if [ "$SKIP_EXISTING_STEPS" = true ] && [ -d "megahit_assembly" ] && [ -f "megahit_assembly/final.contigs.fa" ]; then
    echo "  Skipping: Assembled contigs already exist" | tee -a "$LOG_FILE"
else
    rm -rf megahit_assembly
    
    # Collect R1 and R2 files
    R1_FILES=$(find fastp_output -name "*_trimmed_1.fastq.gz" | sort | paste -sd ',' -)
    R2_FILES=$(find fastp_output -name "*_trimmed_2.fastq.gz" | sort | paste -sd ',' -)
    
    if [ -z "$R1_FILES" ]; then
        echo -e "${RED}[ERROR]${NC} No trimmed fastq files found in fastp_output/" | tee -a "$LOG_FILE"
        exit 1
    fi
    
    echo "  Running MEGAHIT assembly with $THREADS threads..." | tee -a "$LOG_FILE"
    conda run -n "$ENV_MEGAHIT" megahit -1 "$R1_FILES" -2 "$R2_FILES" -o megahit_assembly -t "$THREADS" 2>&1 | tee -a "$LOG_FILE"
    check_command "MEGAHIT assembly"
    
    echo "  Assembly complete. Contig count:" | tee -a "$LOG_FILE"
    grep -c ">" megahit_assembly/final.contigs.fa | tee -a "$LOG_FILE"
fi

#------------------------------------------------------------------------------#
# Step 6: Read Mapping and Depth Calculation                                  #
#------------------------------------------------------------------------------#

echo "" | tee -a "$LOG_FILE"
echo "Step 6: Mapping reads to contigs and calculating coverage depth..." | tee -a "$LOG_FILE"

mkdir -p depth_calculation

CONTIGS="megahit_assembly/final.contigs.fa"

# Index contigs (once per run)
if [ ! -f "${CONTIGS}.bwt" ]; then
    echo "  Indexing contigs with BWA..." | tee -a "$LOG_FILE"
    conda run -n "$ENV_MAPPING" bwa index "$CONTIGS" 2>&1 | tee -a "$LOG_FILE"
    check_command "BWA indexing"
fi

# Map reads for each sample
if [ "$SKIP_EXISTING_STEPS" = true ] && [ -f "depth_calculation/depth.txt" ]; then
    echo "  Skipping: Depth file already calculated" | tee -a "$LOG_FILE"
else
    while IFS= read -r SAMPLE; do
        [ -z "$SAMPLE" ] && continue
        
        READ1="fastp_output/${SAMPLE}_trimmed_1.fastq.gz"
        READ2="fastp_output/${SAMPLE}_trimmed_2.fastq.gz"
        BAM_SORTED="depth_calculation/${SAMPLE}_aligned_reads.sorted.bam"
        
        if [ -f "$BAM_SORTED" ] && [ "$SKIP_EXISTING_STEPS" = true ]; then
            echo "  Skipping $SAMPLE: BAM already aligned" | tee -a "$LOG_FILE"
            continue
        fi
        
        if [ ! -f "$READ1" ]; then
            echo "  Warning: Trimmed reads not found for $SAMPLE, skipping." | tee -a "$LOG_FILE"
            continue
        fi
        
        echo "  Aligning reads for $SAMPLE..." | tee -a "$LOG_FILE"
        conda run -n "$ENV_MAPPING" bash -c "
            bwa mem -t $THREADS '$CONTIGS' '$READ1' '$READ2' 2>/dev/null | \
            samtools view -bS - 2>/dev/null | \
            samtools sort -o '$BAM_SORTED' 2>/dev/null && \
            samtools index '$BAM_SORTED' 2>/dev/null
        " 2>&1 | tee -a "$LOG_FILE"
        check_command "BWA alignment for $SAMPLE"
    done < "$SRA_ID_LIST"
    
    # Calculate depth across all BAM files
    echo "  Calculating coverage depth..." | tee -a "$LOG_FILE"
    conda run -n "$ENV_MAXBIN" jgi_summarize_bam_contig_depths --outputDepth depth_calculation/depth.txt depth_calculation/*.sorted.bam 2>&1 | tee -a "$LOG_FILE"
    check_command "Depth calculation"
fi

#------------------------------------------------------------------------------#
# Step 7: Binning with MaxBin                                                 #
#------------------------------------------------------------------------------#

echo "" | tee -a "$LOG_FILE"
echo "Step 7: Binning contigs into genome assemblies (MaxBin)..." | tee -a "$LOG_FILE"

if [ "$SKIP_EXISTING_STEPS" = true ] && [ -d "maxbin_output" ] && [ "$(ls -A maxbin_output 2>/dev/null)" ]; then
    echo "  Skipping: MaxBin results already exist" | tee -a "$LOG_FILE"
else
    rm -rf maxbin_output
    mkdir -p maxbin_output
    
    R1_READS=$(find fastp_output -name "*_trimmed_1.fastq.gz" | sort | tr '\n' ' ')
    R2_READS=$(find fastp_output -name "*_trimmed_2.fastq.gz" | sort | tr '\n' ' ')
    
    echo "  Running MaxBin with depth file..." | tee -a "$LOG_FILE"
    conda run -n "$ENV_MAXBIN" run_MaxBin.pl -contig "$CONTIGS" -out maxbin_output/maxbin_out \
        -reads $R1_READS -reads2 $R2_READS -abund depth_calculation/depth.txt -thread "$THREADS" 2>&1 | tee -a "$LOG_FILE"
    check_command "MaxBin binning"
    
    BIN_COUNT=$(find maxbin_output -name "*.fasta" -o -name "*.fa" | wc -l)
    echo "  Binning complete. Bins generated: $BIN_COUNT" | tee -a "$LOG_FILE"
fi

#------------------------------------------------------------------------------#
# Step 8: Quality Assessment with CheckM                                      #
#------------------------------------------------------------------------------#

echo "" | tee -a "$LOG_FILE"
echo "Step 8: Assessing bin quality (CheckM)..." | tee -a "$LOG_FILE"

if [ "$SKIP_EXISTING_STEPS" = true ] && [ -d "checkm_output" ] && [ "$(ls -A checkm_output 2>/dev/null)" ]; then
    echo "  Skipping: CheckM assessment already complete" | tee -a "$LOG_FILE"
else
    rm -rf checkm_output
    mkdir -p checkm_output
    
    # Verify bins exist
    if [ -z "$(find maxbin_output -name '*.fasta' -o -name '*.fa' 2>/dev/null)" ]; then
        echo -e "${YELLOW}[WARNING]${NC} No bins found in maxbin_output/. Skipping CheckM analysis." | tee -a "$LOG_FILE"
    else
        echo "  Running CheckM lineage workflow..." | tee -a "$LOG_FILE"
        conda run -n "$ENV_CHECKM" checkm lineage_wf -x fasta -t "$THREADS" maxbin_output checkm_output 2>&1 | tee -a "$LOG_FILE"
        check_command "CheckM lineage workflow"
        
        echo "  Generating CheckM QA report..." | tee -a "$LOG_FILE"
        conda run -n "$ENV_CHECKM" checkm qa -o 2 -f checkm_output/checkm_qa_results.txt \
            checkm_output/lineage.ms checkm_output 2>&1 | tee -a "$LOG_FILE"
        check_command "CheckM QA report"
        
        echo "  Generating quality assessment plots..." | tee -a "$LOG_FILE"
        conda run -n "$ENV_CHECKM" checkm gc_plot checkm_output/ maxbin_output/ checkm_output/gc_plot.pdf 2>&1 | tee -a "$LOG_FILE"
        conda run -n "$ENV_CHECKM" checkm coding_plot checkm_output/ maxbin_output/ checkm_output/coding_plot.pdf 2>&1 | tee -a "$LOG_FILE"
        
        echo -e "${GREEN}Quality plots generated.${NC}" | tee -a "$LOG_FILE"
    fi
fi

#------------------------------------------------------------------------------#
# Pipeline Complete                                                            #
#------------------------------------------------------------------------------#

echo "" | tee -a "$LOG_FILE"
echo "=========================================="
echo -e "${GREEN}PIPELINE COMPLETED SUCCESSFULLY${NC}"
echo "=========================================="
echo "Completion time: $(date)" | tee -a "$LOG_FILE"
echo "" | tee -a "$LOG_FILE"
echo "Output summary:" | tee -a "$LOG_FILE"
echo "  - Trimmed reads: fastp_output/" | tee -a "$LOG_FILE"
echo "  - Assembled contigs: megahit_assembly/final.contigs.fa" | tee -a "$LOG_FILE"
echo "  - Binned genomes: maxbin_output/" | tee -a "$LOG_FILE"
echo "  - Quality assessment: checkm_output/checkm_qa_results.txt" | tee -a "$LOG_FILE"
echo "  - Pipeline log: $LOG_FILE" | tee -a "$LOG_FILE"
echo "=========================================="
