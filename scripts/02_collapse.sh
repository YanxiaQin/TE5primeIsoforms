#!/bin/bash

# --- Configuration ---
FLNC_BAM_DIR="/home/yxqin/neoantigen/Analysis/pacbio/00_data/PRJNA1176011/flncBAM"
ALIGNED_BAM_DIR="/home/yxqin/neoantigen/Analysis/pacbio/02_pbmm2/PRJNA1176011"
COLLAPSE_GFF_DIR="/home/yxqin/neoantigen/Analysis/pacbio/03_collapse/PRJNA1176011"
LOG_DIR="/home/yxqin/neoantigen/Analysis/pacbio/00_log/collapse"

set -euo pipefail

echo "$(date): Starting Iso-Seq collapse tasks."

# Create output and log directories if they don't exist
mkdir -p "$COLLAPSE_GFF_DIR" || { echo "Error: Cannot create directory $COLLAPSE_GFF_DIR."; exit 1; }
mkdir -p "$LOG_DIR" || { echo "Error: Cannot create directory $LOG_DIR."; exit 1; }

# Iterate over all .flnc.bam files in the input directory
find "$FLNC_BAM_DIR" -maxdepth 1 -type f -name "*.flnc.bam" | while read flnc_bam_path; do
    # Extract sample ID (e.g., SRR31072033)
    FLNC_BAM_BASENAME=$(basename "$flnc_bam_path")
    SAMPLE_ID="${FLNC_BAM_BASENAME%.flnc.bam}"

    # Define paths for current sample
    ALIGNED_BAM="${ALIGNED_BAM_DIR}/${SAMPLE_ID}.sorted.aligned.bam"
    OUTPUT_GFF="${COLLAPSE_GFF_DIR}/${SAMPLE_ID}.collapsed.gff"
    LOG_FILE="${LOG_DIR}/${SAMPLE_ID}_collapse.log"

    echo "$(date): Processing sample: $SAMPLE_ID"

    # Check if final output GFF already exists
    if [ -f "$OUTPUT_GFF" ]; then
        echo "$(date): Output GFF '$OUTPUT_GFF' already exists. Skipping $SAMPLE_ID."
        continue
    fi

    # Check if aligned BAM dependency exists
    if [ ! -f "$ALIGNED_BAM" ]; then
        echo "$(date): Error: Aligned BAM '$ALIGNED_BAM' not found for $SAMPLE_ID. Skipping collapse."
        continue
    fi


    # --- Step 2: Run isoseq collapse ---
    echo "$(date):   Running isoseq collapse for $SAMPLE_ID..."
    isoseq collapse \
        --do-not-collapse-extra-5exons \
        "$ALIGNED_BAM" \
        "$flnc_bam_path" \
        "$OUTPUT_GFF" \
        >> "$LOG_FILE" 2>&1 
    
    COLLAPSE_STATUS=$?

    if [ $COLLAPSE_STATUS -eq 0 ]; then
        echo "$(date):   Iso-Seq collapse completed successfully for $SAMPLE_ID. Output: $OUTPUT_GFF"
    else
        echo "$(date):   Error: Iso-Seq collapse failed for $SAMPLE_ID (exit code $COLLAPSE_STATUS). Check $LOG_FILE for details."
    fi

done

echo "$(date): All Iso-Seq collapse tasks finished."
