#!/bin/bash

# --- Configuration ---
DATA_DIR="/home/yxqin/neoantigen/Analysis/pacbio/00_data/PRJNA1176011/flncBAM"
REFERENCE_GENOME="/home/yxqin/neoantigen/reference/GRCh38.primary_assembly.genome.fa"
OUTPUT_BAM_DIR="/home/yxqin/neoantigen/Analysis/pacbio/02_pbmm2/PRJNA1176011"
LOG_DIR="/home/yxqin/neoantigen/Analysis/pacbio/00_log/pbmm2"

set -euo pipefail

echo "$(date): Starting pbmm2 alignment tasks."

mkdir -p "$OUTPUT_BAM_DIR"
mkdir -p "$LOG_DIR"

find "$DATA_DIR" -maxdepth 1 -type f -name "*.flnc.bam" | while read flnc_bam_path; do
    FLNC_BAM_BASENAME=$(basename "$flnc_bam_path")
    SAMPLE_ID="${FLNC_BAM_BASENAME%.flnc.bam}"

    OUTPUT_ALIGNED_BAM="${OUTPUT_BAM_DIR}/${SAMPLE_ID}.sorted.aligned.bam"
    LOG_FILE="${LOG_DIR}/${SAMPLE_ID}_pbmm2.log"

    echo "$(date): Processing sample: $SAMPLE_ID"

    if [ -f "$OUTPUT_ALIGNED_BAM" ]; then
        echo "$(date): Output file '$OUTPUT_ALIGNED_BAM' already exists. Skipping alignment for $SAMPLE_ID."
        continue
    fi

    echo "$(date): Aligning $FLNC_BAM_BASENAME..."
    
    pbmm2 align \
        --preset ISOSEQ \
        --sort \
        "$flnc_bam_path" \
        "$REFERENCE_GENOME" \
        "$OUTPUT_ALIGNED_BAM" \
        > "$LOG_FILE" 2>&1
    
    PBMM2_STATUS=$?

    if [ $PBMM2_STATUS -eq 0 ]; then
        echo "$(date): Alignment for $SAMPLE_ID completed successfully. Output: $OUTPUT_ALIGNED_BAM"
    else
        echo "$(date): Error: Alignment for $SAMPLE_ID failed with exit code $PBMM2_STATUS. Check $LOG_FILE for details."
    fi

done

echo "$(date): All pbmm2 alignment tasks finished."