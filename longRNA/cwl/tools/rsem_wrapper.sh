#!/bin/bash
set -e

# This script is called by the CWL tool
# Arguments are passed as environment variables

echo "=== RSEM Wrapper Script ===" >&2

if [ "$USE_STRICT_FILTERING" = "true" ]; then
  echo "=== Two-step mode: STAR + RSEM with strict filtering ===" >&2
  echo "Step 1: STAR alignment..." >&2

  # Build STAR command
  STAR_CMD="STAR --runThreadN $THREADS --genomeDir $REF_DIR --readFilesIn $R1"
  [ -n "$R2" ] && STAR_CMD="$STAR_CMD $R2"
  [ "$STAR_GZIPPED" = "true" ] && STAR_CMD="$STAR_CMD --readFilesCommand zcat"

  STAR_CMD="$STAR_CMD --outFileNamePrefix ${SAMPLE_NAME}_STARpass/"
  STAR_CMD="$STAR_CMD --outSAMtype BAM Unsorted"
  STAR_CMD="$STAR_CMD --outSAMunmapped Within"
  STAR_CMD="$STAR_CMD --quantMode TranscriptomeSAM"
  STAR_CMD="$STAR_CMD --outSAMattributes NH HI AS NM MD"

  [ "$PERFECT_MATCH_ONLY" = "true" ] && STAR_CMD="$STAR_CMD --outFilterMismatchNmax 0"
  [ "$UNIQUE_MAPPERS_ONLY" = "true" ] && STAR_CMD="$STAR_CMD --outFilterMultimapNmax 1"

  STAR_CMD="$STAR_CMD --outFilterType BySJout"
  STAR_CMD="$STAR_CMD --alignSJoverhangMin 8"
  STAR_CMD="$STAR_CMD --alignSJDBoverhangMin 1"
  STAR_CMD="$STAR_CMD --outFilterIntronMotifs None"
  STAR_CMD="$STAR_CMD --alignEndsType EndToEnd"

  eval $STAR_CMD

  echo "Step 2: RSEM quantification..." >&2

  # Build RSEM command
  RSEM_CMD="rsem-calculate-expression --alignments"
  [ -n "$R2" ] && RSEM_CMD="$RSEM_CMD --paired-end"
  RSEM_CMD="$RSEM_CMD --num-threads $THREADS"
  [ -n "$FORWARD_PROB" ] && RSEM_CMD="$RSEM_CMD --forward-prob $FORWARD_PROB"
  [ "$ESTIMATE_RSPD" = "true" ] && RSEM_CMD="$RSEM_CMD --estimate-rspd"
  [ -n "$SEED" ] && RSEM_CMD="$RSEM_CMD --seed $SEED"
  [ "$CALC_CI" = "true" ] && RSEM_CMD="$RSEM_CMD --calc-ci"
  RSEM_CMD="$RSEM_CMD --no-bam-output"
  RSEM_CMD="$RSEM_CMD ${SAMPLE_NAME}_STARpass/Aligned.toTranscriptome.out.bam"
  RSEM_CMD="$RSEM_CMD $REF_DIR/$REF_NAME"
  RSEM_CMD="$RSEM_CMD $SAMPLE_NAME"

  eval $RSEM_CMD > ${SAMPLE_NAME}.rsem.log 2> ${SAMPLE_NAME}.rsem.err

  # Copy STAR logs
  [ -d ${SAMPLE_NAME}.temp ] && rm -rf ${SAMPLE_NAME}.temp
  mkdir -p ${SAMPLE_NAME}.temp
  cp ${SAMPLE_NAME}_STARpass/Log.* ${SAMPLE_NAME}.temp/ 2>/dev/null || true

else
  echo "=== Standard mode: RSEM integrated STAR ===" >&2

  # Build RSEM command with integrated STAR
  RSEM_CMD="rsem-calculate-expression"
  [ "$USE_STAR" = "true" ] && RSEM_CMD="$RSEM_CMD --star"
  [ "$STAR_GZIPPED" = "true" ] && RSEM_CMD="$RSEM_CMD --star-gzipped-read-file"
  [ -n "$R2" ] && RSEM_CMD="$RSEM_CMD --paired-end"
  [ -n "$FORWARD_PROB" ] && RSEM_CMD="$RSEM_CMD --forward-prob $FORWARD_PROB"
  [ "$ESTIMATE_RSPD" = "true" ] && RSEM_CMD="$RSEM_CMD --estimate-rspd"
  [ "$CALC_CI" = "true" ] && RSEM_CMD="$RSEM_CMD --calc-ci"
  [ -n "$SEED" ] && RSEM_CMD="$RSEM_CMD --seed $SEED"
  RSEM_CMD="$RSEM_CMD --num-threads $THREADS"
  RSEM_CMD="$RSEM_CMD --star-output-genome-bam"
  RSEM_CMD="$RSEM_CMD $R1"
  [ -n "$R2" ] && RSEM_CMD="$RSEM_CMD $R2"
  RSEM_CMD="$RSEM_CMD $REF_DIR/$REF_NAME"
  RSEM_CMD="$RSEM_CMD $SAMPLE_NAME"

  eval $RSEM_CMD > ${SAMPLE_NAME}.rsem.log 2> ${SAMPLE_NAME}.rsem.err

  # Create temp directory with STAR log placeholders for CWL compatibility
  # In standard mode, RSEM runs STAR internally and doesn't output separate STAR logs
  [ -d ${SAMPLE_NAME}.temp ] && rm -rf ${SAMPLE_NAME}.temp
  mkdir -p ${SAMPLE_NAME}.temp
  echo "STAR logs not available in RSEM standard mode (integrated STAR)" > ${SAMPLE_NAME}.temp/Log.final.out
  echo "STAR logs not available in RSEM standard mode (integrated STAR)" > ${SAMPLE_NAME}.temp/Log.progress.out
  echo "STAR logs not available in RSEM standard mode (integrated STAR)" > ${SAMPLE_NAME}.temp/Log.out
fi

echo "Complete!" >&2
