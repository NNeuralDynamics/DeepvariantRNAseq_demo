#!/bin/bash
set -euo pipefail

# =============================================================================
# Usage Information
# =============================================================================
usage() {
  echo "Usage: $0 --sample SAMPLE_NAME --input_bam INPUT_BAM --output_dir OUTPUT_DIR --ref_fasta REF_FASTA [--threads THREADS]"
  echo ""
  echo "  --sample        Sample name"
  echo "  --input_bam     Input BAM file"
  echo "  --output_dir    Directory where output will be stored"
  echo "  --ref_fasta     Reference genome FASTA file"
  echo "  --threads       (Optional) Number of threads for parallel processing [default: 32]"
  exit 1
}

# Check for minimum required arguments
if [ "$#" -lt 8 ]; then
  usage
fi

# =============================================================================
# Command-line Argument Parsing
# =============================================================================
THREADS=32  # Default value
while [[ "$#" -gt 0 ]]; do
  case "$1" in
    --sample) SAMPLE_NAME="$2"; shift 2;;
    --input_bam) INPUT_BAM="$2"; shift 2;;
    --output_dir) OUTPUT_DIR="$2"; shift 2;;
    --ref_fasta) REF_FASTA="$2"; shift 2;;
    --threads) THREADS="$2"; shift 2;;
    *) echo "Error: Unknown parameter $1"; usage;;
  esac
done

# Define file paths
DICT="${OUTPUT_DIR}/${SAMPLE_NAME}_genome.dict"
RG_BAM="${OUTPUT_DIR}/${SAMPLE_NAME}_rg.bam"
DEDUP_BAM="${OUTPUT_DIR}/${SAMPLE_NAME}_dedup.bam"
SPLIT_BAM="${OUTPUT_DIR}/${SAMPLE_NAME}_split.bam"
RAW_VARIANTS="${OUTPUT_DIR}/${SAMPLE_NAME}_raw_variants.vcf"
FILTERED_VARIANTS="${OUTPUT_DIR}/${SAMPLE_NAME}_filtered_variants.vcf"
TMP_DIR="${OUTPUT_DIR}/temp_gatk"
mkdir -p "$TMP_DIR"

# =============================================================================
# Main Pipeline Steps
# =============================================================================

# Step 1: Create Sequence Dictionary
java -jar picard.jar CreateSequenceDictionary -R "$REF_FASTA" -O "$DICT"

# Step 2: Add Read Groups
java -jar picard.jar AddOrReplaceReadGroups \
    -I "$INPUT_BAM" \
    -O "$RG_BAM" \
    -RGID 1 \
    -RGLB lib1 \
    -RGPL illumina \
    -RGPU unit1 \
    -RGSM "$SAMPLE_NAME"

# Step 3: Mark Duplicates
java -jar picard.jar MarkDuplicates \
    -I "$RG_BAM" \
    -O "$DEDUP_BAM" \
    -M "${OUTPUT_DIR}/${SAMPLE_NAME}_dup_metrics.txt" \
    -REMOVE_DUPLICATES true

# Step 4: Split N Cigar Reads
gatk SplitNCigarReads \
    -R "$REF_FASTA" \
    -I "$DEDUP_BAM" \
    -O "$SPLIT_BAM" \
    --tmp-dir "$TMP_DIR"

# Step 5: HaplotypeCaller
gatk HaplotypeCaller \
    -R "$REF_FASTA" \
    -I "$SPLIT_BAM" \
    -O "$RAW_VARIANTS" \
    --native-pair-hmm-threads "$THREADS" \
    --minimum-mapping-quality 0 \
    --disable-read-filter MappingQualityReadFilter \
    --dont-use-soft-clipped-bases \
    --standard-min-confidence-threshold-for-calling 10 \
    --min-base-quality-score 10 \
    --max-reads-per-alignment-start 0 \
    --java-options "-XX:ParallelGCThreads=$THREADS -Dsamjdk.compression_level=2"

# Step 6: Variant Filtration
gatk VariantFiltration \
    -R "$REF_FASTA" \
    -V "$RAW_VARIANTS" \
    -O "$FILTERED_VARIANTS" \
    --filter-expression "FS > 30.0" \
    --filter-name "FS" \
    --filter-expression "QD < 2.0" \
    --filter-name "QD"

echo "Pipeline for $SAMPLE_NAME completed successfully!"