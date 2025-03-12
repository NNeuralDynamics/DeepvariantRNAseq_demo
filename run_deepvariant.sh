#!/bin/bash
set -euo pipefail

# =============================================================================
# Usage Information
# =============================================================================
usage() {
  echo "Usage: $0 --input_dir INPUT_DIR --output_dir OUTPUT_DIR --ref_dir REF_DIR --model_dir MODEL_DIR [--min_coverage MIN] [--annotation ANNOTATION_FILE]"
  echo ""
  echo "  --input_dir      Directory containing input BAM files (and .bai files)"
  echo "  --output_dir     Directory where output will be stored"
  echo "  --ref_dir        Directory containing reference files (e.g. GRCh38.primary_assembly.genome.fa and annotation)"
  echo "  --model_dir      Directory containing the DeepVariant model file (e.g. model.ckpt)"
  echo "  --min_coverage   (Optional) Minimum per-base coverage threshold [default: 3]"
  echo "  --annotation     (Optional) Annotation file (GTF format). Default: \$REF_DIR/gencode.v45.primary_assembly.annotation.gtf"
  exit 1
}

if [ "$#" -lt 8 ]; then
  usage
fi

# =============================================================================
# Command-line Argument Parsing
# =============================================================================
while [[ "$#" -gt 0 ]]; do
  case "$1" in
    --input_dir) INPUT_DIR="$2"; shift 2;;
    --output_dir) OUTPUT_DIR="$2"; shift 2;;
    --ref_dir) REF_DIR="$2"; shift 2;;
    --model_dir) MODEL_DIR="$2"; shift 2;;
    --min_coverage) MIN_COVERAGE="$2"; shift 2;;
    --annotation) ANNOTATION="$2"; shift 2;;
    *) echo "Error: Unknown parameter $1"; usage;;
  esac
done

# Set defaults if not provided
MIN_COVERAGE="${MIN_COVERAGE:-3}"
if [[ -z "${ANNOTATION:-}" ]]; then
  ANNOTATION="${REF_DIR}/gencode.v45.primary_assembly.annotation.gtf"
fi

# =============================================================================
# Configuration and Docker Image Definitions
# =============================================================================
BIN_VERSION="1.4.0"
DEEPVARIANT_IMAGE="google/deepvariant:${BIN_VERSION}"
BEDTOOLS_IMAGE="quay.io/biocontainers/bedtools:2.23.0--h5b5514e_6"
MOSDEPTH_IMAGE="quay.io/biocontainers/mosdepth:0.3.1--h4dc83fb_1"

# Calculate the number of CPUs to use (70% of available cores)
TOTAL_CPUS=$(nproc)
CPU_LIMIT=$(echo "${TOTAL_CPUS} * 0.70 / 1" | bc)

# Create main output directories
LOG_DIR="${OUTPUT_DIR}/logs"
INTERMEDIATE_RESULTS_TOP="${OUTPUT_DIR}/intermediate_results_dir"
mkdir -p "${OUTPUT_DIR}" "${LOG_DIR}" "${INTERMEDIATE_RESULTS_TOP}"

# =============================================================================
# Functions
# =============================================================================

# Check that required Docker images are available; pull if missing.
check_and_download_images() {
  echo "Checking required Docker images..."
  for IMAGE in "${DEEPVARIANT_IMAGE}" "${BEDTOOLS_IMAGE}" "${MOSDEPTH_IMAGE}"; do
    if [[ -z "$(docker images -q ${IMAGE} 2>/dev/null)" ]]; then
      echo "Image ${IMAGE} not found locally. Pulling..."
      docker pull "${IMAGE}"
    else
      echo "Image ${IMAGE} is available."
    fi
  done
}

# Get list of BAM files from INPUT_DIR and create per-sample output directories.
# For a file like D1_NKFeeder_pos.bam, the donor ID is "D1_NKFeeder" and the sample tag is "pos".
get_bam_files() {
  find "${INPUT_DIR}" -type f -name "*.bam" | while read -r BAM; do
    FILENAME=$(basename "${BAM}")
    # Extract donor ID (everything before the last underscore)
    DONOR_ID=$(echo "${FILENAME}" | rev | cut -d'_' -f2- | rev)
    # Determine sample tag based on filename content
    if [[ "${FILENAME}" == *"pos"* ]]; then
      SAMPLE_TAG="pos"
    else
      SAMPLE_TAG="neg"
    fi

    SAMPLE_OUTPUT_DIR="${OUTPUT_DIR}/${DONOR_ID}/${SAMPLE_TAG}"
    mkdir -p "${SAMPLE_OUTPUT_DIR}"
    echo "${BAM} ${SAMPLE_OUTPUT_DIR} ${DONOR_ID} ${SAMPLE_TAG}"
  done
}

# Run DeepVariant workflow per sample.
run_deepvariant() {
  check_and_download_images

  get_bam_files | while read -r BAM SAMPLE_OUTPUT_DIR DONOR_ID SAMPLE_TAG; do
    echo "========================================"
    echo "Processing sample: ${DONOR_ID} - ${SAMPLE_TAG}"
    echo "BAM file: ${BAM}"
    
    # Define file paths for this sample:
    INTERMEDIATE_DIR="${SAMPLE_OUTPUT_DIR}/intermediate_results"
    VCF_FILE="${SAMPLE_OUTPUT_DIR}/${DONOR_ID}_${SAMPLE_TAG}.vcf.gz"
    LOG_FILE="${SAMPLE_OUTPUT_DIR}/${DONOR_ID}_${SAMPLE_TAG}_deepvariant.log"
    CDS_BED="${SAMPLE_OUTPUT_DIR}/all_chromosomes_CDS.bed"
    COVERAGE_BED="${SAMPLE_OUTPUT_DIR}/${DONOR_ID}_${SAMPLE_TAG}_3x.bed"
    FINAL_BED="${SAMPLE_OUTPUT_DIR}/${DONOR_ID}_${SAMPLE_TAG}_final_regions.bed"

    mkdir -p "${INTERMEDIATE_DIR}"

    echo "[$(date)] Starting processing for ${DONOR_ID} - ${SAMPLE_TAG}" | tee "${LOG_FILE}"

    # ---------------------------
    # Step 1: Extract CDS Regions from GTF
    # ---------------------------
    echo "[$(date)] Extracting CDS regions from annotation file: ${ANNOTATION}" | tee -a "${LOG_FILE}"
    docker run --rm \
      -v "$(pwd):$(pwd)" \
      -w "$(pwd)" \
      "${BEDTOOLS_IMAGE}" \
      /bin/bash -c "awk -v OFS='\t' '\$3==\"CDS\" && \$4 < \$5 {print \$1, \$4-1, \$5}' '${ANNOTATION}' | grep -E '^chr' | sort -k1,1 -k2,2n | uniq > '${CDS_BED}'"

    # ---------------------------
    # Step 2: Calculate per-base coverage with mosdepth
    # ---------------------------
    echo "[$(date)] Calculating per-base coverage with mosdepth..." | tee -a "${LOG_FILE}"
    MOSDEPTH_PREFIX="${SAMPLE_OUTPUT_DIR}/${DONOR_ID}_${SAMPLE_TAG}_coverage"
    docker run --rm \
      -v "$(pwd):$(pwd)" \
      -w "$(pwd)" \
      "${MOSDEPTH_IMAGE}" \
      mosdepth --threads "${CPU_LIMIT}" "${MOSDEPTH_PREFIX}" "${BAM}"

    # ---------------------------
    # Step 3: Filter and merge regions with >= ${MIN_COVERAGE}x coverage
    # ---------------------------
    echo "[$(date)] Filtering regions with >= ${MIN_COVERAGE}x coverage..." | tee -a "${LOG_FILE}"
    docker run --rm \
      -v "$(pwd):$(pwd)" \
      -w "$(pwd)" \
      "${BEDTOOLS_IMAGE}" \
      /bin/bash -c "gzip -dc '${MOSDEPTH_PREFIX}.per-base.bed.gz' | egrep -v 'HLA|decoy|random|alt|chrUn|chrEBV' | awk -v OFS='\t' -v min_cov=${MIN_COVERAGE} '\$4>=min_cov {print \$1, \$2, \$3}' > '${SAMPLE_OUTPUT_DIR}/temp_3x.bed'"

    echo "[$(date)] Merging adjacent regions using bedtools..." | tee -a "${LOG_FILE}"
    docker run --rm \
      -v "$(pwd):$(pwd)" \
      -w "$(pwd)" \
      "${BEDTOOLS_IMAGE}" \
      /bin/bash -c "bedtools merge -d 1 -i '${SAMPLE_OUTPUT_DIR}/temp_3x.bed' > '${COVERAGE_BED}'"
    rm -f "${SAMPLE_OUTPUT_DIR}/temp_3x.bed"

    # ---------------------------
    # Step 4: Intersect CDS Regions with Coverage Regions
    # ---------------------------
    echo "[$(date)] Intersecting CDS regions with coverage regions..." | tee -a "${LOG_FILE}"
    docker run --rm \
      -v "$(pwd):$(pwd)" \
      -w "$(pwd)" \
      "${BEDTOOLS_IMAGE}" \
      /bin/bash -c "bedtools intersect -a '${COVERAGE_BED}' -b '${CDS_BED}' > '${FINAL_BED}'"

    # ---------------------------
    # Step 5: Run DeepVariant
    # ---------------------------
    echo "[$(date)] Running DeepVariant..." | tee -a "${LOG_FILE}"
    docker run --rm \
      -v "$(pwd):$(pwd)" \
      -w "$(pwd)" \
      "${DEEPVARIANT_IMAGE}" \
      run_deepvariant \
        --model_type=WES \
        --customized_model="${MODEL_DIR}/model.ckpt" \
        --ref="${REF_DIR}/GRCh38.primary_assembly.genome.fa" \
        --reads="${BAM}" \
        --output_vcf="${VCF_FILE}" \
        --num_shards="${CPU_LIMIT}" \
        --regions="${FINAL_BED}" \
        --make_examples_extra_args="split_skip_reads=true,channels=''" \
        --runtime_report \
        --logging_dir="${SAMPLE_OUTPUT_DIR}" \
        --intermediate_results_dir="${INTERMEDIATE_DIR}" 2>&1 | tee -a "${LOG_FILE}"

    echo "[$(date)] DeepVariant run completed for ${DONOR_ID} - ${SAMPLE_TAG}. VCF saved at ${VCF_FILE}" | tee -a "${LOG_FILE}"
  done
}

# =============================================================================
# Main Script Execution
# =============================================================================
run_deepvariant
