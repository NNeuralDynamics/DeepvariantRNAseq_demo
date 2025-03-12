#!/bin/bash

# Check if required arguments are provided
if [ "$#" -ne 6 ]; then
    echo "Usage: $0 <truth_vcf> <query_vcf> <benchmark_bed> <reference_fasta> <target_regions> <vcfeval_template>"
    exit 1
fi

# Assign variables from arguments
TRUTH_VCF=$1
QUERY_VCF=$2
BENCHMARK_BED=$3
REFERENCE_FASTA=$4
TARGET_REGIONS=$5
VCFEVAL_TEMPLATE=$6

# Run hap.py with user-provided inputs
docker run \
  -v $(pwd):$(pwd) \
  -w $(pwd) \
  jmcdani20/hap.py:v0.3.12 /opt/hap.py/bin/hap.py \
    "$TRUTH_VCF" \
    "$QUERY_VCF" \
    -f "$BENCHMARK_BED" \
    -r "$REFERENCE_FASTA" \
    -o happy/happy.output \
    --engine=vcfeval \
    --pass-only \
    --target-regions="$TARGET_REGIONS" \
    --threads=32 \
    --engine-vcfeval-template "$VCFEVAL_TEMPLATE"