This repository contains resources and instructions for using DeepVariant—a deep learning-based variant caller—for RNAseq analysis. DeepVariant leverages a convolutional neural network to analyze pileup image tensors generated from aligned reads (BAM or CRAM files) and produces high-quality variant calls in VCF or gVCF format.

## What is DeepVariant?

DeepVariant is an innovative variant caller that uses deep learning to accurately identify genetic variants. It works by:
- Taking aligned sequencing reads (BAM/CRAM format).
- Converting read data into pileup image tensors.
- Classifying these tensors using a convolutional neural network.
- Outputting variant calls in standard VCF/gVCF formats.

This approach enables DeepVariant to achieve high precision, even in complex genomic regions.

## Why Use DeepVariant for RNAseq Analysis?

While RNAseq is primarily used to measure gene expression, it can also reveal valuable information about RNA editing events and expressed genetic variants. DeepVariant’s advanced deep learning algorithms help to:
- Address challenges like variable coverage and alternative splicing inherent to RNAseq data.
- Improve the accuracy of variant detection within transcriptomic data.
- Provide a complementary approach to traditional variant callers in the context of RNAseq.

## What is GATK?

The Genome Analysis Toolkit (GATK) is a comprehensive software package developed by the Broad Institute for analyzing high-throughput sequencing data. It is widely used for tasks such as:
- Variant discovery and genotyping.
- Data pre-processing and quality control.
- Providing robust, standardized workflows for genomic data analysis.

GATK is often integrated with other tools to form end-to-end pipelines for variant calling.

## Requirements

- **JAVA:**  
- **Docker:**  
- **Bash:** 

## Installing GATK

Download GATK from the [GATK GitHub releases page](https://github.com/broadinstitute/gatk?tab=readme-ov-file#downloading).

## Installing Picard

Download Picard from the [Picard Tools website](https://broadinstitute.github.io/picard/).

## Steps to run the deepvariant

### Get the variant calling files from gatk


#### **run_gatk.sh**
```bash
bash pipeline.sh --sample D1nonkillers \
                 --input_bam /path/to/input.bam \
                 --output_dir /home/user/output_dir \
                 --ref_fasta GRCh38.primary_assembly.genome.fa \
                 --threads 16
                 
Result:

sample_gatk.vcf

```



- **run_deepvariant.sh**
```bash
bash deepvariant.sh --input_dir INPUT_DIR \
                    --output_dir OUTPUT_DIR \
                    --ref_dir REF_DIR \
                    --model_dir MODEL_DIR \
                    [--min_coverage MIN] \
                    [--annotation ANNOTATION_FILE]

result 
sample_deep_variant.vcf
```

- **happy.sh**
```bash 
bash happy.sh <truth_vcf> <query_vcf> <benchmark_bed> <reference_fasta> <target_regions>
```

### Result 
accuracy of deepvariant
Benchmarking summary:
Type Filter  TRUTH.TOTAL  TRUTH.TP  TRUTH.FN  QUERY.TOTAL  QUERY.FP  QUERY.UNK  FP.gt  FP.al  METRIC.Recall  METRIC.Precision  METRIC.Frac_NA  METRIC.F1_Score  TRUTH.TOTAL.TiTv_ratio  QUERY.TOTAL.TiTv_ratio  TRUTH.TOTAL.het_hom_ratio  QUERY.TOTAL.het_hom_ratio
INDEL    ALL        10628     10585        43        21064        22      10001     18      3       0.995954          0.998011        0.474791         0.996982                     NaN                     NaN                   1.748961                   2.319825
INDEL   PASS        10628     10585        43        21064        22      10001     18      3       0.995954          0.998011        0.474791         0.996982                     NaN                     NaN                   1.748961                   2.319825
  SNP    ALL        70166     69918       248        84834        56      14822     13      3       0.996466          0.999200        0.174718         0.997831                2.296566                2.083842                   1.883951                   1.913523
  SNP   PASS        70166     69918       248        84834        56      14822     13      3       0.996466          0.999200        0.174718         0.997831                2.296566                2.083842                   1.883951                   1.913523
