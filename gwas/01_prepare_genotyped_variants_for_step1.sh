#!/usr/bin/env bash
set -euo pipefail

PROJECT_DX="${PROJECT_DX:-${project:-}}"
if [ -z "${PROJECT_DX}" ]; then
  echo "Set PROJECT_DX to the UKB RAP project, for example project-xxxx"
  exit 1
fi

DX_WORK_DIR="${DX_WORK_DIR:-${PROJECT_DX}:/gwas_work}"
GENOTYPED_PREFIX="${GENOTYPED_PREFIX:-ukb22418_all_v2_merged}"
GENOTYPED_QC_PREFIX="${GENOTYPED_QC_PREFIX:-ukb22418_allQC_v2_merged_ppp_complement}"

MAF_THRESHOLD="${MAF_THRESHOLD:-0.01}"
VARIANT_MISSINGNESS_THRESHOLD="${VARIANT_MISSINGNESS_THRESHOLD:-0.10}"
HWE_P_THRESHOLD="${HWE_P_THRESHOLD:-1e-15}"
SAMPLE_MISSINGNESS_THRESHOLD="${SAMPLE_MISSINGNESS_THRESHOLD:-0.10}"

run_plink_qc="plink2 --bfile ${GENOTYPED_PREFIX} \
  --maf ${MAF_THRESHOLD} \
  --geno ${VARIANT_MISSINGNESS_THRESHOLD} \
  --hwe ${HWE_P_THRESHOLD} \
  --mind ${SAMPLE_MISSINGNESS_THRESHOLD} \
  --make-bed \
  --out ${GENOTYPED_QC_PREFIX}"

dx run swiss-army-knife \
  -iin="${DX_WORK_DIR}/${GENOTYPED_PREFIX}.bed" \
  -iin="${DX_WORK_DIR}/${GENOTYPED_PREFIX}.bim" \
  -iin="${DX_WORK_DIR}/${GENOTYPED_PREFIX}.fam" \
  -icmd="${run_plink_qc}" \
  --tag="GWAS step 1 genotype QC" \
  --instance-type "mem1_ssd1_v2_x16" \
  --destination="${DX_WORK_DIR}" \
  --brief --yes
