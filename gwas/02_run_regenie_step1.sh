#!/usr/bin/env bash
set -euo pipefail

trait="${1:-}"
case "${trait}" in
  fat_mass|fm)
    pheno_col="fat_mass"
    pheno_file="n_fm_ppp_complement_4regenie.phe"
    out_prefix="fm_regenie_37_results_ppp_complement"
    ;;
  fat_free_mass|ffm)
    pheno_col="fat_free_mass"
    pheno_file="n_ffm_ppp_complement_4regenie.phe"
    out_prefix="ffm_regenie_37_results_ppp_complement"
    ;;
  *)
    echo "Usage: bash gwas/02_run_regenie_step1.sh fat_mass|fat_free_mass"
    exit 1
    ;;
esac

PROJECT_DX="${PROJECT_DX:-${project:-}}"
if [ -z "${PROJECT_DX}" ]; then
  echo "Set PROJECT_DX to the UKB RAP project, for example project-xxxx"
  exit 1
fi

DX_WORK_DIR="${DX_WORK_DIR:-${PROJECT_DX}:/gwas_work}"
DX_TEXT_DIR="${DX_TEXT_DIR:-${PROJECT_DX}:/controlled_access/gwas}"
GENOTYPED_QC_PREFIX="${GENOTYPED_QC_PREFIX:-ukb22418_allQC_v2_merged_ppp_complement}"

run_regenie_cmd="regenie --step 1 --out ${out_prefix} \
  --bed ${GENOTYPED_QC_PREFIX} \
  --phenoFile ${pheno_file} \
  --covarFile ${pheno_file} \
  --phenoCol ${pheno_col} \
  --covarCol sex \
  --covarCol age \
  --covarCol age2 \
  --covarCol age_sex \
  --covarCol age2_sex \
  --covarCol pc{1:20} \
  --catCovarList genotype_measurement_batch,centre,genetic_array \
  --bsize 1000 --loocv --gz --threads 16 --maxCatLevels 120"

dx run swiss-army-knife \
  -iin="${DX_WORK_DIR}/${GENOTYPED_QC_PREFIX}.bed" \
  -iin="${DX_WORK_DIR}/${GENOTYPED_QC_PREFIX}.bim" \
  -iin="${DX_WORK_DIR}/${GENOTYPED_QC_PREFIX}.fam" \
  -iin="${DX_TEXT_DIR}/${pheno_file}" \
  -icmd="${run_regenie_cmd}" \
  --tag="GWAS REGENIE step 1 ${pheno_col}" \
  --instance-type "mem2_ssd2_v2_x16" \
  --destination="${DX_WORK_DIR}" \
  --brief --yes
