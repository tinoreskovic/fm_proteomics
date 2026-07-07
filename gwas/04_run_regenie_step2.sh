#!/usr/bin/env bash
set -euo pipefail

trait="${1:-}"
case "${trait}" in
  fat_mass|fm)
    pheno_col="fat_mass"
    pheno_file="n_fm_ppp_complement_4regenie.phe"
    pred_prefix="fm_regenie_37_results_ppp_complement"
    ;;
  fat_free_mass|ffm)
    pheno_col="fat_free_mass"
    pheno_file="n_ffm_ppp_complement_4regenie.phe"
    pred_prefix="ffm_regenie_37_results_ppp_complement"
    ;;
  *)
    echo "Usage: bash gwas/04_run_regenie_step2.sh fat_mass|fat_free_mass"
    exit 1
    ;;
esac

PROJECT_DX="${PROJECT_DX:-${project:-}}"
if [ -z "${PROJECT_DX}" ]; then
  echo "Set PROJECT_DX to the UKB RAP project, for example project-xxxx"
  exit 1
fi

DX_TEXT_DIR="${DX_TEXT_DIR:-${PROJECT_DX}:/controlled_access/gwas}"
DX_WORK_DIR="${DX_WORK_DIR:-${PROJECT_DX}:/gwas_work}"
DX_GWAS_WORK_DIR="${DX_GWAS_WORK_DIR:-${PROJECT_DX}:/gwas_work/ppp_complement}"
IMPUTED_FIELD="${IMPUTED_FIELD:-ukb22828}"

run_regenie_for_chr() {
  local chr="$1"
  local chr_label="$2"
  local run_regenie_cmd

  run_regenie_cmd="regenie --step 2 --out assoc_37_ppp_complement.c${chr_label} \
    --bed ${IMPUTED_FIELD}_c${chr_label}_v3_ppp_complement \
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
    --pred ${pred_prefix}_pred.list \
    --bsize 200 --threads 16 --maxCatLevels 120 --gz"

  dx run swiss-army-knife \
    -iin="${DX_WORK_DIR}/${IMPUTED_FIELD}_c${chr_label}_v3_ppp_complement.bed" \
    -iin="${DX_WORK_DIR}/${IMPUTED_FIELD}_c${chr_label}_v3_ppp_complement.bim" \
    -iin="${DX_WORK_DIR}/${IMPUTED_FIELD}_c${chr_label}_v3_ppp_complement.fam" \
    -iin="${DX_TEXT_DIR}/${pheno_file}" \
    -iin="${DX_WORK_DIR}/${pred_prefix}_pred.list" \
    -iin="${DX_WORK_DIR}/${pred_prefix}_1.loco.gz" \
    -icmd="${run_regenie_cmd}" \
    --tag="GWAS REGENIE step 2 ${pheno_col} chr${chr_label}" \
    --instance-type "mem1_ssd1_v2_x16" \
    --destination="${DX_GWAS_WORK_DIR}" \
    --brief --yes
}

for chr in {1..22}; do
  run_regenie_for_chr "${chr}" "${chr}"
done

run_regenie_for_chr "X" "X"
