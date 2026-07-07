#!/usr/bin/env bash
set -euo pipefail

trait="${1:-}"
case "${trait}" in
  fat_mass|fm)
    pheno_file="n_fm_ppp_complement_4regenie.phe"
    ;;
  fat_free_mass|ffm)
    pheno_file="n_ffm_ppp_complement_4regenie.phe"
    ;;
  *)
    echo "Usage: bash gwas/03_prepare_imputed_variants_for_step2.sh fat_mass|fat_free_mass"
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
DX_IMPUTATION_DIR="${DX_IMPUTATION_DIR:-/Bulk/Imputation/UKB imputation from genotype}"
IMPUTED_FIELD="${IMPUTED_FIELD:-ukb22828}"

MAF_THRESHOLD="${MAF_THRESHOLD:-0.01}"
VARIANT_MISSINGNESS_THRESHOLD="${VARIANT_MISSINGNESS_THRESHOLD:-0.10}"
HWE_P_THRESHOLD="${HWE_P_THRESHOLD:-1e-15}"

run_plink_for_chr() {
  local chr="$1"
  local chr_label="$2"
  local run_plink_imp

  run_plink_imp="plink2 --bgen ${IMPUTED_FIELD}_c${chr}_b0_v3.bgen ref-first \
    --sample ${IMPUTED_FIELD}_c${chr}_b0_v3.sample \
    --no-pheno \
    --keep ${pheno_file} \
    --maf ${MAF_THRESHOLD} \
    --geno ${VARIANT_MISSINGNESS_THRESHOLD} \
    --hwe ${HWE_P_THRESHOLD} \
    --make-bed \
    --out ${IMPUTED_FIELD}_c${chr_label}_v3_ppp_complement"

  dx run swiss-army-knife \
    -iin="${DX_IMPUTATION_DIR}/${IMPUTED_FIELD}_c${chr}_b0_v3.bgen" \
    -iin="${DX_IMPUTATION_DIR}/${IMPUTED_FIELD}_c${chr}_b0_v3.sample" \
    -iin="${DX_TEXT_DIR}/${pheno_file}" \
    -icmd="${run_plink_imp}" \
    --tag="GWAS step 2 imputed QC chr${chr_label}" \
    --instance-type "mem2_ssd2_v2_x16" \
    --destination="${DX_WORK_DIR}" \
    --brief --yes
}

for chr in {1..22}; do
  run_plink_for_chr "${chr}" "${chr}"
done

run_plink_for_chr "X" "X"
