#!/usr/bin/env bash
set -euo pipefail

trait="${1:-}"
case "${trait}" in
  fat_mass|fm)
    pheno_col="fat_mass"
    output_file="fat_mass_gwas_regenie.txt"
    ;;
  fat_free_mass|ffm)
    pheno_col="fat_free_mass"
    output_file="fat_free_mass_gwas_regenie.txt"
    ;;
  *)
    echo "Usage: bash gwas/05_merge_regenie_results.sh fat_mass|fat_free_mass"
    exit 1
    ;;
esac

PROJECT_DX="${PROJECT_DX:-${project:-}}"
if [ -z "${PROJECT_DX}" ]; then
  echo "Set PROJECT_DX to the UKB RAP project, for example project-xxxx"
  exit 1
fi

DX_GWAS_WORK_DIR="${DX_GWAS_WORK_DIR:-${PROJECT_DX}:/gwas_work/ppp_complement}"
DX_OUTPUT_DIR="${DX_OUTPUT_DIR:-${PROJECT_DX}:/controlled_access}"
GWAS_WORK_MOUNT="${GWAS_WORK_MOUNT:-/mnt/project/gwas_work/ppp_complement}"

merge_cmd="set -euo pipefail
out_file=\"${output_file}\"
cp \"${GWAS_WORK_MOUNT}\"/assoc_37_ppp_complement.c*_${pheno_col}.regenie.gz .
gunzip -f ./*.regenie.gz
first_file=\$(ls ./*.regenie | sort | head -n 1)
head -n 1 \"\${first_file}\" | tr -s ' ' '\t' > \"\${out_file}\"
for f in \$(ls ./*.regenie | sort); do
  tail -n +2 \"\${f}\" | tr -s ' ' '\t' >> \"\${out_file}\"
done"

dx run swiss-army-knife \
  -iin="${DX_GWAS_WORK_DIR}/assoc_37_ppp_complement.c1_${pheno_col}.regenie.gz" \
  -icmd="${merge_cmd}" \
  --tag="GWAS merge ${pheno_col}" \
  --instance-type "mem1_ssd1_v2_x16" \
  --destination="${DX_OUTPUT_DIR}" \
  --brief --yes
