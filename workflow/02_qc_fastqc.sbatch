#!/usr/bin/env bash
#SBATCH --job-name=fastqc
#SBATCH --cpus-per-task=4
#SBATCH --mem=20G
#SBATCH --time=04:00:00
#SBATCH --output=%x_%j.log
#SBATCH --error=%x_%j.err

set -euo pipefail
source workflow/config.env

mkdir -p "${QC_DIR}" "${LOG_DIR}"
[[ -n "${MOD_FASTQC}" ]] && module load "${MOD_FASTQC}" || true

for fq1 in "${FASTQ_DIR}"/*_1.fastq; do
  fq2="${fq1/_1.fastq/_2.fastq}"
  [[ -f "${fq2}" ]] || continue
  fastqc -t "${THREADS}" -o "${QC_DIR}" "${fq1}" "${fq2}"
done
