#!/usr/bin/env bash
#SBATCH --job-name=hisat2
#SBATCH --cpus-per-task=4
#SBATCH --mem=40G
#SBATCH --time=12:00:00
#SBATCH --output=%x_%j.log
#SBATCH --error=%x_%j.err

set -euo pipefail
source workflow/config.env

mkdir -p "${BAM_DIR}" "${LOG_DIR}"
[[ -n "${MOD_HISAT2}" ]] && module load "${MOD_HISAT2}" || true
[[ -n "${MOD_SAMTOOLS}" ]] && module load "${MOD_SAMTOOLS}" || true

for fq1 in "${FASTQ_DIR}"/*_1.fastq; do
  fq2="${fq1/_1.fastq/_2.fastq}"
  [[ -f "${fq2}" ]] || continue

  base=$(basename "${fq1}" _1.fastq)
  sam="${BAM_DIR}/${base}.sam"
  bam="${BAM_DIR}/${base}.bam"
  sorted="${BAM_DIR}/${base}.sorted.bam"

  hisat2 -x "${HISAT2_INDEX_PREFIX}" -1 "${fq1}" -2 "${fq2}" -S "${sam}" -p "${THREADS}"
  samtools view -S -b "${sam}" > "${bam}"
  rm -f "${sam}"

  samtools sort "${bam}" -o "${sorted}"
  rm -f "${bam}"
  samtools index "${sorted}"
done
