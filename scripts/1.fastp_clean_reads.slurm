#!/usr/bin/env bash
#SBATCH -J fastp
#SBATCH -o data/cleaned_data/%x_%j.out
#SBATCH -e data/cleaned_data/%x_%j.err
#SBATCH -c 8
#SBATCH -t 08:00:00
#SBATCH -p savio


# Try RAW2 first; if not found, fall back to RAW1
R1="${RAW2}/${sample}_R1.fastq.gz"
R2="${RAW2}/${sample}_R2.fastq.gz"
if [[ ! -f "$R1" || ! -f "$R2" ]]; then
  R1="${RAW1}/${sample}_R1.fastq.gz"
  R2="${RAW1}/${sample}_R2.fastq.gz"
fi

if [[ ! -f "$R1" || ! -f "$R2" ]]; then
  echo "ERROR: FASTQs not found for ${sample}" >&2
  exit 1
fi

module load fastp

fastp \
  -i "$R1" -I "$R2" \
  -o "${CLEAN}/${sample}_R1.fastq.gz" \
  -O "${CLEAN}/${sample}_R2.fastq.gz" \
  --detect_adapter_for_pe \
  -h "${CLEAN}/${sample}_fastp.html" \
  -j "${CLEAN}/${sample}_fastp.json"
