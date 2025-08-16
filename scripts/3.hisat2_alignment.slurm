#!/usr/bin/env bash
#SBATCH -J hisat2_aln
#SBATCH -o data/aligned_data/%x_%j.out
#SBATCH -e data/aligned_data/%x_%j.err
#SBATCH -c 24
#SBATCH -t 24:00:00
#SBATCH -p savio


R1="${CLEAN}/${sample}_R1.fastq.gz"
R2="${CLEAN}/${sample}_R2.fastq.gz"

module load hisat2

hisat2 -p 24 --dta -x "$HISAT2_INDEX" \
  -1 "$R1" -2 "$R2" \
  -S "${ALN}/${sample}.sam"
