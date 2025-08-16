#!/usr/bin/env bash
#SBATCH -J hisat2_index
#SBATCH -o reference/%x_%j.out
#SBATCH -e reference/%x_%j.err
#SBATCH -c 24
#SBATCH -t 04:00:00
#SBATCH -p savio


module load hisat2 

hisat2-build -p 24 "$REF_FASTA" "$HISAT2_INDEX"
