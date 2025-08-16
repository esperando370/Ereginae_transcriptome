#!/usr/bin/env bash
#SBATCH -J htseq
#SBATCH -o data/hiseq_counts/%x_%j.out
#SBATCH -e data/hiseq_counts/%x_%j.err
#SBATCH -c 8
#SBATCH -t 04:00:00
#SBATCH -p savio

module load python
module load htseq

GTF="reference/genomic_fixed_clean.gtf"
SORTBAM="${ALN}/${sample}.sorted.bam"

htseq-count \
  -f bam \
  -r name \
  -s no \
  -t exon \
  -i gene_id \
  "$SORTBAM" \
  "$GTF" > "${COUNTS}/${sample}_counts.txt"
