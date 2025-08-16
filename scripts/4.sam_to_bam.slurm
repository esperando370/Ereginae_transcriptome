#!/usr/bin/env bash
#SBATCH -J bam_qc
#SBATCH -o data/aligned_data/%x_%j.out
#SBATCH -e data/aligned_data/%x_%j.err
#SBATCH -c 24
#SBATCH -t 08:00:00
#SBATCH -p savio

SAM="${ALN}/${sample}.sam"
BAM="${ALN}/${sample}.bam"
SORTBAM="${ALN}/${sample}.sorted.bam"

module load picard

# 1) SAM -> BAM
"$SAMTOOLS_BIN" view -@ 24 -bS "$SAM" -o "$BAM"

# 2) Sort
picard SortSam I="$BAM" O="$SORTBAM" SORT_ORDER=coordinate

# 3) Index
picard BuildBamIndex I="$SORTBAM"

# 4) QC
"$SAMTOOLS_BIN" flagstat "$BAM" > "${ALN}/${sample}_flagstat.txt"
