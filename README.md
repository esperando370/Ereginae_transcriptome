# Ereginae_transcriptome
E. reginae transcriptome pipeline

# RNA-seq pipeline (fastp → HISAT2 → samtools/picard → htseq-count → DESeq2)

This repo contains SLURM job scripts and an R analysis to clean paired-end reads, build a HISAT2 index, align, generate counts with `htseq-count`, and run DESeq2 + figures and GO analyses.

## Requirements

- SLURM cluster
- Conda environment with:
  - `fastp`, `hisat2`, `samtools`, `picard`, `gffread`, `htseq`
  - R (≥4.2) with: `DESeq2`, `tidyverse`, `pheatmap`, `ggpubr`, `ggplot2`, `plotly`, `ggVennDiagram`, `stringr`, `apeglm`, `topGO`, `AnnotationDbi`, `org.Gg.eg.db`, `circlize`, `reshape2`
- Edit `config/paths.env` to match your user paths

## Quick start

```bash
# 0) Edit config/paths.env with your absolute paths
nano config/paths.env

# 1) Clean reads with fastp for all samples
for i in $(cat config/samples/all_samples_names.txt); do
  sbatch --export=sample=$i scripts/1.fastp_clean_reads.slurm
done

# 2) Build HISAT2 index (once)
sbatch scripts/2.hisat2_index.slurm

# 3) Align with HISAT2
for i in $(cat config/samples/all_samples_names.txt); do
  sbatch --export=sample=$i scripts/3.hisat2_alignment.slurm
done

# 4) SAM→BAM, sort+index, QC
for i in $(cat config/samples/all_samples_names.txt); do
  sbatch --export=sample=$i scripts/4.sam_to_bam.slurm
done

# 5) Prepare GTF (only once, when reference changes)
bash scripts/gff_to_gtf_and_fix.sh

# 6) htseq-count (per sample)
for i in $(cat config/samples/all_samples_names.txt); do
  sbatch --export=sample=$i scripts/5.htseq_counts.slurm
done

# 7) Run the R analysis locally or on the cluster login node
Rscript R/6_DESeq2_pipeline.R
