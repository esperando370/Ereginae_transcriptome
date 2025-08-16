#!/usr/bin/env bash

cd reference

# 1) GFF → GTF
gffread ./genomic.gff -T -o ./genomic.gtf

# 2) Validate with reference
gffread ./genomic.gff -g "$REF_FASTA" -T -o ./genomic.gtf

# 3) Ensure exon features have gene_id; promote gene_name to gene_id if missing
cp genomic.gtf genomic_fixed.gtf

awk -F'\t' '{
  if ($3 == "exon" && $9 !~ /gene_id/) {
    match($9, /gene_name "([^"]+)"/, gname)
    if (gname[1] != "") {
      $9 = "gene_id \"" gname[1] "\"; " $9
    }
  }
  print
}' genomic_fixed.gtf > genomic_fixed.tmp && mv genomic_fixed.tmp genomic_fixed.gtf

# 4) Clean attribute formatting to keep 9 fields intact
awk 'BEGIN {OFS="\t"} {if (NF >= 9) print $1,$2,$3,$4,$5,$6,$7,$8,substr($0, index($0,$9));}' \
  genomic_fixed.gtf > genomic_fixed_clean.gtf

# Sanity checks
grep -P '\texon\t' genomic_fixed_clean.gtf | grep -v 'gene_id' | wc -l
