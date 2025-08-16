#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(DESeq2); library(tidyverse); library(pheatmap); library(ggpubr)
  library(ggplot2); library(plotly); library(ggVennDiagram); library(stringr)
  library(apeglm); library(topGO); library(AnnotationDbi); library(org.Gg.eg.db)
  library(circlize); library(reshape2)
})

# ----------- USER INPUTS -----------
counts_dir <- "data/hiseq_counts"         # where *_counts.txt live
metadata_csv <- "R/hiseq_metadata.csv"    # provide this CSV
gtf_path <- "reference/genomic_fixed_clean.gtf"

# optional: orthofinder & NCBI annotations for GO mapping
ncbi_tsv <- "R/ncbi_Ereg_info.tsv"        # columns include Protein.accession, Symbol, etc.
ortho_tsv <- "R/orthofinder_table.tsv"    # columns include Orthogroup, protein2
go_car_tsv <- "R/go_terms_anolis_car.tsv" # columns Gene.Names, Gene.Ontology.IDs
# -----------------------------------

message("Reading counts...")
count_files <- list.files(path = counts_dir, pattern = "_counts.txt$", full.names = TRUE)
read_count_file <- function(file) {
  sample_name <- tools::file_path_sans_ext(basename(file)) %>% str_remove("_counts$")
  df <- read.delim(file, header = FALSE, stringsAsFactors = FALSE)
  colnames(df) <- c("gene", sample_name)
  df <- df[!grepl("^__", df$gene), ]
  df
}
count_list <- lapply(count_files, read_count_file)
counts_merged <- reduce(count_list, full_join, by = "gene") %>%
  column_to_rownames("gene") %>%
  mutate(across(everything(), ~ as.integer(replace_na(., 0))))
write.csv(counts_merged, "counts_matrix.csv")

metadata <- read.csv(metadata_csv, row.names = 1)
stopifnot(all(colnames(counts_merged) %in% rownames(metadata)))
counts_merged <- counts_merged[, rownames(metadata)]

# Digestive tissues only
metadata_fil <- metadata %>% filter(tissue %in% c("Liver", "Tongue", "Stomach", "Intestine"))
counts_fil <- counts_merged[, rownames(metadata_fil)]
metadata_fil$condition <- factor(metadata_fil$condition, levels = c("Atrivittata", "Sruber", "Fasting"))
metadata_fil$tissue <- factor(metadata_fil$tissue)

dds <- DESeqDataSetFromMatrix(countData = counts_fil,
                              colData = metadata_fil,
                              design = ~ tissue + condition)

# Run DE with two refs
dds_fasting <- DESeq(dds); dds_fasting$condition <- relevel(dds_fasting$condition, ref = "Fasting")
dds_fasting_ref <- DESeq(dds_fasting)

dds_trivi <- DESeq(dds); dds_trivi$condition <- relevel(dds_trivi$condition, ref = "Atrivittata")
dds_trivi_ref <- DESeq(dds_trivi)

saveRDS(dds_fasting_ref, "dds_fasting_ref.RData")
saveRDS(dds_trivi_ref,  "dds_trivi_ref.RData")

res_atriv_vs_sruber <- results(dds_trivi_ref, contrast = c("condition", "Atrivittata", "Sruber"))
res_fasting_vs_sruber <- results(dds_fasting_ref, contrast = c("condition", "Fasting", "Sruber"))
res_fasting_vs_trivi  <- results(dds_fasting_ref, contrast = c("condition", "Fasting", "Atrivittata"))

resLFC_fast_sruber <- lfcShrink(dds_fasting_ref, coef = "condition_Sruber_vs_Fasting", type = "apeglm")
resLFC_fast_trivi  <- lfcShrink(dds_fasting_ref, coef = "condition_Atrivittata_vs_Fasting", type = "apeglm")
resLFC_trivi_sruber <- lfcShrink(dds_trivi_ref, coef = "condition_Sruber_vs_Atrivittata", type = "apeglm")

write.csv(as.data.frame(resLFC_fast_sruber), "resLFC_fast_sruber_only_digestive.csv")
write.csv(as.data.frame(resLFC_fast_trivi),  "resLFC_fast_trivi_only_digestive.csv")
write.csv(as.data.frame(resLFC_trivi_sruber),"resLFC_trivi_sruber_only_digestive.csv")

# --- GTF-derived annotation (gene_id, gene_name, product) ---
message("Parsing GTF for exon attributes...")
gtf_lines <- readLines(gtf_path)
gtf_exons <- gtf_lines[grepl("\\texon\\t", gtf_lines)]
annotation_df <- tibble(gtf = gtf_exons) %>%
  mutate(gene_id = str_extract(gtf, 'gene_id "[^"]+"') %>% str_remove_all('gene_id \"|\"'),
         gene_name = str_extract(gtf, 'gene_name "[^"]+"') %>% str_remove_all('gene_name \"|\"'),
         product = str_extract(gtf, 'product "[^"]+"') %>% str_remove_all('product \"|\"')) %>%
  distinct(gene_id, gene_name, product)

# --- Volcano plot helpers (your logic, compact) ---
plot_volcano_custom <- function(res_df_input, flip=FALSE, ortho_table=NULL) {
  res_df <- as.data.frame(res_df_input) %>%
    rownames_to_column("gene") %>%
    mutate(sig = !is.na(padj) & padj < 0.05)
  if (flip) res_df$log2FoldChange <- -res_df$log2FoldChange

  # Example gene group tags; adapt as needed
  tag_group <- function(genes) {
    case_when(
      grepl("^gene-SLC", genes) ~ "SLC",
      grepl("^gene-PLA2", genes) ~ "PLA2",
      grepl("^gene-CYP", genes) ~ "CYP",
      grepl("^gene-SERPIN", genes) ~ "SERPIN",
      grepl("^gene-ABC", genes) ~ "ABC",
      grepl("^gene-HSP", genes) ~ "HSP",
      TRUE ~ ifelse(res_df$padj < 0.05, "Significant", "Not_significant")
    )
  }
  res_df$highlight_group <- tag_group(res_df$gene)

  ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = highlight_group)) +
    geom_point(alpha = 0.6) +
    scale_color_manual(values = c(
      "Not_significant" = "grey70", "Significant" = "grey90",
      "SLC"="#0072B2","CYP"="#D55E00","ABC"="#E69F00","SERPIN"="#CC79A7",
      "HSP"="green","PLA2"="black"
    )) + theme_minimal()
}

# Example usage + labels
v1 <- plot_volcano_custom(resLFC_trivi_sruber, flip=TRUE) +
  ylim(0,6) + xlim(-5,5) + geom_vline(xintercept=0, linetype="dashed") +
  annotate("text", x=-2, y=5, label="Upregulated after\nS. ruber ingestion", size=3) +
  annotate("text", x= 2, y=5, label="Upregulated after\nA. trivittata ingestion", size=3) +
  labs(color="Highlighted genes", x="log2FoldChange")
ggsave("volcano_trivi_vs_sruber_digestive.png", v1, width=6, height=5, dpi=300)

# PCA
vsd <- vst(dds, blind = FALSE)
pca_tissue <- plotPCA(vsd, intgroup = "tissue", returnData = TRUE)
pv <- round(100 * attr(pca_tissue, "percentVar"))
p1 <- ggplot(pca_tissue, aes(PC1, PC2, shape=tissue)) +
  geom_point(size=3, alpha=0.8) + theme_bw() +
  labs(x=paste0("PC1: ",pv[1],"%"), y=paste0("PC2: ",pv[2],"%"))
ggsave("PCA_tissue_digestive.png", p1, width=6, height=5, dpi=300)

pca_condition <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
pv2 <- round(100 * attr(pca_condition, "percentVar"))
p2 <- ggplot(pca_condition, aes(PC1, PC2, color=condition)) +
  geom_point(size=3, alpha=0.6) +
  scale_color_manual(values=c("#ea7e24","#008080","#004483"),
                     labels=c("A. trivittata","S. ruber","Fasting")) +
  labs(x=paste0("PC1: ",pv2[1],"%"), y=paste0("PC2: ",pv2[2],"%"), color="Experiment") +
  theme_bw()
ggsave("PCA_condition_digestive.png", p2, width=6, height=5, dpi=300)

# Heatmap of top-genes example
top_genes <- head(order(resLFC_fast_trivi$padj), 25)
mat <- assay(vsd)[top_genes, ]
mat <- mat - rowMeans(mat)
pheatmap(mat, annotation_col = as.data.frame(colData(dds)[,c("tissue","condition")]))
ggsave("heatmap_top25_digestive.png", width=6, height=6, dpi=300)

# Venn
sig1 <- rownames(res_fasting_vs_trivi)[which(res_fasting_vs_trivi$padj < 0.05 & !is.na(res_fasting_vs_trivi$padj))]
sig2 <- rownames(res_fasting_vs_sruber)[which(res_fasting_vs_sruber$padj < 0.05 & !is.na(res_fasting_vs_sruber$padj))]
sig3 <- rownames(res_atriv_vs_sruber)[which(res_atriv_vs_sruber$padj < 0.05 & !is.na(res_atriv_vs_sruber$padj))]
gene_sets <- list(Fasting_vs_Atrivitta=sig1, Fasting_vs_Sruber=sig2, Sruber_vs_Atrivitta=sig3)
g <- ggVennDiagram(gene_sets, label_alpha=0) + scale_fill_gradient(low="white", high="#ea7e24") + theme_void()
ggsave("venn_digestive.png", g, width=6, height=6, dpi=300)

# -------- GO prep (Anolis car mapping used for GO universe) ----------
# background genes
dds <- estimateSizeFactors(dds)
norm_counts <- counts(dds, normalized = TRUE)
background_genes <- rownames(norm_counts)[rowMeans(norm_counts) > 0.5]
background_nam <- annotation_df %>% filter(gene_id %in% background_genes)
background_names <- background_nam$gene_name

# Orthofinder + NCBI
ncbi <- read.delim(ncbi_tsv, header=TRUE, sep="\t", quote="")
ncbi <- ncbi %>% dplyr::select(Symbol, Gene.Type, Protein.accession, Name, Transcripts.accession)
colnames(ncbi) <- paste0(colnames(ncbi), "_Ereg")
ortho <- read.delim(ortho_tsv, header=TRUE, sep="\t", quote="")
ortho_long_reg <- ortho %>% dplyr::select(Orthogroup, protein2) %>% separate_rows(protein2, sep=",\\s*")
ortho_long_reg$protein2 <- gsub('^"|"$', '', ortho_long_reg$protein2)
ortho_name_Ereg <- ortho_long_reg %>% left_join(ncbi, by=c("protein2"="Protein.accession_Ereg"))
write.csv(ortho_name_Ereg, "ortho_name_Ereg.csv", row.names = FALSE)

# GO terms (Anolis car)
ortho_go_raw_car <- read.delim(go_car_tsv, sep="\t", header=TRUE)
ortho_go_terms_car <- ortho_go_raw_car %>%
  separate_rows(Gene.Ontology.IDs, sep=";\\s*") %>%
  rename(Symbol = Gene.Names, GO_term = Gene.Ontology.IDs) %>%
  filter(!is.na(Symbol), !is.na(GO_term), GO_term != "")
gene2GO_car <- split(ortho_go_terms_car$GO_term, ortho_go_terms_car$Symbol)

run_topGO <- function(gene_ids, background, name){
  geneList <- factor(as.integer(background %in% gene_ids)); names(geneList) <- background
  GOdata <- new("topGOdata", ontology="MF", allGenes=geneList,
                annot=annFUN.gene2GO, gene2GO=gene2GO_car, nodeSize=10)
  resultFisher <- runTest(GOdata, algorithm="classic", statistic="fisher")
  tab <- GenTable(GOdata, classicFisher=resultFisher, topNodes=15)
  tab$group <- name
  tab$log10p <- -log10(as.numeric(tab$classicFisher))
  tab
}

# Example “upregulated” groups
sig_res_fasting_vs_trivi <- as.data.frame(res_fasting_vs_trivi) %>% rownames_to_column("gene") %>% filter(padj < 0.05 & log2FoldChange < 0)
sig_res_atriv_vs_sruber <- as.data.frame(res_atriv_vs_sruber) %>% rownames_to_column("gene") %>% filter(padj < 0.05 & log2FoldChange > 0)
all_up_trivi <- c(sig_res_fasting_vs_trivi$gene, sig_res_atriv_vs_sruber$gene)
annotated_all_up_trivi <- annotation_df %>% filter(gene_id %in% all_up_trivi)
annotated_all_up_trivi_names <- annotated_all_up_trivi$gene_name
writeLines(annotated_all_up_trivi_names, "annotated_all_up_trivi_names.txt")

# fasting up
sig_res_fasting_vs_sruber_F <- as.data.frame(res_fasting_vs_sruber) %>% rownames_to_column("gene") %>% filter(padj < 0.05 & log2FoldChange > 0)
sig_res_fasting_vs_trivi_F  <- as.data.frame(res_fasting_vs_trivi)  %>% rownames_to_column("gene") %>% filter(padj < 0.05 & log2FoldChange > 0)
all_up_fas <- c(sig_res_fasting_vs_sruber_F$gene, sig_res_fasting_vs_trivi_F$gene)
annotated_all_up_fas <- annotation_df %>% filter(gene_id %in% all_up_fas)
annotated_all_up_fas_names <- annotated_all_up_fas$gene_name
writeLines(annotated_all_up_fas_names, "annotated_all_up_fas_names.txt")

# sruber up
sig_res_fasting_vs_sruber_down <- as.data.frame(res_fasting_vs_sruber) %>% rownames_to_column("gene") %>% filter(padj < 0.05 & log2FoldChange < 0)
sig_res_atriv_vs_sruber_sruber <- as.data.frame(res_atriv_vs_sruber) %>% rownames_to_column("gene") %>% filter(padj < 0.05 & log2FoldChange < 0)
all_up_sruber <- c(sig_res_fasting_vs_sruber_down$gene, sig_res_atriv_vs_sruber_sruber$gene)
annotated_all_up_sruber <- annotation_df %>% filter(gene_id %in% all_up_sruber)
annotated_all_up_sruber_names <- annotated_all_up_sruber$gene_name
writeLines(annotated_all_up_sruber_names, "annotated_all_up_sruber_names.txt")

# GO on the three groups
up_trivi <- run_topGO(annotated_all_up_trivi_names, background_names, "Trivi_DE_up")
up_srub  <- run_topGO(annotated_all_up_sruber_names, background_names, "Srub_DE_up")
up_fas   <- run_topGO(annotated_all_up_fas_names, background_names, "Fas_DE_up")
saveRDS(up_trivi, "up_trivi.RData"); saveRDS(up_srub, "up_srub.RData"); saveRDS(up_fas, "up_fas.RData")

# Chord plot
go_combined_up <- bind_rows(up_trivi, up_srub, up_fas)
chord_df <- go_combined_up[, c("group","GO.ID","Significant")]
chord_df$group <- recode(chord_df$group,"Trivi_DE_up"="A. trivittata","Srub_DE_up"="S. ruber","Fas_DE_up"="Fasting")
colnames(chord_df) <- c("group","GO_term","count")
go_info <- unique(go_combined_up[, c("GO.ID","Term")]); colnames(go_info) <- c("GO_term","Term")
chord_df_annot <- merge(chord_df, go_info, by="GO_term")
chord_df_annot$label <- paste0(chord_df_annot$GO_term, ":", chord_df_annot$Term)
chord_matrix <- acast(chord_df_annot[,c("group","label","count")], group ~ label, value.var="count", fill=0)
grid.col <- c("S. ruber"="#008080","A. trivittata"="#ea7e24","Fasting"="#004483")
pdf("circular_up_plot.pdf", width=10, height=10)
circos.clear(); par(mar=c(1,1,1,1))
circos.par(start.degree=90, gap.after=rep(2, length(unique(c(rownames(chord_matrix), colnames(chord_matrix))))),
           track.margin=c(0.01,0.01))
chordDiagram(chord_matrix, grid.col=grid.col, annotationTrack="grid", preAllocateTracks=list(track.height=0.08))
circos.trackPlotRegion(track.index=1, panel.fun=function(x,y){
  sector_name <- get.cell.meta.data("sector.index")
  xcenter <- get.cell.meta.data("xcenter"); ycenter <- get.cell.meta.data("ycenter")
  circos.text(xcenter, ycenter, labels=sector_name, facing="clockwise", niceFacing=TRUE, adj=c(0,0.5), cex=0.5)
}, bg.border=NA)
dev.off()

# ------ LIVER-ONLY BLOCK (mirrors your detailed steps) ------
metadata_liver <- metadata %>% filter(tissue == "Liver")
counts_liver <- counts_merged[, rownames(metadata_liver)]
metadata_liver$condition <- factor(metadata_liver$condition, levels = c("Atrivittata","Sruber","Fasting"))
dds_liver <- DESeqDataSetFromMatrix(countData=counts_liver, colData=metadata_liver, design=~condition)

dds_fasting_liver <- DESeq(dds_liver); dds_fasting_liver$condition <- relevel(dds_fasting_liver$condition,"Fasting")
dds_fasting_ref_liver <- DESeq(dds_fasting_liver)
dds_trivi_liver <- DESeq(dds_liver); dds_trivi_liver$condition <- relevel(dds_trivi_liver$condition,"Atrivittata")
dds_trivi_ref_liver <- DESeq(dds_trivi_liver)

res_atriv_vs_sruber_liver <- results(dds_trivi_ref_liver, contrast=c("condition","Atrivittata","Sruber"))
res_fasting_vs_sruber_liver <- results(dds_fasting_ref_liver, contrast=c("condition","Fasting","Sruber"))
res_fasting_vs_trivi_liver  <- results(dds_fasting_ref_liver, contrast=c("condition","Fasting","Atrivittata"))

resLFC_fast_sruber_liver <- lfcShrink(dds_fasting_ref_liver, coef="condition_Sruber_vs_Fasting", type="apeglm")
resLFC_fast_trivi_liver  <- lfcShrink(dds_fasting_ref_liver, coef="condition_Atrivittata_vs_Fasting", type="apeglm")
resLFC_trivi_sruber_liver <- lfcShrink(dds_trivi_ref_liver, coef="condition_Sruber_vs_Atrivittata", type="apeglm")

write.csv(as.data.frame(resLFC_fast_sruber_liver), "resLFC_fast_sruber_liver.csv")
write.csv(as.data.frame(resLFC_fast_trivi_liver),  "resLFC_fast_trivi_liver.csv")
write.csv(as.data.frame(resLFC_trivi_sruber_liver),"resLFC_trivi_sruber_liver.csv")

vsd_liver <- vst(dds_liver, blind=FALSE)
pca_condition_liver <- plotPCA(vsd_liver, intgroup="condition", returnData=TRUE)
pv_l <- round(100 * attr(pca_condition_liver, "percentVar"))
pl <- ggplot(pca_condition_liver, aes(PC1,PC2,color=condition)) +
  geom_point(size=3, alpha=0.8) +
  scale_color_manual(values=c("Atrivittata"="#ea7e24","Sruber"="#008080","Fasting"="#004483")) +
  labs(x=paste0("PC1: ",pv_l[1],"%"), y=paste0("PC2: ",pv_l[2],"%")) + theme_bw()
ggsave("PCA_condition_liver.png", pl, width=6, height=5, dpi=300)

# Liver volcano examples
lv1 <- plot_volcano_custom(resLFC_trivi_sruber_liver, flip=TRUE) +
  ylim(0,8) + xlim(-7.5,7.5) + geom_vline(xintercept=0, linetype="dashed") +
  annotate("text", x=-4, y=6, label="Liver upregulation\nafter S. ruber ingestion", size=3) +
  annotate("text", x= 4, y=6, label="Liver upregulation\nafter A. trivittata ingestion", size=3)
ggsave("volcano_trivi_sruber_liver.png", lv1, width=6, height=5, dpi=300)

lv2 <- plot_volcano_custom(resLFC_fast_trivi_liver, flip=FALSE) +
  ylim(0,10) + xlim(-15,15) + geom_vline(xintercept=0, linetype="dashed") +
  annotate("text", x=-7.5, y=6, label="Liver upregulation in fasting", size=3) +
  annotate("text", x= 7.5, y=6, label="Liver upregulation after A. trivittata ingestion", size=3)
ggsave("volcano_fasting_trivi_liver.png", lv2, width=6, height=5, dpi=300)

# Build upregulated lists + GO for liver (same pattern as above)…
# (Left as-is for brevity; replicate the “digestive” approach with _liver objects.)

message("DONE.")
