# =============================================================================
# Differential transcriptome analysis of the posterior superior temporal gyrus
# (STGp) in schizophrenia (SZ) vs healthy controls (HC)
#
# Input  (data/):    Metadata_rin_tin_age_sex.csv, HC_counts_an.txt, SZ_counts_an.txt
# Output (results/): figures (.png), tables (.csv), intermediate objects (.rds),
#                    sessionInfo.txt
#
# Run from the repository root:  source("analysis.R")
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(DESeq2)
  library(org.Hs.eg.db)
  library(pheatmap)
  library(ggrepel)
  library(patchwork)
  library(clusterProfiler)
})

data_dir <- "data"
out_dir  <- "results"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
out <- function(f) file.path(out_dir, f)

set.seed(42)


# 1. Load and merge data -------------------------------------------------------

meta <- read.csv(file.path(data_dir, "Metadata_rin_tin_age_sex.csv"), row.names = 1)
hc <- fread(file.path(data_dir, "HC_counts_an.txt"), header = TRUE, sep = "auto",
            data.table = FALSE)
sz <- fread(file.path(data_dir, "SZ_counts_an.txt"), header = TRUE, sep = "auto",
            data.table = FALSE)

counts_all <- cbind(hc, sz[, -1])
rownames(counts_all) <- counts_all$gene_id
counts_all$gene_id <- NULL


# 2. STGp subset and sample quality control -----------------------------------

stgp_meta <- meta[meta$region_acr == "STGp" & meta$group %in% c("SZ", "HC"), ]

# QC threshold: 5th percentile of TIN
tin_threshold <- quantile(stgp_meta$tin, probs = 0.05, na.rm = TRUE)

stgp_meta_before_qc <- stgp_meta %>%
  mutate(QC_status = ifelse(tin < tin_threshold, "Removed", "Retained"))

# Relationship between RIN and TIN (Spearman correlation)
cor_result <- cor.test(stgp_meta_before_qc$rin, stgp_meta_before_qc$tin,
                       method = "spearman")
cor_label <- paste0("Spearman R = ", round(cor_result$estimate, 2),
                    ", p = ", format.pval(cor_result$p.value, digits = 3))

ggplot(stgp_meta_before_qc,
       aes(x = rin, y = tin, color = group, shape = QC_status)) +
  geom_point(size = 3) +
  geom_hline(yintercept = tin_threshold, linetype = "dashed", color = "red") +
  geom_smooth(method = "lm", se = TRUE, color = "grey40", linetype = "solid",
              inherit.aes = FALSE, aes(x = rin, y = tin)) +
  annotate("text", x = max(stgp_meta_before_qc$rin, na.rm = TRUE) * 0.6,
           y = min(stgp_meta_before_qc$tin, na.rm = TRUE) + 1,
           label = cor_label, size = 3.5, hjust = 0) +
  theme_bw() +
  labs(title = "Relationship between RIN and TIN in STGp samples",
       subtitle = "Dashed line = 5th percentile TIN threshold",
       x = "RIN", y = "TIN", color = "Group", shape = "QC status")

ggsave(out("STGp_RIN_TIN_scatter.png"), width = 7, height = 5, dpi = 150)

# Remove samples below the TIN threshold
stgp_meta <- stgp_meta[!is.na(stgp_meta$tin) & stgp_meta$tin >= tin_threshold, ]
cat("Samples removed by TIN filter:",
    sum(stgp_meta_before_qc$QC_status == "Removed"), "\n")
counts_stgp <- counts_all[, rownames(stgp_meta), drop = FALSE]

# TIN and RIN distributions after QC
ggplot(stgp_meta, aes(x = group, y = tin, fill = group)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, alpha = 0.6) +
  labs(title = "TIN distribution in STGp samples after QC",
       y = "TIN (Transcript Integrity Number)") +
  theme_bw()
ggsave(out("STGp_TIN_distribution.png"), width = 6, height = 5)

ggplot(stgp_meta, aes(x = group, y = rin, fill = group)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, alpha = 0.6) +
  labs(title = "RIN distribution in STGp samples after QC",
       y = "RIN (RNA Integrity Number)") +
  theme_bw()
ggsave(out("STGp_RIN_distribution.png"), width = 6, height = 5)

# Covariate balance in the final sample (Welch t-tests, sex table)
print(t.test(age ~ group, data = stgp_meta))
print(t.test(rin ~ group, data = stgp_meta))
print(t.test(tin ~ group, data = stgp_meta))
print(table(stgp_meta$sex, stgp_meta$group))


# 3. Protein-coding gene selection --------------------------------------------

clean_ids <- gsub("\\..*", "", rownames(counts_stgp))

gene_info <- AnnotationDbi::select(org.Hs.eg.db, keys = clean_ids,
                    columns = "GENETYPE", keytype = "ENSEMBL")

pc_genes   <- gene_info$ENSEMBL[gene_info$GENETYPE == "protein-coding"]
counts_pc  <- counts_stgp[clean_ids %in% pc_genes, ]

cat("Protein-coding genes:", nrow(counts_pc), "\n")


# 4. Expression filter: detected in >= 50% of samples in either group ---------

sz_samples <- rownames(stgp_meta)[stgp_meta$group == "SZ"]
hc_samples <- rownames(stgp_meta)[stgp_meta$group == "HC"]

keep_sz <- rowSums(counts_pc[, sz_samples, drop = FALSE] > 0) >= 0.5 * length(sz_samples)
keep_hc <- rowSums(counts_pc[, hc_samples, drop = FALSE] > 0) >= 0.5 * length(hc_samples)

counts_filtered <- counts_pc[keep_sz | keep_hc, ]

cat("Genes before filtering:", nrow(counts_pc), "\n")
cat("Genes after filtering (>= 50% of samples in SZ or HC):", nrow(counts_filtered), "\n")

saveRDS(counts_filtered, out("STGp_counts_filtered.rds"))
saveRDS(stgp_meta,       out("STGp_metadata_filtered.rds"))
print(summary(stgp_meta$tin))


# 5. Normalization and modelling (DESeq2) -------------------------------------

counts_filtered <- counts_filtered[, rownames(stgp_meta)]

stgp_meta$group <- factor(stgp_meta$group, levels = c("HC", "SZ"))
stgp_meta$sex   <- factor(stgp_meta$sex)

# Design: sex, age and TIN as covariates; group is the effect of interest
dds <- DESeqDataSetFromMatrix(countData = counts_filtered,
                              colData   = stgp_meta,
                              design    = ~ sex + age + tin + group)

# Pre-filter: at least 10 reads in total across samples
dds <- dds[rowSums(counts(dds)) >= 10, ]
dds <- DESeq(dds)
vsd <- vst(dds, blind = FALSE)

res         <- results(dds, contrast = c("group", "SZ", "HC"))
res_ordered <- res[order(res$padj), ]
res_df      <- as.data.frame(res_ordered)

write.csv(res_df, out("STGp_Wald_results_all_genes.csv"))


# 6. PCA: VST, rlog and log2(normalized counts + 1) ---------------------------

p_vst <- plotPCA(vsd, intgroup = "group") +
  theme_minimal() +
  labs(title = "PCA on VST-transformed data")
ggsave(out("STGp_PCA_VST.png"), p_vst, width = 6, height = 5, dpi = 150)

rld <- rlog(dds, blind = FALSE)
p_rlog <- plotPCA(rld, intgroup = "group") +
  theme_minimal() +
  labs(title = "PCA on rlog-transformed counts")
ggsave(out("STGp_PCA_rlog.png"), p_rlog, width = 6, height = 5, dpi = 150)

log_counts <- log2(counts(dds, normalized = TRUE) + 1)
pca <- prcomp(t(log_counts), center = TRUE, scale. = FALSE)

pca_df <- data.frame(PC1 = pca$x[, 1],
                     PC2 = pca$x[, 2],
                     group = colData(dds)$group)

percentVar <- pca$sdev^2 / sum(pca$sdev^2)
pc1_lab <- paste0("PC1: ", round(percentVar[1] * 100, 1), "%")
pc2_lab <- paste0("PC2: ", round(percentVar[2] * 100, 1), "%")

p_log <- ggplot(pca_df, aes(PC1, PC2, color = group)) +
  geom_point(size = 3) +
  theme_minimal() +
  labs(title = "PCA on log2(normalized counts + 1)", x = pc1_lab, y = pc2_lab)
ggsave(out("STGp_PCA_log2norm.png"), p_log, width = 6, height = 5, dpi = 150)


# 7. Differential expression: Wald test and LRT -------------------------------

# 7.1 Wald test
res_sig_only <- res_df[which(res_df$padj < 0.05), ]
print(res_sig_only)
res_wald_sig <- as.character(rownames(res_sig_only))
cat("Significant genes (Wald):", length(res_wald_sig), "\n")

# 7.2 Likelihood ratio test (reduced model without the group term)
dds_lrt <- DESeq(dds, test = "LRT", reduced = ~ sex + age + tin)

res_lrt     <- results(dds_lrt)
res_lrt_df  <- as.data.frame(res_lrt[order(res_lrt$padj), ])
res_lrt_sig <- as.character(rownames(res_lrt_df[which(res_lrt_df$padj < 0.05), ]))
cat("Significant genes (LRT):", length(res_lrt_sig), "\n")

# 7.3 Comparison of the two tests
both_sig  <- base::intersect(res_wald_sig, res_lrt_sig)
only_wald <- base::setdiff(res_wald_sig, res_lrt_sig)
only_lrt  <- base::setdiff(res_lrt_sig, res_wald_sig)

cat("Significant in both tests:", length(both_sig), "\n")
cat("Wald only:", length(only_wald), "\n")
cat("LRT only:",  length(only_lrt),  "\n")
print(summary(res_lrt_df$pvalue))

# 7.4 Significant LRT genes
lrt_sig_df <- res_lrt_df[which(res_lrt_df$padj < 0.05), ]
lrt_sig_df <- lrt_sig_df[order(lrt_sig_df$padj), ]

lrt_table <- data.frame(
  ENSEMBL_ID     = rownames(lrt_sig_df),
  baseMean       = round(lrt_sig_df$baseMean, 2),
  log2FoldChange = round(lrt_sig_df$log2FoldChange, 3),
  padj           = formatC(lrt_sig_df$padj, format = "e", digits = 2),
  row.names      = NULL
)
print(lrt_table)
write.csv(lrt_table, out("STGp_LRT_significant_genes.csv"), row.names = FALSE)

# 7.5 Filter out technical artefacts: |log2FC| > 0.5 and baseMean < 100000
lrt_sig_filtered <- lrt_sig_df[abs(lrt_sig_df$log2FoldChange) > 0.5 &
                                 lrt_sig_df$baseMean < 100000, ]
cat("Genes after filtering:", nrow(lrt_sig_filtered), "\n")

lrt_table_filtered <- data.frame(
  ENSEMBL_ID     = rownames(lrt_sig_filtered),
  baseMean       = round(lrt_sig_filtered$baseMean, 2),
  log2FoldChange = round(lrt_sig_filtered$log2FoldChange, 3),
  padj           = formatC(lrt_sig_filtered$padj, format = "e", digits = 2),
  row.names      = NULL
)
print(lrt_table_filtered)
write.csv(lrt_table_filtered, out("STGp_LRT_significant_genes_filtered.csv"),
          row.names = FALSE)

# 7.6 Gene symbols
clean_ids <- gsub("\\..*", "", lrt_table_filtered$ENSEMBL_ID)

gene_symbols <- AnnotationDbi::select(org.Hs.eg.db, keys = clean_ids,
                       columns = c("SYMBOL", "GENENAME"), keytype = "ENSEMBL")

lrt_table_final <- merge(lrt_table_filtered, gene_symbols,
                         by.x = "ENSEMBL_ID", by.y = "ENSEMBL", all.x = TRUE)
lrt_table_final <- lrt_table_final[order(lrt_table_final$padj), ]

print(lrt_table_final[, c("ENSEMBL_ID", "SYMBOL", "log2FoldChange", "baseMean", "padj")])
write.csv(lrt_table_final, out("STGp_LRT_significant_genes_annotated.csv"),
          row.names = FALSE)

# Wald labels
wald_ids <- gsub("\\..*", "", rownames(res_sig_only))
wald_symbols <- AnnotationDbi::select(org.Hs.eg.db, keys = wald_ids,
                       columns = "SYMBOL", keytype = "ENSEMBL")

wald_label_df <- subset(res_df, padj < 0.05)
wald_label_df$ensembl_clean <- gsub("\\..*", "", rownames(wald_label_df))
wald_label_df$gene_label <- wald_symbols$SYMBOL[
  match(wald_label_df$ensembl_clean, wald_symbols$ENSEMBL)]

res_df$significant <- ifelse(res_df$padj < 0.05, "padj<0.05",
                             ifelse(res_df$pvalue < 0.05, "p<0.05", "ns"))

# 7.7 Volcano plots
p_wald <- ggplot(res_df,
                 aes(x = log2FoldChange, y = -log10(pvalue), color = significant)) +
  geom_point(alpha = 0.5, size = 1) +
  geom_text_repel(data = wald_label_df, aes(label = gene_label),
                  size = 3.5, box.padding = 0.5, color = "red") +
  scale_color_manual(values = c("padj<0.05" = "red",
                                "p<0.05"    = "orange",
                                "ns"        = "grey")) +
  theme_bw() +
  labs(title    = "Wald test",
       subtitle = paste0("Significant genes (padj<0.05): ", length(res_wald_sig)),
       x = "log2FoldChange", y = "-log10(p-value)")

# LRT labels (symbol; falls back to the Ensembl ID if no symbol is found)
lrt_label_df <- res_lrt_df[rownames(res_lrt_df) %in% rownames(lrt_sig_filtered), ]
lrt_label_df$ensembl_clean <- gsub("\\..*", "", rownames(lrt_label_df))

lrt_symbols <- AnnotationDbi::select(org.Hs.eg.db, keys = lrt_label_df$ensembl_clean,
                      columns = "SYMBOL", keytype = "ENSEMBL")

lrt_label_df$gene_label <- lrt_symbols$SYMBOL[
  match(lrt_label_df$ensembl_clean, lrt_symbols$ENSEMBL)]
lrt_label_df$gene_label <- ifelse(is.na(lrt_label_df$gene_label),
                                  lrt_label_df$ensembl_clean,
                                  lrt_label_df$gene_label)

res_lrt_df$significant_filtered <- ifelse(
  res_lrt_df$padj < 0.05 &
    abs(res_lrt_df$log2FoldChange) > 0.5 &
    res_lrt_df$baseMean < 100000, "padj<0.05",
  ifelse(res_lrt_df$pvalue < 0.05, "p<0.05", "ns"))

p_lrt <- ggplot(res_lrt_df,
                aes(x = log2FoldChange, y = -log10(pvalue),
                    color = significant_filtered)) +
  geom_point(alpha = 0.5, size = 1) +
  geom_text_repel(data = lrt_label_df, aes(label = gene_label),
                  size = 3, box.padding = 0.5, max.overlaps = 15, color = "red") +
  scale_color_manual(values = c("padj<0.05" = "red",
                                "p<0.05"    = "orange",
                                "ns"        = "grey")) +
  coord_cartesian(ylim = c(0, 10)) +
  theme_bw() +
  labs(title    = "LRT",
       subtitle = paste0("Significant genes after filtering: ", nrow(lrt_sig_filtered)),
       x = "log2FoldChange", y = "-log10(p-value)")

p_wald + p_lrt +
  plot_annotation(
    title = "Differential gene expression in STGp: Wald vs LRT",
    theme = theme(plot.title = element_text(face = "bold", size = 13))
  )
ggsave(out("STGp_volcano_Wald_vs_LRT_final.png"), width = 14, height = 5, dpi = 150)


# 8. Gene set enrichment analysis (GSEA, GO Biological Process) ---------------

# 8.1 Ranked gene list: sign(log2FC) * -log10(p), using LRT p-values
res_lrt_df$ensembl_clean <- gsub("\\..*", "", rownames(res_lrt_df))

res_lrt_clean <- res_lrt_df[!is.na(res_lrt_df$pvalue) &
                              !is.na(res_lrt_df$log2FoldChange), ]

# Replace p = 0 with the smallest non-zero p-value to keep ranks finite
min_pvalue <- min(res_lrt_clean$pvalue[res_lrt_clean$pvalue > 0], na.rm = TRUE)
res_lrt_clean$pvalue_fixed <- ifelse(res_lrt_clean$pvalue == 0,
                                     min_pvalue, res_lrt_clean$pvalue)

res_lrt_clean$rank_score <- sign(res_lrt_clean$log2FoldChange) *
  (-log10(res_lrt_clean$pvalue_fixed))

res_lrt_clean <- res_lrt_clean[!is.na(res_lrt_clean$rank_score) &
                                 !is.infinite(res_lrt_clean$rank_score), ]

ranked_lrt <- res_lrt_clean$rank_score
names(ranked_lrt) <- res_lrt_clean$ensembl_clean
ranked_lrt <- ranked_lrt[!duplicated(names(ranked_lrt))]
ranked_lrt <- sort(ranked_lrt, decreasing = TRUE)

cat("Genes in ranked list:", length(ranked_lrt), "\n")
cat("Any Inf:", any(is.infinite(ranked_lrt)), "\n")
cat("Any NA:",  any(is.na(ranked_lrt)), "\n")
cat("Top 5 upregulated in SZ:\n");   print(head(ranked_lrt, 5))
cat("Top 5 downregulated in SZ:\n"); print(tail(ranked_lrt, 5))

# 8.2 GSEA
gsea_lrt <- gseGO(
  geneList      = ranked_lrt,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENSEMBL",
  ont           = "BP",
  minGSSize     = 15,
  maxGSSize     = 500,
  pvalueCutoff  = 0.05,
  pAdjustMethod = "BH",
  verbose       = FALSE,
  seed          = 42
)

gsea_df <- as.data.frame(gsea_lrt)
cat("Significant GO terms (GSEA):", nrow(gsea_df), "\n")
cat("NES > 0 (activated in SZ):", sum(gsea_df$NES > 0), "\n")
cat("NES < 0 (suppressed in SZ):", sum(gsea_df$NES < 0), "\n")

gsea_df$direction <- ifelse(gsea_df$NES > 0, "Activated in SZ", "Suppressed in SZ")

top_pathways <- gsea_df %>%
  group_by(direction) %>%
  slice_min(p.adjust, n = 10) %>%
  ungroup() %>%
  arrange(NES) %>%
  mutate(Description = stringr::str_wrap(Description, width = 40),
         Description = factor(Description, levels = Description))

# 8.3 Visualization
ggplot(top_pathways,
       aes(x = NES, y = Description, size = setSize, color = p.adjust)) +
  geom_point() +
  scale_color_gradient(low = "red", high = "steelblue", labels = scales::scientific) +
  scale_size_continuous(range = c(4, 10)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
  facet_wrap(~direction, scales = "free_y", ncol = 1) +
  theme_bw(base_size = 12) +
  labs(title    = "GSEA: GO Biological Process in STGp (SZ vs HC)",
       subtitle = "Top 10 terms per direction, LRT-based ranking",
       x = "Normalized Enrichment Score (NES)", y = NULL,
       color = "p.adjust", size = "Gene set size") +
  theme(axis.text.y = element_text(size = 10, color = "black"),
        strip.text  = element_text(face = "bold", size = 12),
        plot.title  = element_text(face = "bold", size = 13))
ggsave(out("STGp_GSEA_LRT_top10.png"), width = 10, height = 9, dpi = 150)

write.csv(gsea_df, out("STGp_GSEA_LRT_results.csv"), row.names = FALSE)

cat("\nTop 5 activated in SZ:\n")
print(gsea_df %>% filter(NES > 0) %>% arrange(p.adjust) %>% head(5) %>%
        dplyr::select(Description, NES, p.adjust, setSize))

cat("\nTop 5 suppressed in SZ:\n")
print(gsea_df %>% filter(NES < 0) %>% arrange(p.adjust) %>% head(5) %>%
        dplyr::select(Description, NES, p.adjust, setSize))


# 9. Session information ------------------------------------------------------

writeLines(capture.output(sessionInfo()), out("sessionInfo.txt"))
