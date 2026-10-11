# =============================================================================
# GSE25624 – Neuroblastoma bone marrow: differential expression analysis
#
# Groups (from GEO sample titles):
#   I_BM  = metastatic NB, infiltrated bone marrow   (n = 9)
#   NI_BM = localized NB, non-infiltrated bone marrow (n = 23)
#   H_BM  = healthy young donors                      (n = 8)
#
# Steps: GEOquery -> log2 (if needed) + quantile normalisation -> limma
#        -> probe-to-gene collapsing -> DEG tables (FDR < 0.05, |log2FC| >= 1)
#        -> heatmaps, volcano plots, Venn diagram + intersections
#
# Run from the repository root:
#   Rscript scripts/deg_analysis.R
# Outputs are written to results/ and figures/.
# =============================================================================

suppressPackageStartupMessages({
  library(GEOquery)
  library(limma)
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
  library(VennDiagram)
  library(grid)
  library(pheatmap)
})

FDR_CUTOFF <- 0.05
LFC_CUTOFF <- 1

dir.create("data", showWarnings = FALSE)
dir.create("results/degs", recursive = TRUE, showWarnings = FALSE)
dir.create("results/venn", recursive = TRUE, showWarnings = FALSE)
dir.create("figures/volcano", recursive = TRUE, showWarnings = FALSE)
dir.create("figures/venn", recursive = TRUE, showWarnings = FALSE)
dir.create("figures/heatmaps", recursive = TRUE, showWarnings = FALSE)

# -----------------------------------------------------------------------------
# 1. Load GEO data (cached in data/ after first download)
# -----------------------------------------------------------------------------
rds_path <- "data/GSE25624_gset.rds"
if (file.exists(rds_path)) {
  gset <- readRDS(rds_path)
} else {
  gset <- getGEO("GSE25624", GSEMatrix = TRUE, AnnotGPL = TRUE)[[1]]
}
fvarLabels(gset) <- make.names(fvarLabels(gset))

# -----------------------------------------------------------------------------
# 2. Sample groups (suffix of the GEO sample title)
# -----------------------------------------------------------------------------
pheno <- pData(gset)
group <- factor(trimws(sub(".*-", "", pheno$title)))
stopifnot(all(levels(group) %in% c("H_BM", "NI_BM", "I_BM")))
gset$group <- group
print(table(group))

# -----------------------------------------------------------------------------
# 3. Normalisation: log2 if data look un-logged, then quantile normalisation
# -----------------------------------------------------------------------------
ex <- exprs(gset)
qx <- quantile(ex, c(0, .25, .5, .75, .99, 1), na.rm = TRUE)
if (qx[5] > 100 || (qx[6] - qx[1] > 50)) {
  ex[ex <= 0] <- NA
  ex <- log2(ex)
}
ex <- normalizeBetweenArrays(ex, method = "quantile")
exprs(gset) <- ex

# -----------------------------------------------------------------------------
# 4. limma: cell-means model + three contrasts
# -----------------------------------------------------------------------------
design <- model.matrix(~0 + group)
colnames(design) <- levels(group)
fit <- lmFit(gset, design)

cont.matrix <- makeContrasts(
  I_vs_H  = I_BM  - H_BM,
  NI_vs_H = NI_BM - H_BM,
  I_vs_NI = I_BM  - NI_BM,
  levels  = design
)
fit2 <- eBayes(contrasts.fit(fit, cont.matrix))

tt_all <- lapply(colnames(cont.matrix), function(ct)
  topTable(fit2, coef = ct, number = Inf, sort.by = "P"))
names(tt_all) <- colnames(cont.matrix)

# -----------------------------------------------------------------------------
# 5. Probe -> gene collapsing
#    FDR filter on probes, drop unannotated probes, keep the most significant
#    probe per gene, then apply the |log2FC| filter to that probe.
# -----------------------------------------------------------------------------
get_gene_DEGs <- function(tt, gene_col = "Gene.symbol",
                          fdr_cutoff = 0.05, lfc_cutoff = NULL) {
  tt$ProbeID <- rownames(tt)
  tt <- tt[tt$adj.P.Val < fdr_cutoff, ]
  tt <- tt[!is.na(tt[[gene_col]]) & tt[[gene_col]] != "", ]
  tt <- tt[order(tt$adj.P.Val), ]
  tt <- tt[!duplicated(tt[[gene_col]]), ]
  if (!is.null(lfc_cutoff)) tt <- tt[abs(tt$logFC) >= lfc_cutoff, ]
  tt
}

gene_DEGs_by_contrast <- lapply(tt_all, get_gene_DEGs,
                                fdr_cutoff = FDR_CUTOFF,
                                lfc_cutoff = LFC_CUTOFF)
stopifnot(all(sapply(gene_DEGs_by_contrast, function(tt)
  length(unique(tt$Gene.symbol)) == nrow(tt))))
print(sapply(gene_DEGs_by_contrast, nrow))

for (ct in names(gene_DEGs_by_contrast)) {
  tt <- gene_DEGs_by_contrast[[ct]]
  write.table(
    data.frame(ProbeID = tt$ProbeID, Gene.symbol = tt$Gene.symbol,
               log2FC = tt$logFC, AveExpr = tt$AveExpr, t = tt$t,
               P.Value = tt$P.Value, adj.P.Val = tt$adj.P.Val),
    file = file.path("results/degs", paste0("GENE_DEGs_", ct, "_FDR05.tsv")),
    sep = "\t", row.names = FALSE, quote = FALSE)
}

saveRDS(gset, rds_path)  # normalised object (not tracked by git)

# -----------------------------------------------------------------------------
# 6. Heatmaps: top 15 DEGs, z-scored, samples of the two groups compared
# -----------------------------------------------------------------------------
make_heatmap <- function(deg_table, title, gset, top_n = 50,
                         show_rownames = TRUE, filename = NA) {
  deg_table <- head(deg_table[order(deg_table$adj.P.Val), ], top_n)
  expr <- exprs(gset)[deg_table$ProbeID, ]
  rownames(expr) <- deg_table$Gene.symbol
  expr <- t(scale(t(expr)))
  expr <- expr[apply(expr, 1, function(x) all(is.finite(x))), ]
  ann_col <- data.frame(Group = gset$group)
  rownames(ann_col) <- colnames(expr)
  pheatmap(expr, annotation_col = ann_col,
           cluster_rows = TRUE, cluster_cols = TRUE,
           show_rownames = show_rownames, show_colnames = TRUE,
           fontsize_row = 8, fontsize_col = 10, border_color = NA,
           main = title, filename = filename, width = 10, height = 5)
}

gset_IH  <- gset[, gset$group %in% c("I_BM",  "H_BM")]; gset_IH$group  <- droplevels(gset_IH$group)
gset_NIH <- gset[, gset$group %in% c("NI_BM", "H_BM")]; gset_NIH$group <- droplevels(gset_NIH$group)

make_heatmap(gene_DEGs_by_contrast$I_vs_H,  "Metastatic vs Healthy (I vs H)",
             gset_IH,  top_n = 15, filename = "figures/heatmaps/heatmap_I_vs_H.png")
make_heatmap(gene_DEGs_by_contrast$NI_vs_H, "Localized vs Healthy (NI vs H)",
             gset_NIH, top_n = 15, filename = "figures/heatmaps/heatmap_NI_vs_H.png")

# -----------------------------------------------------------------------------
# 7. Volcano plots (gene-level DEGs, top 5 labelled per contrast)
#    NOTE: built from the filtered DEG tables, so only significant genes are
#    drawn. For a classic volcano, build volcano_df from collapsed tt_all.
# -----------------------------------------------------------------------------
contrast_map <- c(NI_vs_H = "Localized vs Healthy",
                  I_vs_H  = "Metastatic vs Healthy",
                  I_vs_NI = "Metastatic vs Localized")

volcano_df <- bind_rows(lapply(names(gene_DEGs_by_contrast), function(ct)
  gene_DEGs_by_contrast[[ct]] %>%
    mutate(Gene = Gene.symbol, Contrast = ct,
           negLogFDR = -log10(adj.P.Val),
           Significant = adj.P.Val < FDR_CUTOFF & abs(logFC) >= LFC_CUTOFF))) %>%
  mutate(Contrast_clean = factor(contrast_map[Contrast], levels = contrast_map))

top5_labels <- volcano_df %>% filter(Significant) %>%
  group_by(Contrast_clean) %>% arrange(adj.P.Val) %>%
  slice_head(n = 5) %>% ungroup()

p <- ggplot(volcano_df, aes(logFC, negLogFDR)) +
  geom_point(color = "grey75", alpha = 0.6, size = 1) +
  geom_point(data = subset(volcano_df, Significant), color = "red", size = 1) +
  geom_text_repel(data = top5_labels, aes(label = Gene), size = 3,
                  fontface = "bold", box.padding = 0.6, point.padding = 0.4,
                  force = 5, max.overlaps = 10, segment.size = 0.4) +
  facet_wrap(~ Contrast_clean, nrow = 1, drop = FALSE) +
  geom_vline(xintercept = c(-LFC_CUTOFF, LFC_CUTOFF), linetype = "dashed") +
  geom_hline(yintercept = -log10(FDR_CUTOFF), linetype = "dashed") +
  coord_cartesian(xlim = c(-6, 6), ylim = c(0, 14)) +
  labs(title = "Differential expression volcano plots (FDR < 0.05)",
       x = "log2 Fold Change", y = "-log10(FDR)") +
  theme_bw() +
  theme(strip.text = element_text(face = "bold", size = 11),
        plot.title = element_text(face = "bold"), legend.position = "none")

ggsave("figures/volcano/GSE25624_VOLCANO_GENELEVEL_FDR05.pdf", p, width = 14, height = 5)
ggsave("figures/volcano/GSE25624_VOLCANO_GENELEVEL_FDR05.png", p, width = 14, height = 5, dpi = 300)

# -----------------------------------------------------------------------------
# 8. Venn diagram + intersections
# -----------------------------------------------------------------------------
g <- lapply(gene_DEGs_by_contrast, `[[`, "Gene.symbol")
gene_sets <- list("Metastatic vs Healthy"   = g$I_vs_H,
                  "Localized vs Healthy"    = g$NI_vs_H,
                  "Metastatic vs Localized" = g$I_vs_NI)

futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")
venn.plot <- venn.diagram(x = gene_sets, fill = c("red", "blue", "green"),
                          alpha = 0.5, cex = 1.2, cat.cex = 1.0,
                          cat.dist = c(0.06, 0.06, 0.06), margin = 0.25,
                          main = "GSE25624 – Gene-level DEGs (FDR < 0.05)",
                          filename = NULL)

pdf("figures/venn/venn_genelevel_fdr05.pdf", width = 7, height = 7)
grid.draw(venn.plot); dev.off()
png("figures/venn/venn_genelevel_fdr05.png", width = 7, height = 7, units = "in", res = 300)
grid.draw(venn.plot); dev.off()

write_set <- function(x, f)
  write.table(data.frame(Gene = x), file.path("results/venn", f),
              sep = "\t", row.names = FALSE, quote = FALSE)

write_set(Reduce(intersect, g),                              "CORE_ALL.tsv")
write_set(intersect(g$I_vs_H, g$NI_vs_H),                    "CORE_DISEASE.tsv")
write_set(intersect(g$I_vs_H, g$I_vs_NI),                    "CORE_METASTASIS.tsv")
write_set(setdiff(g$I_vs_H,  union(g$NI_vs_H, g$I_vs_NI)),   "UNIQUE_I_vs_H.tsv")
write_set(setdiff(g$NI_vs_H, union(g$I_vs_H,  g$I_vs_NI)),   "UNIQUE_NI_vs_H.tsv")
write_set(setdiff(g$I_vs_NI, union(g$I_vs_H,  g$NI_vs_H)),   "UNIQUE_I_vs_NI.tsv")

print(sapply(gene_sets, length))
writeLines(capture.output(sessionInfo()), "results/sessionInfo.txt")
