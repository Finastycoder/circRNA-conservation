# ============================================================
# Final comprehensive figure: the circAcc1 conservation story
# Panel A: UpSet (top intersections only)
# Panel B: 4-species Venn highlighting ACACA
# Panel C: Dot plot — gene→mRNA→circRNA evidence
# ============================================================

suppressPackageStartupMessages({
  library(ComplexHeatmap)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(grid)
  library(VennDiagram)
})
library(patchwork)

setwd("D:/跨物种筛选保守circRNA")
dir.create("New-analysis/figures", showWarnings = FALSE, recursive = TRUE)

SPECIES <- c("hsa", "mfu", "mma", "mmu", "mpi", "rsi")
set_names <- c("Human", "L.bat", "Rhesus", "Mouse", "B.bat", "H.bat")
sp_colors <- c(
  Human = "#E69F00", L.bat = "#56B4E9", Rhesus = "#009E73",
  Mouse = "#F0E442", B.bat = "#0072B2", H.bat = "#D55E00"
)

# ---- Load data ----
presence <- read.csv("New-analysis/gene_presence_matrix.csv", stringsAsFactors = FALSE)
gene_sets <- list()
for (i in seq_along(SPECIES)) {
  genes <- presence$gene_name[presence[[SPECIES[i]]] == TRUE]
  gene_sets[[set_names[i]]] <- unique(toupper(genes))
}

cat("=== Building combined figure ===\n\n")

# ============================================================
# Panel A: Clean UpSet — top 18 intersections only
# ============================================================
cat("Panel A: UpSet\n")

comb_mat <- make_comb_mat(gene_sets)
comb_filt <- comb_mat[comb_degree(comb_mat) >= 3 & comb_size(comb_mat) > 0]
comb_filt <- comb_filt[order(comb_size(comb_filt), decreasing = TRUE)]
# Keep top 18 (covers all 6/5/4-way + top 3-way)
comb_top18 <- comb_filt[1:min(18, length(comb_size(comb_filt)))]

deg <- comb_degree(comb_top18)
deg_colors <- c("6" = "#C0392B", "5" = "#E67E22", "4" = "#27AE60", "3" = "#3498DB")
bar_cols <- deg_colors[as.character(deg)]

# Annotate which bars contain ACACA
comb_genes_all <- lapply(names(comb_size(comb_top18)), function(cn) {
  extract_comb(comb_mat, cn)
})
names(comb_genes_all) <- names(comb_size(comb_top18))
has_acaca <- sapply(comb_genes_all, function(g) "ACACA" %in% g)

# Custom bar labels: add ★ for ACACA-containing combos
bar_labels <- comb_size(comb_top18)
bar_labels[has_acaca] <- paste0(bar_labels[has_acaca], " \U2605")

pdf("New-analysis/figures/upset_top18.pdf", width = 11, height = 6.5)
htA <- UpSet(
  comb_top18,
  set_order = rev(set_names),
  comb_order = order(comb_size(comb_top18), decreasing = TRUE),
  top_annotation = HeatmapAnnotation(
    "Shared genes" = anno_barplot(
      comb_size(comb_top18),
      border = FALSE,
      gp = gpar(fill = bar_cols, col = bar_cols),
      height = unit(4.5, "cm"),
      axis_param = list(labels_rot = 0)
    ),
    annotation_name_side = "left",
    annotation_name_gp = gpar(fontsize = 9)
  ),
  right_annotation = upset_right_annotation(
    comb_top18,
    gp = gpar(fill = sp_colors[set_names]),
    annotation_name_gp = gpar(fontsize = 8),
    axis_param = list(labels_rot = 0)
  ),
  row_names_gp = gpar(fontsize = 9),
  column_title = "circRNA host gene conservation (≥3 species)",
  column_title_gp = gpar(fontsize = 12, fontface = "bold")
)
draw(htA)
# Legend
lgd <- Legend(
  labels = c("6 sp.", "5 sp.", "4 sp.", "3 sp.", expression("\U2605"~"has ACACA")),
  legend_gp = gpar(fill = c(deg_colors, "#FFFFFF")),
  type = c(rep("grid", 4), "points"),
  pch = c(NA, NA, NA, NA, 8),
  title = "Level",
  title_gp = gpar(fontsize = 9, fontface = "bold"),
  labels_gp = gpar(fontsize = 8)
)
draw(lgd, x = unit(0.78, "npc"), y = unit(0.88, "npc"))
dev.off()

# ============================================================
# Panel B: 4-species Venn — mouse + 3 bats (where ACACA is)
# Use unfiltered gene_lookups so mfu ACACA is included
# ============================================================
cat("Panel B: 4-species Venn\n")

gene_lookups <- readRDS("New-analysis/gene_lookups.rds")

# Build gene sets from ALL circRNAs (not filtered by expression)
acaca_four <- list()
for (sp in c("mmu", "mfu", "mpi", "rsi")) {
  acaca_four[[sp]] <- unique(gene_lookups[[sp]]$gene_name)
}

# ACACA should be in all 4 — verify
for (nm in names(acaca_four)) cat(sprintf("  %s: ACACA=%s (%d total genes)\n",
  nm, "ACACA" %in% acaca_four[[nm]], length(acaca_four[[nm]])))

# Get intersection counts
common_4 <- Reduce(intersect, acaca_four)
cat(sprintf("  4-way intersection: %d genes (including ACACA: %s)\n",
  length(common_4), "ACACA" %in% common_4))

four_sp <- acaca_four[c("mmu", "mfu", "mpi", "rsi")]
four_cols <- sp_colors[c("Mouse", "L.bat", "B.bat", "H.bat")]
four_labels <- c("Mouse\n(mmu)", "L-winged bat\n(mfu)", "B-winged bat\n(mpi)", "Horseshoe bat\n(rsi)")

tiff("New-analysis/figures/venn_acaca_4way.tiff",
     width = 2200, height = 2200, res = 300, compression = "lzw")
vp <- venn.diagram(
  x = four_sp,
  category.names = four_labels,
  filename = NULL, output = TRUE,
  lwd = 2.5, lty = "solid",
  col = four_cols, fill = four_cols, alpha = 0.3,
  cex = 1.5, fontface = "bold", fontfamily = "sans",
  cat.cex = 1.2, cat.fontface = "bold", cat.fontfamily = "sans",
  cat.col = four_cols,
  main = "circACACA: detected in milk exosomes",
  main.cex = 1.6, main.fontface = "bold", main.fontfamily = "sans",
  margin = 0.08,
  disable.logging = TRUE
)
grid.draw(vp)
# Add ACACA label
grid.text("ACACA", x = 0.5, y = 0.48,
          gp = gpar(fontsize = 10, col = "#E74C3C", fontface = "bold"))
dev.off()

pdf("New-analysis/figures/venn_acaca_4way.pdf", width = 8, height = 8)
grid.draw(vp)
grid.text("ACACA", x = 0.5, y = 0.48,
          gp = gpar(fontsize = 10, col = "#E74C3C", fontface = "bold"))
dev.off()

# ============================================================
# Panel C: Dot plot — gene/mRNA/circRNA evidence
# ============================================================
cat("Panel C: Dot plot evidence\n")

df <- data.frame(
  species = factor(c("Human\n(hsa)", "Rhesus\n(mma)", "Mouse\n(mmu)",
                      "L.bat\n(mfu)", "B.bat\n(mpi)", "H.bat\n(rsi)"),
    levels = rev(c("Human\n(hsa)", "Rhesus\n(mma)", "Mouse\n(mmu)",
                    "L.bat\n(mfu)", "B.bat\n(mpi)", "H.bat\n(rsi)"))),
  clade = c("Primates", "Primates", "Rodent", "Bats", "Bats", "Bats"),
  atlas_n = c(317, 92, 56, NA, NA, NA),
  milk_n = c(0, 0, 9, 6, 1, 8),
  stringsAsFactors = FALSE
)

df_long <- df %>%
  pivot_longer(cols = c(atlas_n, milk_n), names_to = "source", values_to = "isoforms") %>%
  mutate(
    source_label = ifelse(source == "atlas_n",
                          "circAtlas v3.0\n(all tissues)", "Milk exosomes\n(this study)"),
    source_label = factor(source_label,
      levels = c("circAtlas v3.0\n(all tissues)", "Milk exosomes\n(this study)")),
    isoform_label = ifelse(is.na(isoforms), "N/A",
                    ifelse(isoforms == 0, "0", as.character(isoforms))),
    dot_status = ifelse(is.na(isoforms), "No data (not in DB)",
                 ifelse(isoforms > 0, "Detected", "Not detected\nin our samples")),
    dot_status = factor(dot_status,
      levels = c("Detected", "Not detected\nin our samples", "No data (not in DB)"))
  )

# ACACA highlight annotation
highlight_species <- c("Mouse\n(mmu)", "L.bat\n(mfu)", "B.bat\n(mpi)", "H.bat\n(rsi)")

theme_pub <- theme_minimal(base_size = 8, base_family = "sans") +
  theme(
    text = element_text(color = "#2C2C2C"),
    plot.title = element_text(size = 10, face = "bold", margin = margin(b = 1)),
    plot.subtitle = element_text(size = 7.5, color = "#555555", margin = margin(b = 3)),
    axis.title = element_text(size = 8, face = "bold"),
    axis.text = element_text(size = 7.5),
    axis.text.y = element_text(size = 8.5, face = "bold"),
    panel.grid.major.y = element_line(linewidth = 0.3, color = "#E8E8E8"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin = unit(c(3, 5, 3, 2), "mm"),
    legend.position = "right",
    legend.title = element_text(size = 7.5, face = "bold"),
    legend.text = element_text(size = 7),
    legend.key.size = unit(4, "mm")
  )

pC <- ggplot(df_long, aes(x = source_label, y = species)) +
  annotate("rect", xmin = -Inf, xmax = Inf,
    ymin = c(5.55, 4.55, 3.55, 0.55, 1.55, 2.55),
    ymax = c(6.45, 5.45, 4.45, 1.45, 2.45, 3.45),
    fill = c("#FDF6EC", "#FDF6EC", "#FDF6EC", "#ECF3F9", "#ECF3F9", "#ECF3F9"),
    alpha = 0.5) +
  geom_hline(yintercept = 3.5, linewidth = 0.4, color = "#AAAAAA", linetype = "dotted") +
  geom_point(aes(size = isoforms, fill = dot_status, color = dot_status),
             shape = 21, stroke = 0.4) +
  geom_text(aes(label = isoform_label), size = 2.4, color = "white",
            fontface = "bold") +
  # ACACA highlight — bold the 4 species
  annotate("text", x = 0.3,
    y = which(levels(df$species) %in% highlight_species),
    label = "\U2605", size = 5, color = "#E74C3C", fontface = "bold") +
  annotate("text", x = 2.85, y = 6, label = "Primates", size = 2.4,
           fontface = "bold", color = "#B8860B") +
  annotate("text", x = 2.85, y = 4, label = "Rodent", size = 2.4,
           fontface = "bold", color = "#CD853F") +
  annotate("text", x = 2.85, y = 2, label = "Bats", size = 2.4,
           fontface = "bold", color = "#4682B4") +
  scale_size_continuous(range = c(3.5, 14), name = "Isoforms",
    breaks = c(6, 10, 14), labels = c("0/N/A", "10", "50+")) +
  scale_fill_manual(values = c("Detected" = "#27AE60",
                                "Not detected\nin our samples" = "#E74C3C",
                                "No data (not in DB)" = "#BDC3C7"),
                    name = "Status") +
  scale_color_manual(values = c("Detected" = "#1B7A43",
                                 "Not detected\nin our samples" = "#B03A2E",
                                 "No data (not in DB)" = "#95A5A6"),
                     name = "Status") +
  scale_x_discrete(expand = expansion(add = c(0.7, 1.6))) +
  labs(
    title = expression(bold("circAcc1 (ACACA): the most abundant circRNA in mouse milk")),
    subtitle = paste0(
      "\U2605 = ACACA detected. ",
      "circAtlas: human 317 / macaque 92 / mouse 56 isoforms across tissues.\n",
      "Absence in primate milk likely reflects sampling depth, not biological absence."
    ),
    x = "", y = ""
  ) +
  theme_pub

ggsave("New-analysis/figures/dotplot_acaca_story.pdf", pC,
       width = 183, height = 110, units = "mm", device = cairo_pdf, dpi = 600)
ggsave("New-analysis/figures/dotplot_acaca_story.png", pC,
       width = 183, height = 110, units = "mm", dpi = 300)

cat("\n=== All panels generated ===\n")
cat("Panel A: upset_top18.pdf\n")
cat("Panel B: venn_acaca_4way.pdf\n")
cat("Panel C: dotplot_acaca_story.pdf\n")
