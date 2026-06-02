# ============================================================
# Single-panel dot plot: circAcc1 cross-species evidence
# Publication quality — Nature-family journal specs
# ============================================================

library(ggplot2)
library(dplyr)
library(tidyr)

setwd("D:/跨物种筛选保守circRNA")
dir.create("New-analysis/figures", showWarnings = FALSE, recursive = TRUE)

# ---- Build data ----
df <- data.frame(
  species = factor(
    c("Human (hsa)", "Rhesus (mma)", "Mouse (mmu)",
      "L-winged bat (mfu)", "B-winged bat (mpi)", "Horseshoe bat (rsi)"),
    levels = rev(c("Human (hsa)", "Rhesus (mma)", "Mouse (mmu)",
                   "L-winged bat (mfu)", "B-winged bat (mpi)", "Horseshoe bat (rsi)"))
  ),
  clade = c("Primate", "Primate", "Rodent", "Bat", "Bat", "Bat"),
  atlas_n = c(317, 92, 56, NA, NA, NA),
  milk_n = c(0, 0, 9, 6, 1, 8),
  gene = c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE),
  mrna = c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE),
  mrna_label = c("GTEx", "GTEx", "GTEx", "polyA 1,889", "polyA 3,428", "RefSeq"),
  stringsAsFactors = FALSE
)

# Long format
df_long <- df %>%
  pivot_longer(cols = c(atlas_n, milk_n), names_to = "source", values_to = "isoforms") %>%
  mutate(
    source_label = ifelse(source == "atlas_n",
                          "circAtlas v3.0 (all tissues)", "Milk exosomes (this study)"),
    source_label = factor(source_label,
      levels = c("circAtlas v3.0 (all tissues)", "Milk exosomes (this study)")),
    isoform_label = ifelse(is.na(isoforms), "N/A",
                    ifelse(isoforms == 0, "0", as.character(isoforms))),
    dot_status = ifelse(is.na(isoforms), "No data",
                 ifelse(isoforms > 0, "Detected", "Not detected")),
    dot_status = factor(dot_status,
      levels = c("Detected", "Not detected", "No data")),
    dot_size = ifelse(is.na(isoforms), 6, ifelse(isoforms == 0, 5, sqrt(isoforms) * 0.9))
  )

# ---- Theme ----
theme_pub <- theme_minimal(base_size = 8, base_family = "sans") +
  theme(
    text = element_text(color = "#2C2C2C"),
    plot.title = element_text(size = 10, face = "bold", margin = margin(b = 1)),
    plot.subtitle = element_text(size = 7.5, color = "#555555", margin = margin(b = 3)),
    plot.caption = element_text(size = 6.5, color = "#888888", margin = margin(t = 3)),
    axis.title = element_text(size = 8, face = "bold"),
    axis.text = element_text(size = 7.5),
    axis.text.y = element_text(size = 8, face = "bold"),
    panel.grid.major.y = element_line(linewidth = 0.3, color = "#E8E8E8"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin = unit(c(3, 6, 3, 2), "mm"),
    legend.position = "right",
    legend.title = element_text(size = 7, face = "bold"),
    legend.text = element_text(size = 7),
    legend.key.size = unit(4, "mm"),
    legend.spacing.y = unit(1.5, "mm")
  )

# ---- Plot ----
cat("Creating single-panel dot plot...\n")

p <- ggplot(df_long, aes(x = source_label, y = species)) +

  # Background shading for clades
  annotate("rect",
    xmin = -Inf, xmax = Inf,
    ymin = c(5.55, 4.55, 3.55, 0.55, 1.55, 2.55),
    ymax = c(6.45, 5.45, 4.45, 1.45, 2.45, 3.45),
    fill = c("#FDF6EC", "#FDF6EC", "#FDF6EC", "#ECF3F9", "#ECF3F9", "#ECF3F9"),
    alpha = 0.6
  ) +

  # Clade divider
  geom_hline(yintercept = 3.5, linewidth = 0.4, color = "#AAAAAA", linetype = "dotted") +

  # Dots
  geom_point(aes(size = dot_size, fill = dot_status, color = dot_status),
             shape = 21, stroke = 0.4, alpha = 0.95) +

  # Labels inside dots — FIXED size to avoid NA fontsize error
  geom_text(aes(label = isoform_label), size = 2.6, color = "white",
            fontface = "bold", lineheight = 0.85) +

  # Gene checkmark (left side)
  geom_text(data = df, aes(x = 0.38, y = species, label = ""),
            inherit.aes = FALSE, size = 0) +
  annotate("text", x = 0.38,
           y = seq_along(levels(df$species)),
           label = rep("", 6), size = 0) +

  # mRNA annotation (right side)
  annotate("text", x = 2.9,
           y = seq_along(levels(df$species)),
           label = paste0("mRNA: ", df$mrna_label[6:1]),
           size = 2.0, hjust = 0, color = "#777777", fontface = "italic") +

  # Clade labels
  annotate("text", x = 3.5, y = 6, label = "Primates", size = 2.6,
           fontface = "bold", color = "#B8860B") +
  annotate("text", x = 3.5, y = 4, label = "Rodent", size = 2.6,
           fontface = "bold", color = "#CD853F") +
  annotate("text", x = 3.5, y = 2, label = "Bats", size = 2.6,
           fontface = "bold", color = "#4682B4") +

  # Clade brackets
  annotate("segment", x = 3.35, xend = 3.35, y = 5.55, yend = 6.45,
           linewidth = 0.5, color = "#B8860B") +
  annotate("segment", x = 3.35, xend = 3.35, y = 3.55, yend = 4.45,
           linewidth = 0.5, color = "#CD853F") +
  annotate("segment", x = 3.35, xend = 3.35, y = 0.55, yend = 3.45,
           linewidth = 0.5, color = "#4682B4") +

  # Scales
  scale_size_continuous(
    range = c(3.5, 16),
    breaks = c(6, 10, 14, 18),
    labels = c("0 / N/A", "10", "50", "100+"),
    name = "Isoform\ncount"
  ) +
  scale_fill_manual(
    values = c("Detected" = "#27AE60", "Not detected" = "#E74C3C", "No data" = "#BDC3C7"),
    name = "Status",
    drop = FALSE
  ) +
  scale_color_manual(
    values = c("Detected" = "#1B7A43", "Not detected" = "#B03A2E", "No data" = "#95A5A6"),
    name = "Status",
    drop = FALSE
  ) +
  scale_x_discrete(expand = expansion(add = c(0.8, 1.8))) +

  # Labels
  labs(
    title = "circAcc1 (ACACA) — Cross-species Conservation",
    subtitle = paste0(
      "Dot area ∝ isoform count. ",
      "Gene present in all 6 genomes. ",
      "mRNA expressed in all 6 species."
    ),
    caption = paste0(
      "circAtlas v3.0 (ngdc.cncb.ac.cn/circatlas/) | ",
      "Milk exosome data: CIRCexplorer3, this study | ",
      "Human: 317 isoforms across tissues (circAtlas), 0 in milk"
    ),
    x = "", y = ""
  ) +
  guides(
    size = guide_legend(order = 1, override.aes = list(fill = "#555555")),
    fill = guide_legend(order = 2, override.aes = list(size = 5)),
    color = "none"
  ) +
  theme_pub

# ---- Export ----
cat("Exporting...\n")
ggsave("New-analysis/figures/Fig_circAcc1_dotplot.pdf", p,
       width = 183, height = 120, units = "mm", device = cairo_pdf, dpi = 600)
ggsave("New-analysis/figures/Fig_circAcc1_dotplot.png", p,
       width = 183, height = 120, units = "mm", dpi = 300)
cat("Done!\n")
