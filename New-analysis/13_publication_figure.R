# ================================================================
# Publication-quality Figure: circAcc1 cross-species conservation
# ================================================================
# Follows: scientific-visualization skill best practices
# Target: Nature-family journal, 183mm (2-col) width, Arial font
# CVD-safe Okabe-Ito palette, PDF vector output
# ================================================================

library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

# ---- Global style (Nature-family journal specs) ----
theme_publication <- theme_minimal(base_size = 7, base_family = "sans") +
  theme(
    text = element_text(color = "#333333"),
    axis.title = element_text(size = 7, face = "bold"),
    axis.text = element_text(size = 6.5),
    strip.text = element_text(size = 7, face = "bold"),
    legend.text = element_text(size = 6.5),
    legend.title = element_text(size = 7, face = "bold"),
    legend.position = "bottom",
    legend.key.size = unit(3, "mm"),
    plot.title = element_text(size = 9, face = "bold"),
    plot.subtitle = element_text(size = 7, color = "#555555"),
    plot.margin = unit(c(2, 2, 2, 2), "mm"),
    panel.grid.major = element_line(linewidth = 0.2, color = "#E0E0E0"),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(fill = NA, color = "#CCCCCC", linewidth = 0.3),
    axis.ticks = element_line(linewidth = 0.2, color = "#CCCCCC"),
    axis.ticks.length = unit(1, "mm")
  )

# CVD-safe Okabe-Ito palette
okabe_ito <- c(
  "#E69F00", "#56B4E9", "#009E73", "#F0E442",
  "#0072B2", "#D55E00", "#CC79A7", "#000000"
)

# Species colors (consistent across all panels)
sp_colors <- c(
  "Human (hsa)"      = okabe_ito[1],
  "Rhesus (mma)"     = okabe_ito[2],
  "Mouse (mmu)"      = okabe_ito[3],
  "L. bat (mfu)"     = okabe_ito[5],
  "B. bat (mpi)"     = okabe_ito[6],
  "H. bat (rsi)"     = okabe_ito[7]
)

setwd("D:/跨物种筛选保守circRNA")
dir.create("New-analysis/figures", showWarnings = FALSE, recursive = TRUE)

# ---- Prepare data ----
df <- data.frame(
  species = factor(
    c("Human (hsa)", "Rhesus (mma)", "Mouse (mmu)",
      "L. bat (mfu)", "B. bat (mpi)", "H. bat (rsi)"),
    levels = c("Human (hsa)", "Rhesus (mma)", "Mouse (mmu)",
               "L. bat (mfu)", "B. bat (mpi)", "H. bat (rsi)")
  ),
  clade = c("Primate", "Primate", "Rodent", "Bat", "Bat", "Bat"),
  # Our data: milk exosome circACACA isoforms (CIRCexplorer3)
  milk_isoforms = c(0, 0, 9, 6, 1, 8),
  # circAtlas v3.0: all-tissue circACACA isoforms
  atlas_isoforms = c(317, 92, 56, NA, NA, NA),
  # ACACA gene information
  gene_chr = c("chr17", "chr16", "chr11", "scaffold0019", "scaffold0021", "NW_017739011.1"),
  # mRNA detection
  mrna = c("GTEx", "GTEx", "GTEx",
           "polyA-seq\n(mean 1,889)", "polyA-seq\n(mean 3,428)", "RefSeq"),
  stringsAsFactors = FALSE
)

# ---- Panel A: circAtlas isoform counts ----
cat("Creating Panel A: circAtlas public database...\n")

df_atlas <- df[!is.na(df$atlas_isoforms), ]
# Add annotation for total
df_atlas$label_y <- df_atlas$atlas_isoforms + 15
df_atlas$label_text <- paste0(df_atlas$atlas_isoforms, " isoforms")

pA <- ggplot(df_atlas, aes(x = species, y = atlas_isoforms, fill = species)) +
  geom_col(width = 0.55, alpha = 0.92) +
  geom_text(aes(label = label_text, y = label_y),
            size = 2.2, fontface = "bold", color = "#333333") +
  scale_fill_manual(values = sp_colors, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(
    title = "circACACA in public databases (circAtlas v3.0)",
    subtitle = "All tissues combined — human: 317, macaque: 92, mouse: 56 isoforms",
    x = "", y = "Number of circRNA isoforms"
  ) +
  theme_publication +
  theme(axis.text.x = element_text(angle = 25, hjust = 1, size = 6))

# ---- Panel B: Milk exosome circACACA (our data) ----
cat("Creating Panel B: Milk exosome circACACA...\n")

# Add detection indicator
df$detected <- ifelse(df$milk_isoforms > 0, "Detected", "Not detected")

pB <- ggplot(df, aes(x = species, y = milk_isoforms, fill = species)) +
  geom_col(width = 0.55, alpha = 0.92) +
  geom_text(
    aes(label = ifelse(milk_isoforms > 0,
                       paste(milk_isoforms, "isoforms"), "Not detected"),
        y = ifelse(milk_isoforms > 0, milk_isoforms + 0.45, 0.35)),
    size = 2.2, fontface = "bold",
    color = ifelse(df$milk_isoforms > 0, "#333333", "#CC0000")
  ) +
  scale_fill_manual(values = sp_colors, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.2)), limits = c(0, 10.5)) +
  labs(
    title = "circACACA in milk exosomes (this study)",
    subtitle = "CIRCexplorer3: detected in 4/6 species (bats + mouse); not in primates",
    x = "", y = "Number of circRNA isoforms"
  ) +
  theme_publication +
  theme(axis.text.x = element_text(angle = 25, hjust = 1, size = 6))

# ---- Panel C: Combined evidence summary ----
cat("Creating Panel C: Evidence summary matrix...\n")

evidence <- data.frame(
  species = factor(rep(c("Human\n(hsa)", "Rhesus\n(mma)", "Mouse\n(mmu)",
                          "L. bat\n(mfu)", "B. bat\n(mpi)", "H. bat\n(rsi)"), 4),
    levels = rev(c("Human\n(hsa)", "Rhesus\n(mma)", "Mouse\n(mmu)",
                    "L. bat\n(mfu)", "B. bat\n(mpi)", "H. bat\n(rsi)"))
  ),
  tier = factor(rep(c("(i) Gene in genome",
                       "(ii) mRNA expressed",
                       "(iii) circRNA in milk",
                       "(iv) circRNA in circAtlas"), each = 6),
    levels = c("(iv) circRNA in circAtlas",
               "(iii) circRNA in milk",
               "(ii) mRNA expressed",
               "(i) Gene in genome")
  ),
  status = c(
    # Gene
    "Yes", "Yes", "Yes", "Yes", "Yes", "Yes",
    # mRNA
    "Yes", "Yes", "Yes", "Yes", "Yes", "Yes",
    # circRNA milk
    "No", "No", "Yes (9)", "Yes (6)", "Yes (1)", "Yes (8)",
    # circAtlas
    "Yes (317)", "Yes (92)", "Yes (56)", "N/A", "N/A", "N/A"
  ),
  stringsAsFactors = FALSE
)

evidence$status_col <- ifelse(grepl("^Yes", evidence$status), "#27AE60",
                       ifelse(grepl("^No", evidence$status), "#E74C3C",
                       ifelse(evidence$status == "N/A", "#F39C12", "#BDC3C7")))

pC <- ggplot(evidence, aes(x = tier, y = species, fill = status_col)) +
  geom_tile(color = "white", linewidth = 0.8, width = 0.85, height = 0.85) +
  geom_text(aes(label = sub("Yes \\(|\\)", "", status)),
            size = 2.2, color = "white", fontface = "bold") +
  scale_fill_identity() +
  labs(
    title = "Cross-species evidence summary",
    subtitle = "Green = confirmed; Red = not detected; Orange = no data; Gray = not in database",
    x = "", y = ""
  ) +
  theme_minimal(base_size = 7) +
  theme(
    axis.text.x = element_text(size = 6, face = "bold"),
    axis.text.y = element_text(size = 6.5),
    panel.grid = element_blank(),
    plot.title = element_text(size = 9, face = "bold"),
    plot.subtitle = element_text(size = 6.5, color = "#555555"),
    plot.margin = unit(c(2, 2, 2, 2), "mm"),
    legend.position = "none"
  )

# ---- Panel D: Schematic illustration ----
cat("Creating Panel D: Key finding schematic...\n")

# Using ggplot to create a text-based schematic
schematic_text <- data.frame(
  x = 0.5, y = c(0.85, 0.65, 0.45, 0.25, 0.08),
  label = c(
    "ACACA gene: present in all 6 species",
    "mRNA: expressed in all 6 species",
    "circRNA in milk exosomes:\n  4/6 detected (all bats + mouse)",
    "circAtlas public DB:\n  317 (human) / 92 (macaque) / 56 (mouse) isoforms",
    "Key insight: circACACA is conserved at gene level;\nabsence from primate milk may reflect limited\nsampling rather than biological absence"
  ),
  size = c(3.2, 3.2, 3.2, 3.2, 2.8),
  color = c("#27AE60", "#27AE60", "#3498DB", "#F39C12", "#555555"),
  fontface = c("bold", "bold", "bold", "bold", "plain"),
  stringsAsFactors = FALSE
)

pD <- ggplot(schematic_text, aes(x = x, y = y)) +
  geom_text(aes(label = label, size = size, color = color, fontface = fontface),
            hjust = 0.5, vjust = 0.5) +
  scale_size_identity() +
  scale_color_identity() +
  # Connecting arrows
  annotate("segment", x = 0.5, xend = 0.5, y = 0.78, yend = 0.70,
           arrow = arrow(length = unit(1.5, "mm"), type = "closed"), linewidth = 0.5, color = "#999999") +
  annotate("segment", x = 0.5, xend = 0.5, y = 0.58, yend = 0.50,
           arrow = arrow(length = unit(1.5, "mm"), type = "closed"), linewidth = 0.5, color = "#999999") +
  annotate("segment", x = 0.5, xend = 0.5, y = 0.38, yend = 0.30,
           arrow = arrow(length = unit(1.5, "mm"), type = "closed"), linewidth = 0.5, color = "#999999") +
  scale_x_continuous(limits = c(0, 1)) +
  scale_y_continuous(limits = c(0, 0.95)) +
  labs(title = "Summary & Interpretation") +
  theme_void(base_size = 7) +
  theme(
    plot.title = element_text(size = 9, face = "bold", hjust = 0.5,
                              margin = margin(b = 3, unit = "mm"))
  )

# ---- Assemble final figure ----
cat("Assembling multi-panel figure...\n")

# Layout: A (top-left) | B (top-right)
#         C (bottom-left) | D (bottom-right)
final_fig <- (pA | pB) / (pC | pD) +
  plot_layout(heights = c(1, 1.15)) +
  plot_annotation(
    tag_levels = "A",
    theme = theme(
      plot.tag = element_text(size = 10, face = "bold", family = "sans"),
      plot.title = element_text(size = 11, face = "bold", family = "sans",
                               margin = margin(b = 4, unit = "mm")),
      plot.subtitle = element_text(size = 8, family = "sans", color = "#555555")
    )
  )

# ---- Export at publication resolution ----
# Nature 2-column: 183mm × suitable height
# Convert to inches: 183mm = 7.2 inches
cat("Exporting PDF...\n")
ggsave("New-analysis/figures/Fig_circAcc1_publication.pdf",
       final_fig,
       width = 183, height = 200, units = "mm",
       dpi = 600,
       device = cairo_pdf)

cat("Exporting TIFF (600 dpi)...\n")
ggsave("New-analysis/figures/Fig_circAcc1_publication.tiff",
       final_fig,
       width = 183, height = 200, units = "mm",
       dpi = 600,
       compression = "lzw")

cat("Exporting PNG preview...\n")
ggsave("New-analysis/figures/Fig_circAcc1_publication.png",
       final_fig,
       width = 183, height = 200, units = "mm",
       dpi = 300)

cat("\n========== Publication figure complete ==========\n")
cat("Output:\n")
cat("  New-analysis/figures/Fig_circAcc1_publication.pdf  (vector, 600 dpi)\n")
cat("  New-analysis/figures/Fig_circAcc1_publication.tiff (600 dpi, LZW)\n")
cat("  New-analysis/figures/Fig_circAcc1_publication.png  (preview, 300 dpi)\n")
cat("\nFigure spec: 183 × 200 mm (Nature 2-column width)\n")
