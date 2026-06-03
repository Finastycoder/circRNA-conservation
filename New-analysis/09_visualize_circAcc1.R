# ============================================================
# circAcc1 (ACACA) cross-species conservation visualization
# ============================================================

library(dplyr)
library(tidyr)
library(ggplot2)

source(if (file.exists("New-analysis/00_project_setup.R")) "New-analysis/00_project_setup.R" else "00_project_setup.R", encoding = "UTF-8")
dir.create("New-analysis/figures", showWarnings = FALSE, recursive = TRUE)

# ---- Data: circAcc1 detection per species ----
# Compiled from CIRIquant + CIRCexplorer3 analysis

acaca_data <- data.frame(
  species = c("hsa", "mma", "mfu", "mpi", "rsi", "mmu"),
  species_label = c("Human\n(hsa)", "Rhesus\n(mma)", "Long-winged bat\n(mfu)",
                    "Big-footed bat\n(mpi)", "Horseshoe bat\n(rsi)", "Mouse\n(mmu)"),
  category = c("Primate", "Primate", "Bat", "Bat", "Bat", "Rodent"),
  # CIRIquant
  cq_detected   = c(FALSE, FALSE, FALSE, FALSE, TRUE, TRUE),
  cq_isoforms   = c(0, 0, 0, 0, 8, 9),
  cq_max_bsj    = c(0, 0, 0, 0, 4, 442),
  # CIRCexplorer3
  ce_detected   = c(FALSE, FALSE, TRUE, TRUE, TRUE, TRUE),
  ce_isoforms   = c(0, 0, 6, 1, 8, 9),
  ce_max_reads  = c(0, 0, 6, 27, 9, 1683),
  ce_mean_reads = c(0, 0, 3.3, 17.3, 4.5, 307),
  stringsAsFactors = FALSE
)

# Add combined detection status
acaca_data$status <- ifelse(acaca_data$ce_detected & acaca_data$cq_detected, "Both tools",
                     ifelse(acaca_data$ce_detected, "CIRCexplorer3 only",
                     ifelse(acaca_data$cq_detected, "CIRIquant only",
                            "Not detected")))

# ---- Figure 1: Combined detection dot plot ----
cat("Figure 1: Detection overview\n")

acaca_long <- acaca_data %>%
  select(species_label, category, status,
         `CIRCexplorer3` = ce_isoforms,
         `CIRIquant` = cq_isoforms) %>%
  pivot_longer(cols = c("CIRCexplorer3", "CIRIquant"),
               names_to = "tool", values_to = "isoforms")

acaca_long$tool <- factor(acaca_long$tool, levels = c("CIRIquant", "CIRCexplorer3"))
acaca_long$species_label <- factor(acaca_long$species_label,
  levels = c("Human\n(hsa)", "Rhesus\n(mma)", "Long-winged bat\n(mfu)",
             "Big-footed bat\n(mpi)", "Horseshoe bat\n(rsi)", "Mouse\n(mmu)"))

p1 <- ggplot(acaca_long, aes(x = tool, y = species_label)) +
  geom_point(aes(size = isoforms, color = isoforms > 0), alpha = 0.85) +
  geom_text(aes(label = ifelse(isoforms > 0, isoforms, "")), size = 3.5,
            color = "white", fontface = "bold") +
  scale_size_continuous(range = c(4, 14), guide = "none") +
  scale_color_manual(values = c("FALSE" = "gray80", "TRUE" = "#E74C3C"),
                     guide = "none") +
  facet_wrap(~ category, scales = "free_y", ncol = 1, strip.position = "right") +
  labs(title = "circAcc1 (ACACA) Detection Across Species",
       subtitle = "Circle size = number of circRNA isoforms; number inside = isoform count",
       x = "Detection Tool", y = "") +
  theme_minimal(base_size = 14) +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text = element_text(size = 11, color = "gray40"),
        plot.title = element_text(face = "bold"))

ggsave("New-analysis/figures/circAcc1_detection_dotplot.pdf", p1,
       width = 7, height = 5.5)
ggsave("New-analysis/figures/circAcc1_detection_dotplot.png", p1,
       width = 7, height = 5.5, dpi = 150)
cat("  Saved: figures/circAcc1_detection_dotplot.pdf/png\n")

# ---- Figure 2: Expression level comparison ----
cat("Figure 2: Expression comparison\n")

expr_data <- acaca_data %>%
  filter(ce_detected) %>%
  select(species_label, ce_isoforms, ce_max_reads, ce_mean_reads) %>%
  pivot_longer(cols = c(ce_max_reads, ce_mean_reads),
               names_to = "metric", values_to = "reads")

expr_data$metric <- factor(expr_data$metric,
  levels = c("ce_max_reads", "ce_mean_reads"),
  labels = c("Max reads (most expressed isoform)", "Mean reads (across isoforms)"))

p2 <- ggplot(expr_data, aes(x = species_label, y = reads, fill = species_label)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  geom_text(aes(label = round(reads, 1)), vjust = -0.5, size = 3.5) +
  facet_wrap(~ metric, scales = "free_y", ncol = 1) +
  scale_fill_manual(values = c("#3498DB", "#2ECC71", "#E74C3C", "#F39C12"),
                    guide = "none") +
  labs(title = "circAcc1 Expression Level (CIRCexplorer3 readNumber)",
       subtitle = paste0("Mouse (mmu) shows strongest expression by far"),
       x = "", y = "Read count") +
  theme_minimal(base_size = 13) +
  theme(axis.text.x = element_text(size = 10),
        panel.grid.major.x = element_blank())

ggsave("New-analysis/figures/circAcc1_expression.pdf", p2,
       width = 8, height = 6)
ggsave("New-analysis/figures/circAcc1_expression.png", p2,
       width = 8, height = 6, dpi = 150)
cat("  Saved: figures/circAcc1_expression.pdf/png\n")

# ---- Figure 3: Isoform count + detection status ----
cat("Figure 3: Isoform landscape\n")

p3 <- ggplot(acaca_data, aes(x = species_label, y = ce_isoforms)) +
  geom_col(aes(fill = status), width = 0.65, alpha = 0.9) +
  geom_text(aes(label = paste0(ce_isoforms, " isoforms"),
                y = ce_isoforms + 0.3), size = 4, fontface = "bold") +
  scale_fill_manual(values = c("Both tools" = "#27AE60",
                                "CIRCexplorer3 only" = "#F39C12",
                                "Not detected" = "#BDC3C7"),
                    name = "Detection") +
  labs(title = "circAcc1 Isoforms Detected per Species",
       subtitle = "CIRCexplorer3 data | mmu (mouse) = reference species",
       x = "", y = "Number of circRNA isoforms") +
  theme_minimal(base_size = 13) +
  theme(axis.text.x = element_text(size = 10),
        legend.position = "bottom",
        panel.grid.major.x = element_blank())

ggsave("New-analysis/figures/circAcc1_isoforms.pdf", p3,
       width = 8, height = 5)
ggsave("New-analysis/figures/circAcc1_isoforms.png", p3,
       width = 8, height = 5, dpi = 150)
cat("  Saved: figures/circAcc1_isoforms.pdf/png\n")

# ---- Figure 4: Conservation summary panel ----
cat("Figure 4: Combined panel\n\n")

# Build a presence matrix for ACACA
sp_levels <- c("hsa", "mma", "mfu", "mpi", "rsi", "mmu")
sp_labels <- c("Human", "Rhesus", "Long-winged\nbat", "Big-footed\nbat",
               "Horseshoe\nbat", "Mouse")

conservation_text <- paste0(
  "circAcc1 (Acetyl-CoA Carboxylase 1)\n",
  "Cross-species Conservation Summary\n\n",
  "━━━━━━━━━━━━━━━━━━━━━━━━━\n\n",
  "Detected in ", sum(acaca_data$ce_detected), " / 6 species\n\n",
  "  ✅ Mouse    (mmu) — 9 isoforms, strongest expression\n",
  "  ✅ R.bat    (rsi) — 8 isoforms, consistent across 4 reps\n",
  "  ✅ B.bat    (mpi) — 1 isoform, moderate expression\n",
  "  ✅ L.bat    (mfu) — 6 isoforms, low expression\n",
  "  ❌ Human   (hsa) — not detected in milk exosomes\n",
  "  ❌ Rhesus  (mma) — not detected in milk exosomes\n\n",
  "━━━━━━━━━━━━━━━━━━━━━━━━━\n\n",
  "Key finding:\n",
  "circAcc1 is conserved in all\n",
  "3 bat species + mouse, but\n",
  "absent from primate milk exosomes.\n\n",
  "This may reflect lactation biology\n",
  "differences rather than genomic\n",
  "absence of the ACACA gene."
)

# Text-only figure
pdf("New-analysis/figures/circAcc1_conservation_summary.pdf", width = 6, height = 8)
par(mar = c(1, 1, 1, 1))
plot.new()
text(0.5, 0.5, conservation_text, cex = 0.9, family = "mono",
     adj = c(0.5, 0.5))
dev.off()
cat("  Saved: figures/circAcc1_conservation_summary.pdf\n")

cat("\n========== All figures generated ==========\n")
cat("Output: New-analysis/figures/\n")
cat("  circAcc1_detection_dotplot.pdf\n")
cat("  circAcc1_expression.pdf\n")
cat("  circAcc1_isoforms.pdf\n")
cat("  circAcc1_conservation_summary.pdf\n")
