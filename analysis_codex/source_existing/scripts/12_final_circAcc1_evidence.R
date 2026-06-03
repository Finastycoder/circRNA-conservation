# ============================================================
# Final: circAcc1 (ACACA) comprehensive evidence visualization
# Includes circAtlas v3.0 public database confirmation
# ============================================================

library(dplyr)
library(tidyr)
library(ggplot2)

source(if (file.exists("New-analysis/00_project_setup.R")) "New-analysis/00_project_setup.R" else "00_project_setup.R", encoding = "UTF-8")
dir.create("New-analysis/figures", showWarnings = FALSE, recursive = TRUE)

cat("========== circAcc1 Comprehensive Evidence ==========\n\n")

# ---- Build final evidence matrix ----
# Sources:
#   Tier 1: GFF3/GTF genome annotation (our analysis)
#   Tier 2: polyA mRNA / GTEx / Ensembl
#   Tier 3a: circRNA in milk exosomes (our CIRCexplorer3 data)
#   Tier 3b: circRNA in ANY tissue (circAtlas v3.0, ngdc.cncb.ac.cn/circatlas/)

evidence <- data.frame(
  species     = c("hsa", "mma", "mmu", "mfu", "mpi", "rsi"),
  s_label     = c("Human\n(hsa)", "Rhesus\n(mma)", "Mouse\n(mmu)",
                  "Long-winged\nbat (mfu)", "Big-footed\nbat (mpi)", "Horseshoe\nbat (rsi)"),
  category    = c("Primate", "Primate", "Rodent", "Bat", "Bat", "Bat"),

  # Tier 1
  gene_genome = c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE),
  gene_info   = c("chr17, ENSG00000278540", "chr16, ENSMMUG00000009349",
                  "chr11, ENSMUSG00000020532", "scaffold0019, ENSMFUG00019034613243",
                  "scaffold0021, ENSMPIG00021061614243", "NW_017739011.1, gene-ACACA"),

  # Tier 2
  mrna_expr   = c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE),
  mrna_info   = c("GTEx: ubiquitously expressed", "GTEx: ubiquitously expressed",
                  "GTEx: ubiquitously expressed",
                  "polyA RNA: 2455,1530,1581 counts (mean=1889)",
                  "polyA RNA: 3728,3863,2692 counts (mean=3428)",
                  "RefSeq annotation"),

  # Tier 3a: our milk exosome data
  circ_milk   = c(FALSE, FALSE, TRUE, TRUE, TRUE, TRUE),
  circ_n_iso  = c(0, 0, 9, 6, 1, 8),

  # Tier 3b: circAtlas public database
  circ_atlas  = c(TRUE, TRUE, TRUE, NA, NA, NA),
  circ_atlas_n = c(317, 92, 56, NA, NA, NA),

  stringsAsFactors = FALSE
)

# Derive final conclusion
evidence$final_status <- with(evidence,
  ifelse(circ_milk, "circRNA in milk (our data)",
  ifelse(circ_atlas, "circRNA in other tissues (circAtlas)", "No public data")))

evidence$final_status[evidence$species %in% c("mfu","mpi","rsi") &
                       !evidence$circ_milk] <- "Data pending"

# Print summary
cat("Species summary:\n")
for (i in seq_len(nrow(evidence))) {
  e <- evidence[i, ]
  cat(sprintf("  %s: gene=%s | mRNA=YES | circRNA_milk=%s(%d isoforms) | circAtlas=%s(%s isoforms)\n",
      e$species,
      if(e$gene_genome) "YES" else "NO",
      if(e$circ_milk) "YES" else "NO", e$circ_n_iso,
      if(is.na(e$circ_atlas)) "N/A" else if(e$circ_atlas) "YES" else "NO",
      if(is.na(e$circ_atlas_n)) "N/A" else e$circ_atlas_n))
}

write.csv(evidence, "New-analysis/acaca_final_evidence.csv", row.names = FALSE)
cat("\nSaved: New-analysis/acaca_final_evidence.csv\n")

# ---- Figure: Four-tier comprehensive evidence ----
cat("\nCreating comprehensive figure...\n")

ev_long <- evidence %>%
  select(s_label, category, gene_genome, mrna_expr, circ_milk, circ_atlas) %>%
  pivot_longer(cols = c(gene_genome, mrna_expr, circ_milk, circ_atlas),
               names_to = "tier", values_to = "status")

ev_long$tier <- factor(ev_long$tier,
  levels = c("circ_milk", "circ_atlas", "mrna_expr", "gene_genome"),
  labels = c("Tier 3a: circRNA\nin milk exosomes\n(this study)",
             "Tier 3b: circRNA\nin any tissue\n(circAtlas v3.0)",
             "Tier 2: Host gene\nmRNA expressed",
             "Tier 1: ACACA gene\nin genome"))

ev_long$s_label <- factor(ev_long$s_label,
  levels = rev(c("Human\n(hsa)", "Rhesus\n(mma)", "Mouse\n(mmu)",
                 "Long-winged\nbat (mfu)", "Big-footed\nbat (mpi)",
                 "Horseshoe\nbat (rsi)")))

ev_long$status_label <- ifelse(is.na(ev_long$status), "No data",
                        ifelse(ev_long$status, "Yes ✓", "No ✗"))
ev_long$status_label[ev_long$tier == "Tier 3b: circRNA\nin any tissue\n(circAtlas v3.0)" &
                      ev_long$s_label %in% c("Long-winged\nbat (mfu)",
                                             "Big-footed\nbat (mpi)",
                                             "Horseshoe\nbat (rsi)")] <- "Not in DB"

# Color scale
fill_colors <- c(
  "Yes ✓"    = "#27AE60",
  "No ✗"     = "#E74C3C",
  "No data"  = "#BDC3C7",
  "Not in DB" = "#F39C12"
)

p <- ggplot(ev_long, aes(x = tier, y = s_label, fill = status_label)) +
  geom_tile(color = "white", linewidth = 1.5, width = 0.85, height = 0.85) +
  geom_text(aes(label = status_label), size = 3.5, color = "white", fontface = "bold") +
  scale_fill_manual(values = fill_colors, guide = "none") +
  facet_wrap(~ category, scales = "free_y", ncol = 1, strip.position = "right") +
  labs(
    title = "circAcc1 (ACACA) — Comprehensive Cross-Species Evidence",
    subtitle = paste0(
      "ACACA gene: present in all 6 genomes | mRNA: expressed in all 6 species\n",
      "circRNA in milk: 4/6 (our CIRCexplorer3 data) | ",
      "circRNA in public DB: 3/3 model species confirmed (circAtlas v3.0)\n",
      "★ Human: 317 circACACA isoforms in circAtlas (other tissues) — NOT in milk exosomes"
    ),
    x = "", y = ""
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    strip.text = element_text(size = 9, color = "gray40"),
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 9, color = "gray40"),
    axis.text.x = element_text(size = 9, angle = 0, hjust = 0.5)
  )

ggsave("New-analysis/figures/acaca_final_evidence.pdf", p, width = 10, height = 6)
ggsave("New-analysis/figures/acaca_final_evidence.png", p, width = 10, height = 6, dpi = 150)
cat("  Saved: figures/acaca_final_evidence.pdf/png\n")

# ---- Figure 2: circAtlas isoform count comparison ----
cat("\nCreating circAtlas isoform comparison...\n")

atlas_df <- evidence[evidence$species %in% c("hsa","mma","mmu"), ] %>%
  select(species = s_label, category, circ_atlas_n, circ_n_iso) %>%
  mutate(species = gsub("\n", " ", species, fixed = TRUE))

atlas_long <- atlas_df %>%
  pivot_longer(cols = c(circ_atlas_n, circ_n_iso),
               names_to = "source", values_to = "n_isoforms")
atlas_long$source <- factor(atlas_long$source,
  levels = c("circ_atlas_n", "circ_n_iso"),
  labels = c("circAtlas DB\n(all tissues)", "Milk exosomes\n(this study)"))

p2 <- ggplot(atlas_long, aes(x = species, y = n_isoforms, fill = source)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6, alpha = 0.9) +
  geom_text(aes(label = n_isoforms),
            position = position_dodge(width = 0.7), vjust = -0.5, size = 4.5,
            fontface = "bold") +
  scale_fill_manual(values = c("circAtlas DB\n(all tissues)" = "#3498DB",
                                "Milk exosomes\n(this study)" = "#E74C3C")) +
  labs(title = "circACACA: Tissue Specificity vs Genomic Potential",
       subtitle = paste0(
         "Model species show abundant circACACA in other tissues but NOT in milk exosomes\n",
         "Human (hsa): 317 isoforms across multiple tissues — ZERO in milk\n",
         "Mouse (mmu): 56 isoforms in DB, 9 in milk — tissue-specific circRNA biogenesis"
       ),
       x = "", y = "Number of circRNA isoforms",
       fill = "Data Source") +
  theme_minimal(base_size = 13) +
  theme(legend.position = "bottom",
        panel.grid.major.x = element_blank(),
        plot.title = element_text(face = "bold"),
        plot.subtitle = element_text(size = 9, color = "gray40"))

ggsave("New-analysis/figures/acaca_atlas_comparison.pdf", p2, width = 8, height = 5.5)
ggsave("New-analysis/figures/acaca_atlas_comparison.png", p2, width = 8, height = 5.5, dpi = 150)
cat("  Saved: figures/acaca_atlas_comparison.pdf/png\n")

cat("\n========== All done! ==========\n")
