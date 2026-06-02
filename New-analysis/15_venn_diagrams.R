# ============================================================
# Multi-species conservation — Venn diagrams by level
# 2-way to 5-way Venn + UpSet for 6-species overview
# ============================================================

library(dplyr)
library(ggplot2)
library(grid)
library(VennDiagram)
library(UpSetR)

setwd("D:/跨物种筛选保守circRNA")
dir.create("New-analysis/figures", showWarnings = FALSE, recursive = TRUE)

# ---- Load data ----
presence <- read.csv("New-analysis/gene_presence_matrix.csv", stringsAsFactors = FALSE)
SPECIES <- c("hsa", "mfu", "mma", "mmu", "mpi", "rsi")

cat("Loaded:", nrow(presence), "genes\n")
for (sp in SPECIES) {
  n <- sum(presence[[sp]])
  cat(sprintf("  %s: %d genes\n", sp, n))
}

# Species sets — CASE NORMALIZED (toupper) to fix Acaca/ACACA mismatch
gene_sets <- list()
for (sp in SPECIES) {
  raw_genes <- presence$gene_name[presence[[sp]] == TRUE]
  gene_sets[[sp]] <- unique(toupper(raw_genes))
  n_merged <- length(raw_genes) - length(gene_sets[[sp]])
  cat(sprintf("  %s: %d genes (merged %d by toupper)\n", sp, length(gene_sets[[sp]]), n_merged))
}

# Verify ACACA fix
cat("\n  ACACA present:")
for (sp in SPECIES) cat(sprintf(" %s=%s", sp, "ACACA" %in% gene_sets[[sp]]))
cat("\n")

# ---- Okabe-Ito species colors ----
sp_colors <- c(
  hsa = "#E69F00", mfu = "#56B4E9", mma = "#009E73",
  mmu = "#F0E442", mpi = "#0072B2", rsi = "#D55E00"
)
sp_labels <- c(
  hsa = "Human\n(hsa)", mfu = "L-winged bat\n(mfu)", mma = "Rhesus\n(mma)",
  mmu = "Mouse\n(mmu)", mpi = "B-winged bat\n(mpi)", rsi = "Horseshoe bat\n(rsi)"
)

# ---- draw_venn helper ----
draw_venn <- function(species_set, main_title, filename) {
  n <- length(species_set)
  gene_list <- gene_sets[species_set]
  fill_cols <- sp_colors[species_set]
  disp_labels <- sp_labels[species_set]

  cat(sprintf("  %d-way: %s\n", n, paste(species_set, collapse=", ")))

  tiff(paste0("New-analysis/figures/venn_", filename, ".tiff"),
       width = 2400, height = 2400, res = 300, compression = "lzw")

  vp <- venn.diagram(
    x = gene_list,
    category.names = disp_labels,
    filename = NULL, output = TRUE,

    lwd = 2.5, lty = "solid",
    col = fill_cols,
    fill = fill_cols, alpha = 0.35,

    cex = 1.4, fontface = "bold", fontfamily = "sans",

    cat.cex = 1.1, cat.fontface = "bold", cat.fontfamily = "sans",
    cat.col = fill_cols,

    main = main_title, main.cex = 1.5,
    main.fontface = "bold", main.fontfamily = "sans",

    margin = 0.08,
    disable.logging = TRUE
  )
  grid.draw(vp)
  dev.off()

  pdf(paste0("New-analysis/figures/venn_", filename, ".pdf"),
      width = 9, height = 9)
  grid.draw(vp)
  dev.off()
}

# ============================
# 6-species: UpSet plot (more readable than 6-circle Venn)
# ============================
cat("\n=== 6-species: UpSet plot ===\n")

gene_list_upset <- list(
  `Human (hsa)` = gene_sets$hsa,
  `Rhesus (mma)` = gene_sets$mma,
  `Mouse (mmu)` = gene_sets$mmu,
  `L.bat (mfu)` = gene_sets$mfu,
  `B.bat (mpi)` = gene_sets$mpi,
  `H.bat (rsi)` = gene_sets$rsi
)

pdf("New-analysis/figures/venn_6way_upset.pdf", width = 14, height = 8)
upset(
  fromList(gene_list_upset),
  nsets = 6,
  nintersects = 25,
  order.by = "freq",
  main.bar.color = "#3498DB",
  sets.bar.color = c(sp_colors["rsi"], sp_colors["mpi"], sp_colors["mmu"],
                      sp_colors["mma"], sp_colors["mfu"], sp_colors["hsa"]),
  matrix.color = "#2C3E50",
  mainbar.y.label = "Number of shared genes",
  sets.x.label = "Genes per species",
  text.scale = c(1.5, 1.3, 1.3, 1.2, 1.5, 1.3)
)
dev.off()

tiff("New-analysis/figures/venn_6way_upset.tiff",
     width = 3600, height = 2000, res = 300, compression = "lzw")
upset(
  fromList(gene_list_upset),
  nsets = 6, nintersects = 25, order.by = "freq",
  main.bar.color = "#3498DB",
  sets.bar.color = c(sp_colors["rsi"], sp_colors["mpi"], sp_colors["mmu"],
                      sp_colors["mma"], sp_colors["mfu"], sp_colors["hsa"]),
  matrix.color = "#2C3E50",
  mainbar.y.label = "Number of shared genes",
  sets.x.label = "Genes per species",
  text.scale = c(1.5, 1.3, 1.3, 1.2, 1.5, 1.3)
)
dev.off()
cat("  Saved: venn_6way_upset.pdf/tiff\n")

# ============================
# 5-way Venn
# ============================
cat("\n=== 5-way Venn diagrams ===\n")

# 5-way (all except mmu — as found in conservation analysis)
five_1 <- c("hsa", "mfu", "mma", "mpi", "rsi")
draw_venn(five_1,
  "5 species shared (all except mouse — 10 genes)",
  "5way_no_mmu")

# 5-way (all except hsa)
five_2 <- c("mma", "mfu", "mmu", "mpi", "rsi")
draw_venn(five_2,
  "5 species shared (all except human)",
  "5way_no_hsa")

# ============================
# 4-way Venn (selected subsets)
# ============================
cat("\n=== 4-way Venn diagrams ===\n")

# Primates + Bats (hsa, mma, mfu, rsi)
draw_venn(c("hsa", "mma", "mfu", "rsi"),
  "Primates + Bats (4 species — 105 genes at >=4-way)",
  "4way_primates_bats")

# Mouse + 3 Bats
draw_venn(c("mmu", "mfu", "mpi", "rsi"),
  "Mouse + 3 Bats",
  "4way_mouse_bats")

# Model organisms + rsi
draw_venn(c("hsa", "mma", "mmu", "rsi"),
  "Model organisms + Horseshoe bat",
  "4way_models_rsi")

# ============================
# 3-way Venn
# ============================
cat("\n=== 3-way Venn diagrams ===\n")

# 3 model organisms
draw_venn(c("hsa", "mma", "mmu"),
  "Model organisms (human, macaque, mouse)",
  "3way_models")

# 3 bats
draw_venn(c("mfu", "mpi", "rsi"),
  "Three bat species",
  "3way_bats")

# Human + Mouse + R.bat (widest evolutionary span)
draw_venn(c("hsa", "mmu", "rsi"),
  "Human · Mouse · Horseshoe bat",
  "3way_hsa_mmu_rsi")

# ============================
# 2-way Venn (key pairs)
# ============================
cat("\n=== 2-way Venn diagrams ===\n")

pairs <- list(
  list(c("hsa", "mma"), "Primates (Human vs Rhesus)"),
  list(c("hsa", "mmu"), "Human vs Mouse"),
  list(c("mma", "mmu"), "Rhesus vs Mouse"),
  list(c("mfu", "mpi"), "Two bats: Long-winged vs Big-footed"),
  list(c("hsa", "rsi"), "Human vs Horseshoe bat"),
  list(c("mmu", "rsi"), "Mouse vs Horseshoe bat"),
  list(c("mpi", "rsi"), "Big-footed vs Horseshoe bat")
)

for (pair_info in pairs) {
  draw_venn(pair_info[[1]], pair_info[[2]],
    paste0("2way_", paste(pair_info[[1]], collapse="_")))
}

# ============================
# Done
# ============================
cat("\n========== Complete ==========\n")
cat("Output: New-analysis/figures/venn_*.pdf (17 diagrams)\n")
