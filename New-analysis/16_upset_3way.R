# ============================================================
# UpSet plot: ALL 3-species combinations with case-normalized genes
# Fixes: ACACA/Acaca/acaca → ACACA (same gene, different case)
# ============================================================

library(dplyr)
library(tidyr)
library(ComplexHeatmap)
library(grid)

setwd("D:/跨物种筛选保守circRNA")
dir.create("New-analysis/figures", showWarnings = FALSE, recursive = TRUE)

# ---- Load gene lookups and build case-NORMALIZED gene sets ----
gene_lookups <- readRDS("New-analysis/gene_lookups.rds")
SPECIES <- c("hsa", "mfu", "mma", "mmu", "mpi", "rsi")

cat("=== Step 1: Build case-normalized gene sets ===\n")

gene_sets_norm <- list()
for (sp in SPECIES) {
  genes <- unique(gene_lookups[[sp]]$gene_name)
  # Normalize to UPPERCASE
  genes_upper <- toupper(genes)
  gene_sets_norm[[sp]] <- unique(genes_upper)
  cat(sprintf("  %s: %d raw → %d after toupper (merged %d)\n",
              sp, length(genes), length(genes_upper),
              length(genes) - length(genes_upper)))
}

# ---- Verify ACACA case fix ----
cat("\n=== ACACA case-fix verification ===\n")
for (sp in SPECIES) {
  has_acaca <- "ACACA" %in% gene_sets_norm[[sp]]
  cat(sprintf("  %s: ACACA present = %s\n", sp, ifelse(has_acaca, "YES", "NO")))
}

# Check overlap before vs after
raw_mmu <- gene_lookups$mmu$gene_name
raw_hsa <- gene_lookups$hsa$gene_name
cat(sprintf("\nRaw mmu has 'Acaca': %s\n", "Acaca" %in% raw_mmu))
cat(sprintf("Raw hsa has 'ACACA': %s\n", "ACACA" %in% raw_hsa))
cat(sprintf("After toupper, 'ACACA' shared by hsa+mmu: %s\n",
            "ACACA" %in% gene_sets_norm$hsa && "ACACA" %in% gene_sets_norm$mmu))

# ---- Build combination matrix ----
cat("\n=== Step 2: Build combination matrix ===\n")

all_genes <- unique(unlist(gene_sets_norm))
cat(sprintf("Total unique genes (case-normalized): %d\n", length(all_genes)))

# Make binary matrix
comb_mat <- make_comb_mat(gene_sets_norm[SPECIES])
cat(sprintf("Combinations: %d\n", length(comb_size(comb_mat))))

# ---- Step 3: Plot focused on 3-way intersections ----
cat("\n=== Step 3: UpSet plot ===\n")

# Extract all 3-way combinations
three_way_idx <- comb_degree(comb_mat) == 3
three_way_sizes <- comb_size(comb_mat)[three_way_idx]
three_way_sizes <- sort(three_way_sizes, decreasing = TRUE)

cat(sprintf("3-way combinations (C(6,3)=20): %d non-empty\n",
            sum(three_way_sizes > 0)))
cat(sprintf("  Range: %d to %d genes\n", min(three_way_sizes), max(three_way_sizes)))
for (nm in names(head(three_way_sizes, 10))) {
  sp_in <- SPECIES[as.logical(as.integer(strsplit(nm, "")[[1]]))]
  cat(sprintf("  %s (%s): %d genes\n", nm, paste(sp_in, collapse="+"), three_way_sizes[nm]))
}

# Filter comb_mat to only 3-way
upset_3way <- comb_mat[comb_degree(comb_mat) == 3 & comb_size(comb_mat) > 0]

# Species colors
sp_colors <- c(
  hsa = "#E69F00", mfu = "#56B4E9", mma = "#009E73",
  mmu = "#F0E442", mpi = "#0072B2", rsi = "#D55E00"
)

# ---- Plot ----
pdf("New-analysis/figures/upset_3way.pdf", width = 14, height = 8)

ht <- UpSet(
  upset_3way,
  set_order = SPECIES,
  comb_order = order(comb_size(upset_3way), decreasing = TRUE),
  top_annotation = HeatmapAnnotation(
    "Shared genes" = anno_barplot(
      comb_size(upset_3way),
      border = FALSE,
      gp = gpar(fill = "#3498DB", col = "#3498DB"),
      height = unit(4, "cm"),
      axis_param = list(labels_rot = 0)
    ),
    annotation_name_side = "left",
    annotation_name_gp = gpar(fontsize = 10)
  ),
  right_annotation = upset_right_annotation(
    upset_3way,
    gp = gpar(fill = sp_colors[SPECIES]),
    annotation_name_gp = gpar(fontsize = 10),
    axis_param = list(labels_rot = 0)
  ),
  row_names_gp = gpar(fontsize = 10),
  column_title = paste0(
    "circRNA host genes: 3-species shared combinations\n",
    "(case-normalized: Acaca→ACACA; ",
    sum(comb_size(upset_3way) > 0), "/20 non-empty, ",
    "max=", max(comb_size(upset_3way)), " genes)"
  ),
  column_title_gp = gpar(fontsize = 13, fontface = "bold")
)
draw(ht)
decorate_annotation("Shared genes", {
  grid.text("Shared genes", x = unit(0.5, "npc"), y = unit(0.5, "npc"),
            gp = gpar(fontsize = 11, fontface = "bold"), rot = 90)
})
dev.off()

# ---- Also save as TIFF ----
tiff("New-analysis/figures/upset_3way.tiff", width = 3600, height = 2000, res = 300, compression = "lzw")
ht2 <- UpSet(
  upset_3way,
  set_order = SPECIES,
  comb_order = order(comb_size(upset_3way), decreasing = TRUE),
  top_annotation = HeatmapAnnotation(
    "Shared genes" = anno_barplot(
      comb_size(upset_3way),
      border = FALSE,
      gp = gpar(fill = "#3498DB", col = "#3498DB"),
      height = unit(4, "cm"),
      axis_param = list(labels_rot = 0)
    ),
    annotation_name_side = "left",
    annotation_name_gp = gpar(fontsize = 10)
  ),
  right_annotation = upset_right_annotation(
    upset_3way,
    gp = gpar(fill = sp_colors[SPECIES]),
    annotation_name_gp = gpar(fontsize = 10),
    axis_param = list(labels_rot = 0)
  ),
  row_names_gp = gpar(fontsize = 10),
  column_title = "circRNA host genes: 3-species shared (case-normalized)",
  column_title_gp = gpar(fontsize = 14, fontface = "bold")
)
draw(ht2)
dev.off()

# ---- Save top 3-way gene lists ----
cat("\n=== Step 4: Save top 3-way gene lists ===\n")

extract_comb_genes <- function(comb_mat_obj, comb_name) {
  idx <- which(names(comb_size(comb_mat_obj)) == comb_name)
  if (length(idx) == 0) return(character(0))
  genes <- extract_comb(comb_mat_obj, comb_name)
  sort(genes)
}

# Extract top 5 three-way intersections
top_3ways <- names(head(three_way_sizes, 5))
for (tw in top_3ways) {
  sp_in <- SPECIES[as.logical(as.integer(strsplit(tw, "")[[1]]))]
  genes <- extract_comb_genes(comb_mat, tw)
  cat(sprintf("\n%s (%s): %d genes\n", tw, paste(sp_in, collapse=" + "), length(genes)))
  cat(sprintf("  Top 20: %s\n", paste(head(genes, 20), collapse=", ")))
}

cat("\n=== Done ===\n")
cat("Output: New-analysis/figures/upset_3way.pdf\n")
