# ============================================================
# UpSet plot: 6-species intersections — degree >= 3 only
# Case-normalized gene names. No 1-way or 2-way clutter.
# ============================================================

suppressPackageStartupMessages({
  library(ComplexHeatmap)
  library(dplyr)
})
library(grid)

setwd("D:/跨物种筛选保守circRNA")
dir.create("New-analysis/figures", showWarnings = FALSE, recursive = TRUE)

# ---- Load and normalize ----
presence <- read.csv("New-analysis/gene_presence_matrix.csv", stringsAsFactors = FALSE)
SPECIES <- c("hsa", "mfu", "mma", "mmu", "mpi", "rsi")

gene_sets <- list()
for (sp in SPECIES) {
  genes <- presence$gene_name[presence[[sp]] == TRUE]
  gene_sets[[sp]] <- unique(toupper(genes))
}

cat("Genes per species (case-normalized):\n")
for (sp in SPECIES) cat(sprintf("  %s: %d\n", sp, length(gene_sets[[sp]])))

# ---- Build combination matrix ----
set_names <- c("Human (hsa)", "L.bat (mfu)", "Rhesus (mma)",
               "Mouse (mmu)", "B.bat (mpi)", "H.bat (rsi)")
names(gene_sets) <- set_names

comb_mat <- make_comb_mat(gene_sets)

# Filter: only degree >= 3
comb_filt <- comb_mat[comb_degree(comb_mat) >= 3 & comb_size(comb_mat) > 0]

n_total <- sum(comb_size(comb_filt) > 0)
cat(sprintf("\nIntersections with >= 3 species: %d (total genes shown)\n", n_total))

# Print summary
deg_sizes <- sort(comb_size(comb_filt), decreasing = TRUE)
for (nm in names(deg_sizes)) {
  sp_in <- set_names[as.logical(as.integer(strsplit(nm, "")[[1]]))]
  cat(sprintf("  %s (%d-way, %s): %d genes\n",
      nm, sum(as.integer(strsplit(nm, "")[[1]])),
      paste(sp_in, collapse = " + "),
      deg_sizes[nm]))
}

# ---- Colors ----
sp_colors <- c(
  "Human (hsa)" = "#E69F00", "L.bat (mfu)" = "#56B4E9",
  "Rhesus (mma)" = "#009E73", "Mouse (mmu)" = "#F0E442",
  "B.bat (mpi)" = "#0072B2", "H.bat (rsi)" = "#D55E00"
)

# Color by degree: 6=dark red, 5=orange, 4=green, 3=blue
deg_colors <- c("6" = "#C0392B", "5" = "#E67E22", "4" = "#27AE60", "3" = "#3498DB")
deg <- comb_degree(comb_filt)
bar_colors <- deg_colors[as.character(deg)]

# ---- Plot ----
pdf("New-analysis/figures/upset_3plus.pdf", width = 12, height = 7)

ht <- UpSet(
  comb_filt,
  set_order = rev(set_names),
  comb_order = order(comb_size(comb_filt), decreasing = TRUE),
  top_annotation = HeatmapAnnotation(
    "Shared genes" = anno_barplot(
      comb_size(comb_filt),
      border = FALSE,
      gp = gpar(fill = bar_colors, col = bar_colors),
      height = unit(5, "cm"),
      axis_param = list(labels_rot = 0)
    ),
    annotation_name_side = "left",
    annotation_name_gp = gpar(fontsize = 10)
  ),
  right_annotation = upset_right_annotation(
    comb_filt,
    gp = gpar(fill = sp_colors[set_names]),
    annotation_name_gp = gpar(fontsize = 9),
    axis_param = list(labels_rot = 0)
  ),
  row_names_gp = gpar(fontsize = 9),
  column_title = paste0(
    "circRNA host gene conservation (≥3 species)\n",
    "Case-normalized (Acaca→ACACA). ",
    "6-way: 0 genes. 5-way: 1 gene. 4-way: 14 genes. 3-way: 56 genes."
  ),
  column_title_gp = gpar(fontsize = 13, fontface = "bold")
)
draw(ht)

# Legend for degree colors
lgd <- Legend(
  labels = c("6 species", "5 species", "4 species", "3 species"),
  legend_gp = gpar(fill = deg_colors),
  title = "Conservation level",
  title_gp = gpar(fontsize = 10, fontface = "bold"),
  labels_gp = gpar(fontsize = 9)
)
draw(lgd, x = unit(0.82, "npc"), y = unit(0.85, "npc"))

dev.off()

# TIFF
tiff("New-analysis/figures/upset_3plus.tiff", width = 3000, height = 1800, res = 300, compression = "lzw")
ht <- UpSet(
  comb_filt,
  set_order = rev(set_names),
  comb_order = order(comb_size(comb_filt), decreasing = TRUE),
  top_annotation = HeatmapAnnotation(
    "Shared genes" = anno_barplot(
      comb_size(comb_filt),
      border = FALSE,
      gp = gpar(fill = bar_colors, col = bar_colors),
      height = unit(5, "cm"),
      axis_param = list(labels_rot = 0)
    ),
    annotation_name_side = "left",
    annotation_name_gp = gpar(fontsize = 10)
  ),
  right_annotation = upset_right_annotation(
    comb_filt,
    gp = gpar(fill = sp_colors[set_names]),
    annotation_name_gp = gpar(fontsize = 9),
    axis_param = list(labels_rot = 0)
  ),
  row_names_gp = gpar(fontsize = 9),
  column_title = paste0(
    "circRNA host gene conservation (>=3 species) | ",
    "Case-normalized | 6-way=0, 5-way=1, 4-way=14, 3-way=56"
  ),
  column_title_gp = gpar(fontsize = 13, fontface = "bold")
)
draw(ht)
draw(lgd, x = unit(0.82, "npc"), y = unit(0.85, "npc"))
dev.off()

cat("\nDone! upset_3plus.pdf / upset_3plus.tiff\n")
