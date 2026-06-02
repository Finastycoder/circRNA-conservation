suppressPackageStartupMessages(library(ComplexHeatmap))
setwd("D:/跨物种筛选保守circRNA")

presence <- read.csv("New-analysis/gene_presence_matrix.csv", stringsAsFactors = FALSE)
SPECIES <- c("hsa", "mfu", "mma", "mmu", "mpi", "rsi")

genesets <- list()
for (sp in SPECIES) {
  genesets[[sp]] <- unique(toupper(presence[["gene_name"]][presence[[sp]] == TRUE]))
}
names(genesets) <- c("Human", "L.bat", "Rhesus", "Mouse", "B.bat", "H.bat")

cm <- make_comb_mat(genesets)

# 6-way genes
g6 <- extract_comb(cm, "111111")
cat("=== 6-way shared genes (", length(g6), ") ===\n", sep = "")
cat(paste(sort(g6), collapse = "\n"), "\n")

# 5-way (all except human)
g5 <- extract_comb(cm, "011111")
cat("\n=== 5-way (all except human, ", length(g5), " genes) ===\n", sep = "")
cat(paste(sort(g5), collapse = "\n"), "\n")
