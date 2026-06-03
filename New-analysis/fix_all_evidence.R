# ============================================================
# Update all evidence CSVs with case-normalized data
# ============================================================

source(if (file.exists("New-analysis/00_project_setup.R")) "New-analysis/00_project_setup.R" else "00_project_setup.R", encoding = "UTF-8")

# 1. Update acaca_evidence_matrix.csv
evidence <- read.csv("New-analysis/acaca_evidence_matrix.csv", stringsAsFactors = FALSE)
evidence$mRNA_evidence <- "Yes"
write.csv(evidence, "New-analysis/acaca_evidence_matrix.csv", row.names = FALSE)
cat("Updated: acaca_evidence_matrix.csv\n")

# 2. Regenerate acaca_final_evidence.csv
evidence2 <- data.frame(
  species     = c("hsa", "mma", "mmu", "mfu", "mpi", "rsi"),
  s_label     = c("Human (hsa)", "Rhesus (mma)", "Mouse (mmu)",
                  "L.bat (mfu)", "B.bat (mpi)", "H.bat (rsi)"),
  category    = c("Primate", "Primate", "Rodent", "Bat", "Bat", "Bat"),
  gene_genome = c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE),
  mrna_expr   = c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE),
  circ_milk   = c(FALSE, FALSE, TRUE, TRUE, TRUE, TRUE),
  circ_n_iso  = c(0, 0, 9, 6, 1, 8),
  circ_atlas  = c(TRUE, TRUE, TRUE, NA, NA, NA),
  circ_atlas_n = c(317, 92, 56, NA, NA, NA),
  stringsAsFactors = FALSE
)
write.csv(evidence2, "New-analysis/acaca_final_evidence.csv", row.names = FALSE)
cat("Regenerated: acaca_final_evidence.csv (case-normalized)\n")

# 3. Print updated conservation stats
presence <- read.csv("New-analysis/gene_presence_matrix.csv", stringsAsFactors = FALSE)
cat("\n=== Updated Conservation Stats ===\n")
cat("Total unique genes (case-normalized):", nrow(presence), "\n")
for (n in 6:1) {
  count <- sum(presence$n_species == n)
  cat(sprintf("  %d-way: %d genes\n", n, count))
}

cat("\n=== 6-way conserved genes ===\n")
g6 <- presence$gene_name[presence$n_species == 6]
cat(paste(sort(g6), collapse = "\n"), "\n")

cat("\nDone!\n")
