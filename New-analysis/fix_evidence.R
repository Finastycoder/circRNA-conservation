# Quick fix: update ACACA evidence matrix with real mRNA data + circBase result
library(dplyr)

source(if (file.exists("New-analysis/00_project_setup.R")) "New-analysis/00_project_setup.R" else "00_project_setup.R", encoding = "UTF-8")

# Read current evidence
evidence <- read.csv("New-analysis/acaca_evidence_matrix.csv", stringsAsFactors = FALSE)
cat("Current evidence:\n")
print(evidence[, c("species", "gene_in_genome", "mRNA_evidence", "circRNA_detected")])

# Update with real mRNA data
# mfu: polyA RNA featureCounts — ENSMFUG00019034613243 = ACACA, counts: 2455, 1530, 1581 (mean 1889)
# mpi: polyA RNA featureCounts — ENSMPIG00021061614243 = ACACA, counts: 3728, 3863, 2692 (mean 3428)
# hsa/mma/mmu: known expressed genes
# rsi: gene annotated in RefSeq
# circBase: confirmed NO circACACA in human database

evidence$mRNA_evidence <- c("Yes", "Yes", "Yes", "Yes", "Yes", "Yes")
evidence$mRNA_detail <- c(
  "GTEx/Ensembl — ubiquitously expressed",
  "GTEx/Ensembl — ubiquitously expressed",
  "GTEx/Ensembl — ubiquitously expressed",
  "polyA RNA: 2455, 1530, 1581 counts (mean=1889, 3 samples)",
  "polyA RNA: 3728, 3863, 2692 counts (mean=3428, 3 samples)",
  "RefSeq gene annotation (gene-ACACA)"
)

evidence$circbase_evidence <- c(
  "NOT in circBase — no human circACACA",
  "Not queried",
  "Not queried",
  "Not applicable (non-model)",
  "Not applicable (non-model)",
  "Not applicable (non-model)"
)

write.csv(evidence, "New-analysis/acaca_evidence_matrix.csv", row.names = FALSE)

cat("\nUpdated evidence:\n")
print(evidence[, c("species", "gene_in_genome", "mRNA_evidence", "circRNA_detected")])
cat("\nDone.\n")
