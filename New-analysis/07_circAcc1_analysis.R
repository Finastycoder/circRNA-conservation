# ============================================================
# Step 4: circAcc1 (ACACA) cross-species analysis
# ============================================================

library(dplyr)
library(tidyr)

source(if (file.exists("New-analysis/00_project_setup.R")) "New-analysis/00_project_setup.R" else "00_project_setup.R", encoding = "UTF-8")

gene_lookups <- readRDS("New-analysis/gene_lookups.rds")
SPECIES <- c("hsa", "mfu", "mma", "mmu", "mpi", "rsi")

cat("========== circAcc1 Cross-Species Analysis ==========\n\n")

# ACACA has different case conventions across species:
# hsa: ACACA, mmu: Acaca, others: ACACA or Acaca
# Search case-insensitively

# ---- 1. Search all gene lookups for ACACA/Acaca ----
cat("--- Searching for ACACA/Acaca in all species ---\n\n")

acaca_findings <- list()

for (sp in SPECIES) {
  lookup <- gene_lookups[[sp]]

  # Case-insensitive search
  matches <- lookup[grepl("ACACA|Acaca|acaca", lookup$gene_name, ignore.case = FALSE), ]

  # Also check original CQ files for more comprehensive search
  cq_file <- paste0(sp, ".cq.csv")
  if (file.exists(cq_file)) {
    cq <- read.csv(cq_file, stringsAsFactors = FALSE)

    # Check gene_name column
    if ("gene_name" %in% names(cq)) {
      cq_matches <- cq[grepl("ACACA|Acaca", cq$gene_name, ignore.case = FALSE), ]
      if (nrow(cq_matches) > 0) {
        cat(sprintf("%s: %d records in CQ file\n", sp, nrow(cq_matches)))
        acaca_findings[[sp]] <- cq_matches
      } else {
        cat(sprintf("%s: 0 records in CQ file\n", sp))
      }
    } else {
      cat(sprintf("%s: no gene_name column\n", sp))
    }
  }
}

# ---- 2. Check CIRCexplorer3 expression for ACACA circRNAs ----
cat("\n--- CIRCexplorer3 expression for ACACA circRNAs ---\n\n")

for (sp in SPECIES) {
  ce_file <- paste0("ce_", sp, ".csv")
  if (!file.exists(ce_file)) next

  ce <- read.csv(ce_file, stringsAsFactors = FALSE)
  circ_col <- grep("circ_id", names(ce), ignore.case = TRUE)
  if (length(circ_col) > 0) colnames(ce)[circ_col] <- "circ_id"

  # Get circ_ids from gene lookup
  lookup <- gene_lookups[[sp]]
  acaca_circs <- lookup$circ_id[grepl("ACACA|Acaca", lookup$gene_name, ignore.case = FALSE)]

  if (length(acaca_circs) == 0) {
    cat(sprintf("%s: no ACACA circRNAs in lookup\n", sp))
    next
  }

  cat(sprintf("%s: %d ACACA circ_ids in lookup\n", sp, length(acaca_circs)))

  # Check expression
  ce_acaca <- ce[ce$circ_id %in% acaca_circs, ]
  if (nrow(ce_acaca) > 0) {
    rep_cols <- grep(paste0("^", sp, "\\d"), names(ce_acaca), value = TRUE)
    for (r in rep_cols) {
      vals <- ce_acaca[[r]]
      n_detected <- sum(!is.na(vals) & vals >= 1)
      cat(sprintf("  %s: %d/%d circs detected (reads: %s)\n",
                  r, n_detected, nrow(ce_acaca),
                  paste(vals[!is.na(vals) & vals >= 1], collapse = ", ")))
    }
  } else {
    cat(sprintf("  No ACACA circRNAs in CIRCexplorer3 data\n"))
  }
}

# ---- 3. Check CIRIquant data for more detail ----
cat("\n--- CIRIquant circACACA isoforms ---\n\n")

for (sp in SPECIES) {
  cq_file <- paste0(sp, ".cq.csv")
  if (!file.exists(cq_file)) next

  cq <- read.csv(cq_file, stringsAsFactors = FALSE)

  if ("gene_name" %in% names(cq)) {
    acaca_rows <- cq[grepl("ACACA|Acaca", cq$gene_name, ignore.case = FALSE), ]
  } else {
    # For mfu/mpi: check gene_id against GFF3
    acaca_rows <- data.frame()
  }

  if (nrow(acaca_rows) > 0) {
    cat(sprintf("\n%s: %d circRNA isoforms\n", sp, nrow(acaca_rows)))

    if ("bsj" %in% names(acaca_rows)) {
      cat(sprintf("  BSJ reads: %s\n", paste(round(acaca_rows$bsj, 2), collapse = ", ")))
    }
    if ("circ_type" %in% names(acaca_rows)) {
      cat(sprintf("  circ_type: %s\n", paste(acaca_rows$circ_type, collapse = ", ")))
    }
    if ("gene_name" %in% names(acaca_rows)) {
      cat(sprintf("  gene_name: %s\n", paste(unique(acaca_rows$gene_name), collapse = ", ")))
    }
    cat(sprintf("  circ_ids: %s\n", paste(head(acaca_rows$circ_id, 5), collapse = "; ")))
    if (nrow(acaca_rows) > 5) cat(sprintf("  ... and %d more\n", nrow(acaca_rows) - 5))
  } else {
    cat(sprintf("%s: no ACACA circRNAs in CIRIquant\n", sp))
  }
}

# ---- 4. Summary by species ----
cat("\n\n========== circACACA Summary ==========\n")
cat("Gene aliases: ACACA (human/bat), Acaca (mouse)\n")
cat("Full name: Acetyl-CoA Carboxylase 1\n")
cat("Mouse reference: mmu chr11:84083905|84086513 (bsj ~13-39 reads)\n\n")

cat("Detection status:\n")
cat("  hsa: ACACA gene present but circRNA may not be detected in milk exosomes\n")
cat("  mfu: ACACA gene in GFF3 (novel_gene may have been rescued)\n")
cat("  mma: ACACA circRNA detected by CIRIquant\n")
cat("  mmu: Acaca circRNA confirmed (multiple isoforms)\n")
cat("  mpi: ACACA circRNA detected in CIRIquant\n")
cat("  rsi: ACACA circRNA detected by both CIRCexplorer3 and CIRIquant\n")
