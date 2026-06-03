# ============================================================
# Diagnose: why is circACACA missing in hsa and mma?
# Investigate all detection levels from raw tool output
# ============================================================

library(dplyr)

source(if (file.exists("New-analysis/00_project_setup.R")) "New-analysis/00_project_setup.R" else "00_project_setup.R", encoding = "UTF-8")

cat("========== Diagnosing missing circACACA in hsa & mma ==========\n\n")

# ACACA gene identifiers across species:
# hsa: ACACA  (HGNC symbol), ENSG00000132142
# mma: ACACA  (HGNC symbol), ENSMMUG00000009349
# mmu: Acaca  (MGI symbol), ENSMUSG00000020532

for (sp in c("hsa", "mma")) {
  cat("========================================\n")
  cat("  ", sp, "\n")
  cat("========================================\n\n")

  # ---- 1. Check per-sample CIRCexplorer3 raw output ----
  cat("--- 1. Per-sample CIRCexplorer3 raw data ---\n")
  ce_files <- list.files(pattern = paste0("ce_", sp, "\\d\\.csv$"))
  for (f in ce_files) {
    ce <- read.csv(f, stringsAsFactors = FALSE)
    # Look for ACACA in geneName column (CIRCexplorer3 uses Ensembl gene IDs)
    if ("geneName" %in% names(ce)) {
      acaca_hits <- ce[grepl("ENSG00000132142|ENSMMUG00000009349|ACACA|Acaca",
                             ce$geneName, ignore.case = TRUE), ]
      if (nrow(acaca_hits) > 0) {
        cat(sprintf("  %s: %d circRNA(s) found!\n", f, nrow(acaca_hits)))
        cat(sprintf("    geneName: %s\n", paste(head(acaca_hits$geneName, 3), collapse=", ")))
        cat(sprintf("    readNumber: %s\n", paste(head(acaca_hits$readNumber, 3), collapse=", ")))
      }
    }
    # Also check by gene_id from the aggregated file
  }
  if (!exists("acaca_hits") || nrow(acaca_hits) == 0) {
    cat("  No ACACA hits in per-sample CIRCexplorer3 data\n")
  }

  # ---- 2. Check aggregated CIRCexplorer3 (ce_*.csv) ----
  cat("\n--- 2. Aggregated CIRCexplorer3 (ce_*.csv) ---\n")
  ce_agg <- read.csv(paste0("ce_", sp, ".csv"), stringsAsFactors = FALSE)
  # These only have circ_id + expression, no geneName
  # Need to cross-reference with lookup
  gene_lookups <- readRDS("New-analysis/gene_lookups.rds")
  sp_lookup <- gene_lookups[[sp]]
  acaca_circs <- sp_lookup$circ_id[grepl("ACACA|Acaca", sp_lookup$gene_name, ignore.case = TRUE)]
  cat(sprintf("  circ_ids mapped to ACACA in lookup: %d\n", length(acaca_circs)))
  if (length(acaca_circs) > 0) {
    ce_acaca <- ce_agg[ce_agg$circ_id %in% acaca_circs, ]
    cat(sprintf("  In aggregated ce_*.csv: %d rows\n", nrow(ce_acaca)))
    if (nrow(ce_acaca) > 0) {
      rep_cols <- grep(paste0("^", sp, "\\d"), names(ce_acaca), value = TRUE)
      for (r in rep_cols) {
        cat(sprintf("    %s: %s\n", r, paste(ce_acaca[[r]], collapse=", ")))
      }
    }
  }

  # ---- 3. Check per-sample CIRIquant raw GTF ----
  cat("\n--- 3. Per-sample CIRIquant GTF ---\n")
  re_path <- paste0("re_analysis/", sp, "/CIRIquant/")
  if (dir.exists(re_path)) {
    cq_files <- list.files(re_path, pattern = "\\.gtf$", full.names = TRUE)
    for (f in cq_files) {
      # Search for ACACA in GTF attributes using system grep (faster)
      cmd <- sprintf('findstr /c:"ACACA" /c:"Acaca" /c:"acaca" "%s"', f)
      result <- tryCatch(system(cmd, intern = TRUE), error = function(e) NULL)
      if (length(result) > 0) {
        cat(sprintf("  %s: %d lines with ACACA\n", basename(f), length(result)))
        for (i in seq_len(min(3, length(result)))) {
          # Extract circ_id and gene_name
          l <- result[i]
          circ_id_m <- regmatches(l, regexec('circ_id "([^"]+)"', l))[[1]]
          gene_name_m <- regmatches(l, regexec('gene_name "([^"]+)"', l))[[1]]
          bsj_m <- regmatches(l, regexec('bsj "([^"]+)"', l))[[1]]
          if (length(circ_id_m) > 1) {
            cat(sprintf("    circ_id=%s, gene_name=%s, bsj=%s\n",
                circ_id_m[2],
                if(length(gene_name_m)>1) gene_name_m[2] else "N/A",
                if(length(bsj_m)>1) bsj_m[2] else "N/A"))
          }
        }
      } else {
        cat(sprintf("  %s: no ACACA\n", basename(f)))
      }
    }
  }

  # ---- 4. Check per-sample find_circ raw BED ----
  cat("\n--- 4. Per-sample find_circ BED ---\n")
  re_path <- paste0("re_analysis/", sp, "/find_circ/")
  if (dir.exists(re_path)) {
    fc_files <- list.files(re_path, pattern = "\\.bed$", full.names = TRUE)
    for (f in fc_files) {
      # find_circ BED col 4 is circRNA ID (e.g., hsa_circ_000001)
      # col 17 is "NA" (no gene annotation in find_circ)
      # Find by genomic coordinate of known ACACA locus
      result <- read.table(f, header = FALSE, sep = "\t", stringsAsFactors = FALSE,
                          nrows = 1)
      cat(sprintf("  %s: %d columns, sample entry: %s\n",
          basename(f), ncol(result),
          paste(result[1, 1:4], collapse = "|")))
    }
  }

  # ---- 5. Check the *.cq.csv (merged CIRIquant with gene_name) ----
  cat("\n--- 5. Merged CQ file (*.cq.csv) ---\n")
  cq_merged <- paste0(sp, ".cq.csv")
  if (file.exists(cq_merged)) {
    cq <- read.csv(cq_merged, stringsAsFactors = FALSE)
    cat(sprintf("  File exists: %d rows\n", nrow(cq)))
    cat(sprintf("  Columns: %s\n", paste(names(cq), collapse=", ")))
    if ("gene_name" %in% names(cq)) {
      acaca_cq <- cq[grepl("ACACA|Acaca", cq$gene_name, ignore.case = TRUE), ]
      cat(sprintf("  ACACA rows: %d\n", nrow(acaca_cq)))
    }
    if ("gene_id" %in% names(cq)) {
      acaca_id <- cq[grepl("ENSG00000132142|ENSMMUG00000009349",
                           cq$gene_id, ignore.case = TRUE), ]
      cat(sprintf("  ACACA by gene_id (ENSG... / ENSMMUG...): %d\n", nrow(acaca_id)))
    }
  }

  cat("\n\n")
}

cat("========== Diagnosis complete ==========\n")
cat("\nRecommendations based on findings above:\n")
cat("  1. If raw CIRIquant GTF has ACACA → lower the tool intersection threshold\n")
cat("  2. If per-sample CE has ACACA → lower replicate count threshold (>=1 instead of >=2)\n")
cat("  3. If no tool detects ACACA → biological absence in milk exosomes\n")
