# ============================================================
# Step 2-5: Expression loading + Conservation analysis + circAcc1
# ============================================================

library(dplyr)
library(tidyr)

setwd("D:/跨物种筛选保守circRNA")

# ---- Load gene lookups from Part 1 ----
gene_lookups <- readRDS("New-analysis/gene_lookups.rds")
SPECIES <- c("hsa", "mfu", "mma", "mmu", "mpi", "rsi")

cat("========== Step 2: Load CIRCexplorer3 expression data ==========\n")

ce_data <- list()
for (sp in SPECIES) {
  f <- paste0("ce_", sp, ".csv")
  ce <- read.csv(f, stringsAsFactors = FALSE)
  # Fix column names: circ_id is the first named column
  circ_col <- grep("circ_id", names(ce), ignore.case = TRUE)
  if (length(circ_col) > 0) {
    colnames(ce)[circ_col] <- "circ_id"
  }
  # Get replicate columns
  rep_cols <- grep(paste0("^", sp, "\\d"), names(ce), value = TRUE)
  cat(sprintf("%s: %d circRNAs, %d replicates (%s)\n",
              sp, nrow(ce), length(rep_cols),
              paste(rep_cols, collapse = ", ")))
  ce_data[[sp]] <- list(df = ce, rep_cols = rep_cols)
}

# ---- Filter: expressed in >= 2 replicates with read >= 1 ----
cat("\n--- Filtering expressed circRNAs ---\n")

MIN_REPLICATES <- 2
MIN_READS <- 1

expressed_circs <- list()
expressed_genes <- list()

for (sp in SPECIES) {
  ce <- ce_data[[sp]]$df
  reps <- ce_data[[sp]]$rep_cols

  # Calculate n_reps with read >= MIN_READS
  expr_matrix <- ce[, reps, drop = FALSE]
  n_expressed <- rowSums(!is.na(expr_matrix) & expr_matrix >= MIN_READS, na.rm = TRUE)

  expressed <- ce[n_expressed >= MIN_REPLICATES, ]

  cat(sprintf("%s: %d/%d circRNAs expressed (>= %d reps with >= %d reads)\n",
              sp, nrow(expressed), nrow(ce), MIN_REPLICATES, MIN_READS))

  # Map to gene_names
  sp_lookup <- gene_lookups[[sp]]
  expr_ids <- expressed$circ_id

  sp_genes <- sp_lookup %>%
    filter(circ_id %in% expr_ids) %>%
    pull(gene_name) %>%
    unique()

  cat(sprintf("  → %d unique gene_names\n", length(sp_genes)))

  expressed_circs[[sp]] <- expressed
  expressed_genes[[sp]] <- sp_genes
}

# ---- Step 3: Gene-level conservation matrix ----
cat("\n========== Step 3: Gene conservation analysis ==========\n")

all_genes <- unique(unlist(expressed_genes))
cat(sprintf("Total unique gene_names across all species: %d\n", length(all_genes)))

presence_matrix <- data.frame(
  gene_name = all_genes,
  hsa = all_genes %in% expressed_genes$hsa,
  mfu = all_genes %in% expressed_genes$mfu,
  mma = all_genes %in% expressed_genes$mma,
  mmu = all_genes %in% expressed_genes$mmu,
  mpi = all_genes %in% expressed_genes$mpi,
  rsi = all_genes %in% expressed_genes$rsi,
  stringsAsFactors = FALSE
)

presence_matrix$n_species <- rowSums(presence_matrix[, SPECIES])

cat("\n=== Gene-level Conservation Summary ===\n")
for (n in 6:1) {
  count <- sum(presence_matrix$n_species == n)
  pct <- round(100 * count / nrow(presence_matrix), 2)
  cat(sprintf("%d species: %d genes (%.2f%%)\n", n, count, pct))
}

# Save full matrix
write.csv(presence_matrix, "New-analysis/gene_presence_matrix.csv", row.names = FALSE)
cat("\nSaved: New-analysis/gene_presence_matrix.csv\n")

# ---- Export conserved gene lists ----
for (n in 6:2) {
  genes_n <- presence_matrix %>%
    filter(n_species == n) %>%
    pull(gene_name)

  if (length(genes_n) > 0) {
    # Build detailed report
    report <- data.frame(gene_name = genes_n, stringsAsFactors = FALSE)
    for (sp in SPECIES) {
      report[[sp]] <- sapply(genes_n, function(g) {
        lookup <- gene_lookups[[sp]]
        circs <- lookup$circ_id[lookup$gene_name == g]
        expr <- expressed_circs[[sp]]
        circs_expr <- intersect(circs, expr$circ_id)
        paste(circs_expr, collapse = ";")
      })
      # Count
      report[[paste0(sp, "_n_circs")]] <- sapply(genes_n, function(g) {
        lookup <- gene_lookups[[sp]]
        circs <- lookup$circ_id[lookup$gene_name == g]
        expr <- expressed_circs[[sp]]
        length(intersect(circs, expr$circ_id))
      })
    }
    write.csv(report, paste0("New-analysis/conserved_genes_", n, "way.csv"), row.names = FALSE)
    cat(sprintf("Saved: New-analysis/conserved_genes_%dway.csv (%d genes)\n", n, length(genes_n)))
  } else {
    cat(sprintf("%d-way conserved: 0 genes\n", n))
  }
}

cat("\n========== Conservation analysis Part 1-2 complete ==========\n")
