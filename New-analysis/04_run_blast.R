# ============================================================
# Step 0 Phase 5: BLAST novel proteins against reference + annotation rescue
# ============================================================
# Header format: >sp|gene_symbol|gene_id|prot_id
# Query format:  >gene_id|novel_name=XXX|aa_len=N|cds_len=N
# ============================================================

library(dplyr)

setwd("D:/跨物种筛选保守circRNA")

BLAST_BIN <- "tools/ncbi-blast-2.17.0+/bin"
BLAST_DIR <- "New-analysis/blast"
BLAST_DB   <- "D:/blast_work/all_ref_db"   # Use ASCII path (no Chinese chars)

# ---- Step 1: Run BLASTp ----
cat("========== Running BLASTp ==========\n")

blastp <- file.path(BLAST_BIN, "blastp.exe")

for (sp in c("mfu", "mpi")) {
  query_file <- file.path(BLAST_DIR, paste0(sp, "_novel_proteins.faa"))
  if (!file.exists(query_file)) {
    cat("Query file not found:", query_file, "\n")
    next
  }

  cat("\nBLASTing", sp, "novel proteins...\n")
  out_file <- file.path(BLAST_DIR, paste0(sp, "_blast_results.txt"))

  cmd <- sprintf(
    '"%s" -query "%s" -db "%s" -out "%s" -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore" -evalue 1e-5 -num_threads 4 -max_target_seqs 5',
    blastp, query_file, BLAST_DB, out_file
  )

  ret <- system(cmd)
  if (ret == 0) {
    cat(sprintf("  BLAST complete: %s\n", out_file))
  } else {
    cat(sprintf("  BLAST ERROR: return code %d\n", ret))
  }
}

# ---- Step 2: Parse BLAST results ----
cat("\n========== Parsing BLAST results ==========\n")

parse_blast <- function(sp) {
  blast_file <- file.path(BLAST_DIR, paste0(sp, "_blast_results.txt"))
  if (!file.exists(blast_file)) {
    cat("  BLAST results not found:", blast_file, "\n")
    return(NULL)
  }

  col_names <- c("qseqid", "sseqid", "pident", "length", "mismatch",
                 "gapopen", "qstart", "qend", "sstart", "send",
                 "evalue", "bitscore")

  results <- read.table(blast_file, header = FALSE, sep = "\t",
                        stringsAsFactors = FALSE, comment.char = "",
                        col.names = col_names)
  cat(sprintf("  %s: %d raw BLAST hits\n", sp, nrow(results)))

  # Parse query ID: gene_id (before first |)
  results$query_gene_id <- sub("\\|.*", "", results$qseqid)

  # Parse subject: sp|gene_symbol|gene_id|prot_id
  s_parts <- strsplit(results$sseqid, "\\|")
  results$subject_species   <- sapply(s_parts, `[`, 1)
  results$subject_gene_symbol <- sapply(s_parts, `[`, 2)
  results$subject_gene_id   <- sapply(s_parts, `[`, 3)

  # Take best hit per query (max bitscore)
  results <- results %>%
    group_by(query_gene_id) %>%
    slice_max(order_by = bitscore, n = 1) %>%
    ungroup()

  cat(sprintf("  %d unique queries with top hit\n", nrow(results)))
  return(results)
}

blast_mfu <- parse_blast("mfu")
blast_mpi <- parse_blast("mpi")

# ---- Step 3: Build novel_gene → real_gene_name mapping ----
cat("\n========== Building annotation rescue mapping ==========\n")

IDENTITY_THRESHOLD <- 50  # minimum % identity

build_mapping <- function(blast_df, sp_name, novel_genes_file) {
  if (is.null(blast_df)) return(NULL)

  novel_genes <- read.csv(novel_genes_file, stringsAsFactors = FALSE)

  mapping <- blast_df %>%
    filter(pident >= IDENTITY_THRESHOLD) %>%
    select(
      gene_id          = query_gene_id,
      new_gene_name    = subject_gene_symbol,
      identity         = pident,
      evalue           = evalue,
      bitscore         = bitscore,
      ref_species      = subject_species,
      ref_gene_id      = subject_gene_id
    )

  # Merge with original novel_gene names
  mapping <- mapping %>%
    left_join(novel_genes %>% select(gene_id, original_name = gene_name),
              by = "gene_id")

  cat(sprintf("  %s: %d novel genes rescued (identity >= %d%%)\n",
              sp_name, nrow(mapping), IDENTITY_THRESHOLD))
  return(mapping)
}

mapping_mfu <- build_mapping(blast_mfu, "mfu", "New-analysis/mfu_novel_gene_list.csv")
mapping_mpi <- build_mapping(blast_mpi, "mpi", "New-analysis/mpi_novel_gene_list.csv")

# Save
if (!is.null(mapping_mfu)) {
  write.csv(mapping_mfu, "New-analysis/blast/mfu_novel_gene_mapping.csv", row.names = FALSE)
}
if (!is.null(mapping_mpi)) {
  write.csv(mapping_mpi, "New-analysis/blast/mpi_novel_gene_mapping.csv", row.names = FALSE)
}

# ---- Step 4: Summary ----
cat("\n========== Summary ==========\n")
for (sp in c("mfu", "mpi")) {
  mapping <- get(paste0("mapping_", sp))
  if (is.null(mapping)) next
  novel_total <- nrow(read.csv(paste0("New-analysis/", sp, "_novel_gene_list.csv")))
  rescued <- nrow(mapping)
  cat(sprintf("\n%s: %d/%d novel genes rescued (%.1f%%)\n",
              sp, rescued, novel_total, 100*rescued/novel_total))
  cat(sprintf("  Top rescued genes:\n    %s\n",
              paste(head(unique(mapping$new_gene_name), 20), collapse = ", ")))
}

cat("\n=== Annotation rescue complete ===\n")
