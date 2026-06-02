# ============================================================
# Prepare combined BLAST database from Ensembl proteomes
# Reformats headers to include species prefix for easy parsing
# ============================================================

library(dplyr)

setwd("D:/跨物种筛选保守circRNA")
dir.create("New-analysis/blast", showWarnings = FALSE, recursive = TRUE)

# ---- Preprocess each reference FASTA ----
cat("========== Preprocessing reference proteomes ==========\n")

ref_files <- list(
  hsa = "New-analysis/blast/hsa_reference_proteins.faa",
  mma = "New-analysis/blast/mma_reference_proteins.faa",
  mmu = "New-analysis/blast/mmu_reference_proteins.faa"
)

all_out <- "New-analysis/blast/all_reference_proteins.faa"
file.create(all_out)

stats <- data.frame(sp = character(), n_proteins = integer(),
                    n_with_symbol = integer(), stringsAsFactors = FALSE)

for (sp in names(ref_files)) {
  f <- ref_files[[sp]]
  if (!file.exists(f)) {
    cat("WARNING: File not found:", f, "\n")
    next
  }

  cat("\nProcessing:", sp, "(", f, ")\n")

  con_in <- file(f, "r")
  con_out <- file(all_out, "a")

  n_prot <- 0
  n_sym <- 0

  while (TRUE) {
    lines <- readLines(con_in, n = 10000)
    if (length(lines) == 0) break

    for (l in lines) {
      if (grepl("^>", l)) {
        n_prot <- n_prot + 1

        # Extract gene_symbol from Ensembl header
        gs_match <- regmatches(l, regexec("gene_symbol:([^ ]+)", l))[[1]]
        gene_symbol <- if (length(gs_match) > 1) gs_match[2] else "UNKNOWN"

        # Extract gene ID
        gid_match <- regmatches(l, regexec("gene:([^ ]+)", l))[[1]]
        gene_id <- if (length(gid_match) > 1) gid_match[2] else "UNKNOWN"

        # Extract protein ID (first token without >)
        pid_match <- regmatches(l, regexec("^>([^ ]+)", l))[[1]]
        prot_id <- if (length(pid_match) > 1) pid_match[2] else "UNKNOWN"

        if (gene_symbol != "UNKNOWN") n_sym <- n_sym + 1

        # New header: >sp|gene_symbol|gene_id|prot_id
        new_header <- paste0(">", sp, "|", gene_symbol, "|", gene_id, "|", prot_id)
        writeLines(new_header, con_out)
      } else {
        # Sequence line — write as-is
        writeLines(l, con_out)
      }
    }
  }

  close(con_in)
  close(con_out)

  cat(sprintf("  %d proteins, %d with gene_symbol (%.1f%%)\n",
              n_prot, n_sym, 100*n_sym/max(1,n_prot)))

  stats <- rbind(stats, data.frame(
    sp = sp, n_proteins = n_prot, n_with_symbol = n_sym,
    stringsAsFactors = FALSE
  ))
}

cat("\n========== Combined reference FASTA ==========\n")
all_size <- file.info(all_out)$size
cat(sprintf("File: %s (%.1f MB)\n", all_out, all_size / 1e6))

total_prots <- sum(stats$n_proteins)
total_sym <- sum(stats$n_with_symbol)
cat(sprintf("Total: %d proteins, %d with gene_symbol (%.1f%%)\n",
            total_prots, total_sym, 100*total_sym/total_prots))

# ---- Build BLAST database ----
cat("\n========== Building BLAST database ==========\n")

blast_bin <- "tools/ncbi-blast-2.17.0+/bin"
db_out <- "New-analysis/blast/all_ref_db"

makeblastdb <- file.path(blast_bin, "makeblastdb.exe")
cmd <- sprintf('"%s" -in "%s" -dbtype prot -out "%s" -title "RefProteome"',
               makeblastdb, all_out, db_out)
cat("Running makeblastdb...\n")
system(cmd)

cat("\n=== Database preparation complete ===\n")
cat("DB files: New-analysis/blast/all_ref_db.phr/.pin/.psq\n")
