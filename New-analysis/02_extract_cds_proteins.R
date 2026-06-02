# ============================================================
# Step 0 Phase 2-3: Efficient CDS extraction + translation
# Uses GRanges + getSeq (batch mode) instead of gene-by-gene loop
# ============================================================

library(Biostrings)
library(GenomicRanges)
library(rtracklayer)
library(dplyr)

setwd("D:/跨物种筛选保守circRNA")
dir.create("New-analysis/blast", showWarnings = FALSE, recursive = TRUE)

# ---- Load novel_gene lists from Phase 1 ----
mfu_novel_genes <- read.csv("New-analysis/mfu_novel_gene_list.csv", stringsAsFactors = FALSE)
mpi_novel_genes <- read.csv("New-analysis/mpi_novel_gene_list.csv", stringsAsFactors = FALSE)

cat("mfu novel genes:", nrow(mfu_novel_genes), "\n")
cat("mpi novel genes:", nrow(mpi_novel_genes), "\n")

# ============================================================
# Efficient approach: Import GFF3 CDS features + filter by gene
# ============================================================

process_species <- function(sp, novel_genes) {
  cat("\n========== Processing", sp, "==========\n")

  gff3_path <- paste0("genome and annotation/", sp, ".gff3")
  genome_path <- paste0("genome and annotation/", sp, ".fa")

  novel_ids <- unique(novel_genes$gene_id)
  cat("Novel gene IDs:", length(novel_ids), "\n")

  # ---- Step 1: Get CDS GRanges from GFF3 for novel genes ----
  # Parse CDS lines from GFF3, keep only those linked to novel genes
  cat("Parsing GFF3 for CDS features linked to novel genes...\n")

  con <- file(gff3_path, "r")
  novel_cds_gr <- GRanges()
  mrna_to_gene <- list()

  # First pass: build mRNA→gene mapping
  cat("  Pass 1: Building mRNA→gene mapping...\n")
  while (TRUE) {
    lines <- readLines(con, n = 100000)
    if (length(lines) == 0) break
    mrna_lines <- grep("\tmRNA\t", lines, value = TRUE)
    for (l in mrna_lines) {
      parts <- strsplit(l, "\t")[[1]]
      attrs <- parts[9]
      id_m <- regmatches(attrs, regexec("ID=([^;]+)", attrs))[[1]]
      parent_m <- regmatches(attrs, regexec("Parent=([^;]+)", attrs))[[1]]
      if (length(id_m) > 1 && length(parent_m) > 1) {
        mrna_to_gene[[id_m[2]]] <- parent_m[2]
      }
    }
  }
  close(con)
  cat(sprintf("  Built %d mRNA→gene mappings\n", length(mrna_to_gene)))

  # Second pass: extract CDS whose parent mRNA maps to a novel gene
  cat("  Pass 2: Extracting CDS for novel genes...\n")
  con <- file(gff3_path, "r")
  cds_list <- list()
  cds_count <- 0
  while (TRUE) {
    lines <- readLines(con, n = 100000)
    if (length(lines) == 0) break
    cds_lines <- grep("\tCDS\t", lines, value = TRUE)
    if (length(cds_lines) == 0) next

    for (l in cds_lines) {
      parts <- strsplit(l, "\t")[[1]]
      attrs <- parts[9]
      parent_m <- regmatches(attrs, regexec("Parent=([^;]+)", attrs))[[1]]
      if (length(parent_m) <= 1) next

      mrna_id <- parent_m[2]
      gene_id <- mrna_to_gene[[mrna_id]]
      if (is.null(gene_id)) next
      if (!gene_id %in% novel_ids) next

      cds_count <- cds_count + 1
      cds_list[[cds_count]] <- GRanges(
        seqnames = parts[1],
        ranges   = IRanges(start = as.integer(parts[4]), end = as.integer(parts[5])),
        strand   = parts[7],
        gene_id  = gene_id
      )
    }
  }
  close(con)
  cat(sprintf("  Extracted %d CDS fragments for novel genes\n", cds_count))

  if (cds_count == 0) {
    cat("  WARNING: No CDS found for novel genes!\n")
    return(NULL)
  }

  novel_cds_gr <- do.call(c, cds_list)

  # ---- Step 2: Extract sequences from genome ----
  cat("Loading genome FASTA...\n")
  genome <- readDNAStringSet(genome_path)
  names(genome) <- sub(" .*", "", names(genome))
  cat(sprintf("  %d contigs\n", length(genome)))

  cat("Extracting CDS sequences (batch mode)...\n")
  # Match contig names
  valid_idx <- which(as.character(seqnames(novel_cds_gr)) %in% names(genome))
  cat(sprintf("  %d/%d CDS have matching contig names\n", length(valid_idx), length(novel_cds_gr)))
  novel_cds_gr <- novel_cds_gr[valid_idx]

  # Get sequences in batch
  cds_seqs <- getSeq(genome, novel_cds_gr)

  # Reverse complement for negative strand
  neg_idx <- which(as.character(strand(novel_cds_gr)) == "-")
  if (length(neg_idx) > 0) {
    cds_seqs[neg_idx] <- reverseComplement(cds_seqs[neg_idx])
  }

  # ---- Step 3: Group by gene, concatenate, translate ----
  cat("Grouping CDS by gene and translating...\n")

  gene_ids <- mcols(novel_cds_gr)$gene_id
  unique_genes <- unique(gene_ids)
  cat(sprintf("  %d unique novel genes with CDS\n", length(unique_genes)))

  protein_list <- list()
  skipped <- 0

  for (i in seq_along(unique_genes)) {
    gid <- unique_genes[i]
    idx <- which(gene_ids == gid)

    g_seqs <- cds_seqs[idx]

    # Concatenate all CDS fragments
    concat <- DNAString(paste(as.character(g_seqs), collapse = ""))

    if (nchar(concat) < 30) {
      skipped <- skipped + 1
      next
    }

    # Check if length is multiple of 3, trim if needed
    remainder <- nchar(concat) %% 3
    if (remainder != 0) {
      concat <- subseq(concat, 1, nchar(concat) - remainder)
    }
    if (nchar(concat) < 30) {
      skipped <- skipped + 1
      next
    }

    # Translate (try multiple frames, pick the one with fewest stops)
    proteins <- list()
    best_protein <- NULL
    best_len <- 0

    for (frame in 0:2) {
      if (nchar(concat) <= frame + 3) next
      s <- subseq(concat, frame + 1, nchar(concat))
      p <- suppressWarnings(translate(s, if.fuzzy.codon = "solve"))
      p_char <- as.character(p)
      # Count stop codons (*) — fewer is better
      n_stops <- nchar(p_char) - nchar(gsub("\\*", "", p_char))
      p_len <- nchar(p_char) - n_stops  # length without stops
      if (n_stops == 0 && p_len > best_len) {
        best_protein <- p
        best_len <- p_len
      }
    }

    # If no perfect frame found, just use frame 0
    if (is.null(best_protein)) {
      best_protein <- proteins[[1]]
      if (is.null(best_protein)) {
        best_protein <- suppressWarnings(translate(
          subseq(concat, 1, nchar(concat)), if.fuzzy.codon = "solve"))
      }
    }

    gn <- novel_genes$gene_name[novel_genes$gene_id == gid]
    if (length(gn) == 0) gn <- gid
    gn <- gn[1]

    protein_list[[gn]] <- list(
      gene_id    = gid,
      novel_name = gn,
      protein    = as.character(best_protein),
      aa_len     = width(best_protein),
      cds_len    = nchar(concat),
      n_cds      = length(idx)
    )

    if (i %% 200 == 0) cat(sprintf("  Progress: %d/%d genes processed\n", i, length(unique_genes)))
  }

  cat(sprintf("  Translated: %d proteins, skipped: %d short CDS\n",
              length(protein_list), skipped))

  # ---- Step 4: Write protein FASTA ----
  fasta_out <- paste0("New-analysis/blast/", sp, "_novel_proteins.faa")
  con_out <- file(fasta_out, "w")
  for (gp in protein_list) {
    header <- paste0(">", gp$gene_id, "|novel_name=", gp$novel_name,
                     "|aa_len=", gp$aa_len, "|cds_len=", gp$cds_len)
    writeLines(header, con_out)
    writeLines(gp$protein, con_out)
  }
  close(con_out)
  cat(sprintf("  Wrote: %s (%d proteins)\n", fasta_out, length(protein_list)))

  return(protein_list)
}

# ---- Run for both species ----
mfu_proteins <- process_species("mfu", mfu_novel_genes)
mpi_proteins <- process_species("mpi", mpi_novel_genes)

cat("\n=== CDS extraction and translation complete ===\n")
cat("Output files:\n")
cat("  New-analysis/blast/mfu_novel_proteins.faa\n")
cat("  New-analysis/blast/mpi_novel_proteins.faa\n")
