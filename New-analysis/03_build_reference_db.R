# ============================================================
# Step 0 Phase 4: Build reference protein database from hsa/mma/mmu
# ============================================================
# Use Ensembl GTFs to extract protein-coding gene CDS,
# translate, and create a BLAST database for annotation rescue.
# ============================================================

library(Biostrings)
library(GenomicRanges)
library(rtracklayer)
library(dplyr)

source(if (file.exists("New-analysis/00_project_setup.R")) "New-analysis/00_project_setup.R" else "00_project_setup.R", encoding = "UTF-8")
dir.create("New-analysis/blast", showWarnings = FALSE, recursive = TRUE)

# ---- Species to process ----
MODEL_SPECIES <- c("hsa", "mma", "mmu")

for (sp in MODEL_SPECIES) {
  cat("\n========== Processing", sp, "==========\n")

  gtf_path   <- paste0("genome and annotation/", sp, ".gtf")
  genome_path <- paste0("genome and annotation/", sp, ".fa")

  # ---- Step 1: Parse GTF for CDS features ----
  cat("  Parsing GTF for CDS features...\n")

  con <- file(gtf_path, "r")
  cds_records <- list()
  count <- 0

  # Note: Ensembl GTF has gene_id, gene_name, transcript_id in column 9
  hsa_has_gene_name <- TRUE  # Ensembl GTF with gene_name

  while (TRUE) {
    lines <- readLines(con, n = 100000)
    if (length(lines) == 0) break

    # Skip comments
    lines <- lines[!grepl("^#", lines)]
    if (length(lines) == 0) next

    # Keep only CDS lines
    cds_lines <- grep("\tCDS\t", lines, value = TRUE)
    if (length(cds_lines) == 0) next

    for (l in cds_lines) {
      parts <- strsplit(l, "\t")[[1]]
      if (length(parts) < 9) next

      attrs <- parts[9]

      # Extract gene_id, gene_name, transcript_id
      gid_m  <- regmatches(attrs, regexec('gene_id "([^"]+)"', attrs))[[1]]
      gname_m <- regmatches(attrs, regexec('gene_name "([^"]+)"', attrs))[[1]]
      tid_m  <- regmatches(attrs, regexec('transcript_id "([^"]+)"', attrs))[[1]]

      gene_id   <- if(length(gid_m) > 1) gid_m[2] else NA
      gene_name <- if(length(gname_m) > 1) gname_m[2] else NA
      tx_id     <- if(length(tid_m) > 1) tid_m[2] else NA

      count <- count + 1
      cds_records[[count]] <- data.frame(
        seqnames = parts[1],
        start    = as.integer(parts[4]),
        end      = as.integer(parts[5]),
        strand   = parts[7],
        gene_id  = gene_id,
        gene_name = gene_name,
        tx_id    = tx_id,
        stringsAsFactors = FALSE
      )
    }
  }
  close(con)
  cat(sprintf("  Found %d CDS features\n", count))

  cds_df <- bind_rows(cds_records)

  # Filter for protein-coding (exclude pseudogenes, etc. — we just take all CDS)
  unique_genes <- unique(cds_df$gene_id)
  cat(sprintf("  %d unique genes with CDS\n", length(unique_genes)))

  # ---- Step 2: Create GRanges and extract sequences ----
  cat("  Loading genome FASTA...\n")
  genome <- readDNAStringSet(genome_path)
  names(genome) <- sub(" .*", "", names(genome))
  cat(sprintf("  %d contigs loaded\n", length(genome)))

  cds_gr <- makeGRangesFromDataFrame(cds_df, keep.extra.columns = TRUE)

  # Filter to valid contigs
  valid_idx <- which(as.character(seqnames(cds_gr)) %in% names(genome))
  cat(sprintf("  %d/%d CDS match genome contigs\n", length(valid_idx), nrow(cds_df)))
  cds_gr <- cds_gr[valid_idx]
  cds_df <- cds_df[valid_idx, ]

  # Extract in batch
  cat("  Extracting CDS sequences...\n")
  cds_seqs <- getSeq(genome, cds_gr)

  # Reverse complement negative strand
  neg_idx <- which(as.character(strand(cds_gr)) == "-")
  if (length(neg_idx) > 0) {
    cds_seqs[neg_idx] <- reverseComplement(cds_seqs[neg_idx])
  }

  # ---- Step 3: Group by gene, translate ----
  cat("  Grouping by gene and translating...\n")

  # Group by gene_id to get the longest protein per gene
  gene_groups <- split(seq_along(cds_seqs), cds_df$gene_id)
  unique_gids <- names(gene_groups)
  cat(sprintf("  %d gene groups\n", length(unique_gids)))

  protein_list <- list()

  for (i in seq_along(unique_gids)) {
    gid <- unique_gids[i]
    idx <- gene_groups[[gid]]

    g_seqs <- cds_seqs[idx]

    # Concatenate all CDS for this gene
    full_cds <- DNAString(paste(as.character(g_seqs), collapse = ""))

    if (nchar(full_cds) < 30) next

    # Try multiple frames, pick the best (fewest stop codons)
    best_aa <- ""
    best_len <- 0
    best_nstops <- Inf

    for (frame in 0:2) {
      if (nchar(full_cds) <= frame + 3) next
      s <- subseq(full_cds, frame + 1, nchar(full_cds))
      remainder <- nchar(s) %% 3
      if (remainder > 0) s <- subseq(s, 1, nchar(s) - remainder)
      if (nchar(s) < 30) next

      p <- suppressWarnings(translate(s, if.fuzzy.codon = "solve"))
      p_char <- as.character(p)
      n_stops <- nchar(p_char) - nchar(gsub("\\*", "", p_char))
      p_len <- nchar(p_char) - n_stops

      if (n_stops < best_nstops ||
          (n_stops == best_nstops && p_len > best_len)) {
        best_aa <- p_char
        best_len <- p_len
        best_nstops <- n_stops
      }
    }

    if (best_len < 10) next

    # Get gene_name (from the GTF, or fall back to gene_id)
    gn <- unique(cds_df$gene_name[cds_df$gene_id == gid])
    gn <- gn[!is.na(gn) & gn != ""]
    if (length(gn) == 0) gn <- gid
    gn <- gn[1]

    protein_list[[gid]] <- list(
      gene_id   = gid,
      gene_name = gn,
      protein   = best_aa,
      aa_len    = best_len
    )

    if (i %% 1000 == 0) cat(sprintf("  Progress: %d/%d genes\n", i, length(unique_gids)))
  }

  cat(sprintf("  Translated %d proteins\n", length(protein_list)))

  # ---- Step 4: Write FASTA ----
  fasta_out <- paste0("New-analysis/blast/", sp, "_reference_proteins.faa")
  con_out <- file(fasta_out, "w")
  for (gp in protein_list) {
    header <- paste0(">", sp, "|", gp$gene_name, "|", gp$gene_id, "|len=", gp$aa_len)
    writeLines(header, con_out)
    writeLines(gp$protein, con_out)
  }
  close(con_out)
  cat(sprintf("  Wrote: %s\n", fasta_out))
}

cat("\n=== Reference protein database complete ===\n")
cat("Files:\n")
for (sp in MODEL_SPECIES) {
  cat("  New-analysis/blast/", sp, "_reference_proteins.faa\n", sep = "")
}
