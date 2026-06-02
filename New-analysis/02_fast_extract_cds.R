# ============================================================
# Fast CDS extraction using GTF (not GFF3!)
# GTF CDS lines directly have gene_id — no mRNA hierarchy needed
# ============================================================

library(Biostrings)
library(dplyr)

setwd("D:/跨物种筛选保守circRNA")

# ---- Load novel gene lists ----
mfu_novel <- read.csv("New-analysis/mfu_novel_gene_list.csv", stringsAsFactors = FALSE)
mpi_novel <- read.csv("New-analysis/mpi_novel_gene_list.csv", stringsAsFactors = FALSE)

# ---- Fast extraction function ----
extract_novel_proteins <- function(sp, novel_genes, gtf_path, genome_path) {
  cat("\n========== Processing", sp, "==========\n")

  novel_ids <- unique(novel_genes$gene_id)
  cat("Novel gene IDs:", length(novel_ids), "\n")

  # ---- Step 1: Extract CDS coordinates from GTF ----
  # GTF format: CDS lines have gene_id "XXX" in the attributes
  # This is much faster: we just read all CDS lines, extract gene_id,
  # and filter for novel genes in one pass

  cat("Reading GTF CDS lines...\n")
  con <- file(gtf_path, "r")
  cds_records <- list()
  count <- 0
  matched <- 0

  while (TRUE) {
    lines <- readLines(con, n = 200000)
    if (length(lines) == 0) break

    # Keep only CDS lines (skip comments)
    cds_lines <- lines[grepl("\tCDS\t", lines)]
    if (length(cds_lines) == 0) next

    for (l in cds_lines) {
      # Parse GTF line: tab-separated
      parts <- strsplit(l, "\t")[[1]]
      if (length(parts) < 9) next

      # Extract gene_id from attributes (column 9)
      attrs <- parts[9]
      gid_match <- regmatches(attrs, regexec('gene_id "([^"]+)"', attrs))[[1]]

      if (length(gid_match) <= 1) next
      gene_id <- gid_match[2]

      # Check if this gene_id is a novel gene
      if (!gene_id %in% novel_ids) next

      matched <- matched + 1
      count <- count + 1
      cds_records[[count]] <- list(
        seqnames = parts[1],
        start    = as.integer(parts[4]),
        end      = as.integer(parts[5]),
        strand   = parts[7],
        gene_id  = gene_id
      )
    }
  }
  close(con)

  cat(sprintf("Found %d CDS fragments for novel genes\n", matched))

  if (length(cds_records) == 0) {
    cat("WARNING: No CDS found!\n")
    return(NULL)
  }

  # Convert to data frame
  cds_df <- bind_rows(cds_records)
  cat(sprintf("Unique novel genes with CDS: %d\n",
              length(unique(cds_df$gene_id))))

  # ---- Step 2: Load genome and extract sequences ----
  cat("Loading genome FASTA...\n")
  genome <- readDNAStringSet(genome_path)
  names(genome) <- sub(" .*", "", names(genome))
  cat(sprintf("%d contigs\n", length(genome)))

  # Filter to contigs that exist in genome
  valid <- which(cds_df$seqnames %in% names(genome))
  cat(sprintf("%d/%d CDS match genome contigs\n", length(valid), nrow(cds_df)))
  if (length(valid) == 0) return(NULL)

  # Use DNAStringSet for extraction (faster than GRanges for DNAStringSet)
  cat("Extracting sequences by contig...\n")

  # Group CDS by contig for batch extraction
  cds_df <- cds_df[valid, ]
  contig_groups <- split(seq_len(nrow(cds_df)), cds_df$seqnames)
  cat(sprintf("  Processing %d contigs...\n", length(contig_groups)))

  # Pre-allocate list for CDS sequences
  cds_seqs_list <- vector("list", nrow(cds_df))
  len_counter <- 0

  for (cg_name in names(contig_groups)) {
    idx <- contig_groups[[cg_name]]
    contig_seq <- genome[[cg_name]]

    for (j in idx) {
      tryCatch({
        frag <- subseq(contig_seq, start = cds_df$start[j], end = cds_df$end[j])
        strand <- cds_df$strand[j]
        if (!is.na(strand) && strand == "-") {
          frag <- reverseComplement(frag)
        }
        cds_seqs_list[[j]] <- frag
        len_counter <- len_counter + 1
      }, error = function(e) {
        # Silently skip coordinates outside contig bounds
      })
    }
    if (len_counter %% 20000 == 0 && len_counter > 0) {
      cat(sprintf("  Extracted %d fragments...\n", len_counter))
    }
  }

  # Remove NULL entries and combine
  valid_seqs <- !sapply(cds_seqs_list, is.null)
  cat(sprintf("  Successfully extracted %d/%d CDS fragments\n", sum(valid_seqs), nrow(cds_df)))
  cds_seqs <- DNAStringSet(cds_seqs_list[valid_seqs])
  cds_df <- cds_df[valid_seqs, ]

  if (length(cds_seqs) == 0) return(NULL)

  # ---- Step 3: Group by gene, concatenate CDS, translate ----
  cat("Translating proteins...\n")

  gene_groups <- split(seq_len(nrow(cds_df)), cds_df$gene_id)
  gene_ids <- names(gene_groups)
  cat(sprintf("%d gene groups to translate\n", length(gene_ids)))

  proteins <- list()
  skipped <- 0

  for (i in seq_along(gene_ids)) {
    gid <- gene_ids[i]
    idx <- gene_groups[[gid]]

    # Concatenate CDS fragments
    concat <- DNAString(paste(as.character(cds_seqs[idx]), collapse = ""))
    if (nchar(concat) < 30) { skipped <- skipped + 1; next }

    # Trim to multiple of 3
    rem <- nchar(concat) %% 3
    if (rem > 0) concat <- subseq(concat, 1, nchar(concat) - rem)
    if (nchar(concat) < 30) { skipped <- skipped + 1; next }

    # Try 3 frames, pick best (fewest stop codons)
    best_aa <- ""
    best_len <- 0
    best_stops <- Inf

    for (frame in 0:2) {
      if (nchar(concat) <= frame + 3) next
      s <- subseq(concat, frame + 1, nchar(concat))
      p <- suppressWarnings(translate(s, if.fuzzy.codon = "solve"))
      p_char <- as.character(p)
      n_stops <- nchar(p_char) - nchar(gsub("\\*", "", p_char))
      p_len <- nchar(p_char) - n_stops

      if (n_stops < best_stops || (n_stops == best_stops && p_len > best_len)) {
        best_aa <- p_char
        best_len <- p_len
        best_stops <- n_stops
      }
    }

    if (best_len < 10) { skipped <- skipped + 1; next }

    gn <- novel_genes$gene_name[novel_genes$gene_id == gid]
    gn <- if(length(gn) > 0) gn[1] else gid

    proteins[[gid]] <- list(
      gene_id    = gid,
      novel_name = gn,
      protein    = best_aa,
      aa_len     = best_len,
      cds_len    = nchar(concat)
    )

    if (i %% 200 == 0) cat(sprintf("  %d/%d genes translated\n", i, length(gene_ids)))
  }

  cat(sprintf("Translated: %d proteins, skipped: %d\n", length(proteins), skipped))

  # ---- Step 4: Write FASTA ----
  fasta_out <- paste0("New-analysis/blast/", sp, "_novel_proteins.faa")
  con_out <- file(fasta_out, "w")
  for (gp in proteins) {
    header <- paste0(">", gp$gene_id, "|novel_name=", gp$novel_name,
                     "|aa_len=", gp$aa_len, "|cds_len=", gp$cds_len)
    writeLines(header, con_out)
    writeLines(gp$protein, con_out)
  }
  close(con_out)
  cat(sprintf("Wrote: %s (%d proteins)\n", fasta_out, length(proteins)))

  return(proteins)
}

# ---- Run ----
# mfu
mfu_prots <- extract_novel_proteins("mfu", mfu_novel,
  "genome and annotation/mfu.gtf",
  "genome and annotation/mfu.fa")

# mpi
mpi_prots <- extract_novel_proteins("mpi", mpi_novel,
  "genome and annotation/mpi.gtf",
  "genome and annotation/mpi.fa")

cat("\n=== All done! ===\n")
cat("Output:\n")
cat("  New-analysis/blast/mfu_novel_proteins.faa\n")
cat("  New-analysis/blast/mpi_novel_proteins.faa\n")
