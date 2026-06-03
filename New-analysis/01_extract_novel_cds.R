# ============================================================
# Step 0: Extract CDS from novel_genes & BLAST annotation rescue
# ============================================================
# mfu/mpi bat species have many genes annotated as "novel_gene_*"
# Goal: use protein homology to assign real gene names
# ============================================================

library(Biostrings)
library(GenomicFeatures)
library(rtracklayer)
library(dplyr)
library(tidyr)

source(if (file.exists("New-analysis/00_project_setup.R")) "New-analysis/00_project_setup.R" else "00_project_setup.R", encoding = "UTF-8")
dir.create("New-analysis", showWarnings = FALSE)
dir.create("New-analysis/blast", showWarnings = FALSE)

# ---- Phase 1: Extract novel_gene → gene_id from GFF3 (text parsing, fast) ----

cat("========== Phase 1: Parse GFF3 for novel_gene entries ==========\n")

for (sp in c("mfu", "mpi")) {
  cat("\n---", sp, "---\n")

  gff3_path <- paste0("genome and annotation/", sp, ".gff3")

  # Parse GFF3 line by line, keep only gene features
  cat("  Reading GFF3 (gene lines only)...\n")
  con <- file(gff3_path, "r")
  gene_lines <- c()
  while (TRUE) {
    lines <- readLines(con, n = 50000)
    if (length(lines) == 0) break
    gene_lines <- c(gene_lines, grep("\tgene\t", lines, value = TRUE))
  }
  close(con)
  cat(sprintf("  Found %d gene features\n", length(gene_lines)))

  # Extract ID and Name from attribute column (col 9)
  extract_attr <- function(line) {
    attrs <- strsplit(line, "\t")[[1]][9]
    id_match <- regmatches(attrs, regexec("ID=([^;]+)", attrs))[[1]]
    name_match <- regmatches(attrs, regexec("Name=([^;]+)", attrs))[[1]]
    data.frame(
      gene_id = if(length(id_match) > 1) id_match[2] else NA,
      gene_name = if(length(name_match) > 1) name_match[2] else NA,
      stringsAsFactors = FALSE
    )
  }

  gene_info <- do.call(rbind, lapply(gene_lines, extract_attr))
  gene_info <- gene_info[!is.na(gene_info$gene_id), ]

  novel_genes <- gene_info[grepl("^novel_gene_", gene_info$gene_name), ]
  named_genes <- gene_info[!grepl("^novel_gene_", gene_info$gene_name) &
                           !is.na(gene_info$gene_name) &
                           gene_info$gene_name != "", ]

  cat(sprintf("  Named genes: %d\n", nrow(named_genes)))
  cat(sprintf("  Novel genes: %d (%.1f%%)\n", nrow(novel_genes),
              100 * nrow(novel_genes) / nrow(gene_info)))

  assign(paste0(sp, "_novel_genes"), novel_genes)
  assign(paste0(sp, "_named_genes"), named_genes)

  write.csv(novel_genes, paste0("New-analysis/", sp, "_novel_gene_list.csv"), row.names = FALSE)
  write.csv(named_genes, paste0("New-analysis/", sp, "_named_gene_list.csv"), row.names = FALSE)
}

cat("\n=== Phase 1 complete ===\n\n")

# ---- Phase 2: Extract CDS coordinates for novel genes from GFF3 ----

cat("========== Phase 2: Extract CDS for novel genes ==========\n")

for (sp in c("mfu", "mpi")) {
  cat("\n---", sp, "---\n")

  gff3_path <- paste0("genome and annotation/", sp, ".gff3")
  novel_genes <- get(paste0(sp, "_novel_genes"))

  # Parse GFF3 for mRNA and CDS lines
  cat("  Parsing mRNA and CDS features...\n")
  con <- file(gff3_path, "r")
  mrna_lines <- c()
  cds_lines <- c()
  while (TRUE) {
    lines <- readLines(con, n = 50000)
    if (length(lines) == 0) break
    mrna_lines <- c(mrna_lines, grep("\tmRNA\t", lines, value = TRUE))
    cds_lines <- c(cds_lines, grep("\tCDS\t", lines, value = TRUE))
  }
  close(con)
  cat(sprintf("  %d mRNA, %d CDS features\n", length(mrna_lines), length(cds_lines)))

  # Build mRNA_id -> gene_id mapping
  extract_mrna_attrs <- function(line) {
    parts <- strsplit(line, "\t")[[1]]
    attrs <- parts[9]
    id_match <- regmatches(attrs, regexec("ID=([^;]+)", attrs))[[1]]
    parent_match <- regmatches(attrs, regexec("Parent=([^;]+)", attrs))[[1]]
    data.frame(
      mrna_id = if(length(id_match) > 1) id_match[2] else NA,
      gene_id = if(length(parent_match) > 1) parent_match[2] else NA,
      stringsAsFactors = FALSE
    )
  }

  mrna_map <- do.call(rbind, lapply(mrna_lines, extract_mrna_attrs))
  mrna_map <- mrna_map[!is.na(mrna_map$mrna_id) & !is.na(mrna_map$gene_id), ]
  cat(sprintf("  %d mRNA-gene mappings\n", nrow(mrna_map)))

  # Novel gene IDs
  novel_ids <- novel_genes$gene_id
  novel_mrnas <- mrna_map[mrna_map$gene_id %in% novel_ids, ]
  cat(sprintf("  %d mRNAs belong to novel genes\n", nrow(novel_mrnas)))

  # Extract CDS features linked to novel gene mRNAs
  extract_cds_attrs <- function(line) {
    parts <- strsplit(line, "\t")[[1]]
    chrom   <- parts[1]
    start   <- as.integer(parts[4])
    end     <- as.integer(parts[5])
    strand  <- parts[7]
    attrs   <- parts[9]
    parent_match <- regmatches(attrs, regexec("Parent=([^;]+)", attrs))[[1]]
    data.frame(
      chrom  = chrom,
      start  = start,
      end    = end,
      strand = strand,
      parent = if(length(parent_match) > 1) parent_match[2] else NA,
      stringsAsFactors = FALSE
    )
  }

  cds_df <- do.call(rbind, lapply(cds_lines, extract_cds_attrs))
  cds_df <- cds_df[!is.na(cds_df$parent), ]

  # Filter CDS to novel gene mRNAs
  novel_mrna_ids <- unique(novel_mrnas$mrna_id)
  novel_cds <- cds_df[cds_df$parent %in% novel_mrna_ids, ]
  cat(sprintf("  %d CDS features from novel genes\n", nrow(novel_cds)))

  # Merge with gene_id
  novel_cds <- merge(novel_cds, novel_mrnas,
                     by.x = "parent", by.y = "mrna_id", all.x = TRUE)

  assign(paste0(sp, "_novel_cds"), novel_cds)
  assign(paste0(sp, "_mrna_map"), mrna_map)

  write.csv(novel_cds, paste0("New-analysis/", sp, "_novel_cds_coords.csv"), row.names = FALSE)
}

cat("\n=== Phase 2 complete ===\n\n")

# ---- Phase 3: Extract CDS sequences from genome FASTA ----

cat("========== Phase 3: Extract CDS sequences ==========\n")

for (sp in c("mfu", "mpi")) {
  cat("\n---", sp, "---\n")

  genome_path <- paste0("genome and annotation/", sp, ".fa")
  novel_cds <- get(paste0(sp, "_novel_cds"))
  novel_genes <- get(paste0(sp, "_novel_genes"))

  # Load genome as DNAStringSet
  cat("  Loading genome FASTA...\n")
  genome <- readDNAStringSet(genome_path)

  # Build contig name index (FASTA headers are often just the contig name after >)
  contig_names <- names(genome)
  # Simplify: remove anything after first space
  names(genome) <- sub(" .*", "", contig_names)
  cat(sprintf("  %d contigs loaded\n", length(genome)))

  # Group CDS by gene_id
  novel_cds$gene_id <- ifelse(is.na(novel_cds$gene_id), novel_cds$parent, novel_cds$gene_id)

  # Get unique novel gene IDs
  gene_ids <- unique(novel_cds$gene_id)
  cat(sprintf("  Extracting CDS for %d novel genes...\n", length(gene_ids)))

  # For each novel gene, concatenate CDS and translate
  gene_proteins <- list()

  for (i in seq_along(gene_ids)) {
    gid <- gene_ids[i]
    g_cds <- novel_cds[novel_cds$gene_id == gid, ]

    if (nrow(g_cds) == 0) next

    # Get the contig(s) these CDS are on
    contigs_used <- unique(g_cds$chrom)

    # Extract each CDS fragment
    cds_seqs <- DNAStringSet()
    for (j in seq_len(nrow(g_cds))) {
      row <- g_cds[j, ]
      contig <- row$chrom
      if (!contig %in% names(genome)) {
        # Try alternative contig name
        alt_contig <- grep(contig, names(genome), value = TRUE, fixed = TRUE)
        if (length(alt_contig) == 0) next
        contig <- alt_contig[1]
      }
      tryCatch({
        seq_frag <- subseq(genome[[contig]], start = row$start, end = row$end)
        if (row$strand == "-") {
          seq_frag <- reverseComplement(seq_frag)
        }
        cds_seqs <- c(cds_seqs, seq_frag)
      }, error = function(e) NULL)
    }

    if (length(cds_seqs) == 0) next

    # Concatenate all CDS fragments for this gene
    full_cds <- unlist(cds_seqs)
    if (length(full_cds) == 0) next

    full_cds_str <- as.character(full_cds)
    if (nchar(full_cds_str) < 30) next  # Skip very short CDS

    concatenated <- DNAString(paste(full_cds_str, collapse = ""))

    # Translate to protein (using standard genetic code)
    protein <- suppressWarnings(translate(concatenated, if.fuzzy.codon = "solve"))

    # Get the novel_gene name
    gn <- novel_genes$gene_name[novel_genes$gene_id == gid]
    if (length(gn) == 0) gn <- gid
    gn <- gn[1]

    gene_proteins[[gn]] <- list(
      gene_id = gid,
      novel_name = gn,
      protein = as.character(protein),
      aa_length = width(protein),
      cds_length = nchar(full_cds_str),
      n_cds_frags = nrow(g_cds)
    )
  }

  cat(sprintf("  Successfully translated %d proteins\n", length(gene_proteins)))

  # Save as FASTA
  fasta_out <- paste0("New-analysis/blast/", sp, "_novel_proteins.faa")
  con_out <- file(fasta_out, "w")
  for (i in seq_along(gene_proteins)) {
    gp <- gene_proteins[[i]]
    header <- paste0(">", gp$novel_name, "|", gp$gene_id, "|len=", gp$aa_length)
    writeLines(header, con_out)
    writeLines(gp$protein, con_out)
  }
  close(con_out)
  cat(sprintf("  Wrote: %s\n", fasta_out))

  assign(paste0(sp, "_novel_proteins"), gene_proteins)
}

cat("\n=== Phase 3 complete ===\n\n")
