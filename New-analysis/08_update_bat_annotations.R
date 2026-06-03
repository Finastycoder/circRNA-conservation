# ============================================================
# Update bat GFF3/GTF annotations with BLAST-rescued gene names
# Output to: bat_annotation/
# ============================================================

library(dplyr)

source(if (file.exists("New-analysis/00_project_setup.R")) "New-analysis/00_project_setup.R" else "00_project_setup.R", encoding = "UTF-8")
dir.create("bat_annotation", showWarnings = FALSE)

for (sp in c("mfu", "mpi")) {
  cat("\n========== Processing", sp, "==========\n")

  # ---- 1. Build gene_id → best_gene_name lookup ----
  # Start with GFF3 names
  cat("  Building gene_id → gene_name lookup...\n")

  gff3_path <- paste0("genome and annotation/", sp, ".gff3")
  con <- file(gff3_path, "r")
  gff_names <- list()
  while (TRUE) {
    lines <- readLines(con, n = 50000)
    if (length(lines) == 0) break
    gene_lines <- grep("\tgene\t", lines, value = TRUE)
    for (l in gene_lines) {
      parts <- strsplit(l, "\t")[[1]]
      attrs <- parts[9]
      id_m   <- regmatches(attrs, regexec("ID=([^;]+)", attrs))[[1]]
      name_m <- regmatches(attrs, regexec("Name=([^;]+)", attrs))[[1]]
      if (length(id_m) > 1 && length(name_m) > 1) {
        gff_names[[id_m[2]]] <- name_m[2]
      }
    }
  }
  close(con)
  cat(sprintf("  %d gene names from GFF3\n", length(gff_names)))

  # Override with BLAST-rescued names
  blast_file <- paste0("New-analysis/blast/", sp, "_novel_gene_mapping.csv")
  if (file.exists(blast_file)) {
    blast_map <- read.csv(blast_file, stringsAsFactors = FALSE)
    blast_map <- blast_map[!is.na(blast_map$new_gene_name) &
                           blast_map$new_gene_name != "UNKNOWN" &
                           blast_map$new_gene_name != "", ]

    # For each gene_id, take the best (highest identity) BLAST hit
    blast_best <- blast_map %>%
      group_by(gene_id) %>%
      slice_max(order_by = identity, n = 1) %>%
      ungroup()

    n_updated <- 0
    for (i in seq_len(nrow(blast_best))) {
      gid <- blast_best$gene_id[i]
      new_name <- blast_best$new_gene_name[i]
      old_name <- gff_names[[gid]]
      if (!is.null(old_name) && grepl("^novel_gene_", old_name)) {
        gff_names[[gid]] <- new_name
        n_updated <- n_updated + 1
      }
    }
    cat(sprintf("  %d novel_gene names replaced with BLAST results\n", n_updated))
  }

  # ---- 2. Update GFF3 file ----
  cat("  Updating GFF3...\n")

  gff_in  <- file(gff3_path, "r")
  gff_out <- file(paste0("bat_annotation/", sp, ".gff3"), "w")

  while (TRUE) {
    lines <- readLines(gff_in, n = 50000)
    if (length(lines) == 0) break

    for (l in lines) {
      if (grepl("^#", l)) {
        writeLines(l, gff_out)
        next
      }

      parts <- strsplit(l, "\t")[[1]]
      if (length(parts) < 9) {
        writeLines(l, gff_out)
        next
      }

      attrs <- parts[9]

      # Check if this line has a Name that is novel_gene and we have a replacement
      name_m <- regmatches(attrs, regexec("Name=([^;]+)", attrs))[[1]]
      if (length(name_m) > 1) {
        old_name <- name_m[2]
        # Get the gene ID from this line or its parent
        id_m <- regmatches(attrs, regexec("ID=([^;]+)", attrs))[[1]]
        if (length(id_m) > 1) {
          gid <- id_m[2]
        } else {
          parent_m <- regmatches(attrs, regexec("Parent=([^;]+)", attrs))[[1]]
          gid <- if (length(parent_m) > 1) parent_m[2] else NULL
        }

        if (!is.null(gid)) {
          new_name <- gff_names[[gid]]
          if (!is.null(new_name) && new_name != old_name) {
            attrs <- gsub(paste0("Name=", old_name),
                         paste0("Name=", new_name), attrs, fixed = TRUE)
            parts[9] <- attrs
            l <- paste(parts, collapse = "\t")
          }
        }
      }

      writeLines(l, gff_out)
    }
  }
  close(gff_in)
  close(gff_out)
  cat(sprintf("  Wrote: bat_annotation/%s.gff3\n", sp))

  # ---- 3. Update GTF file (add gene_name attribute) ----
  cat("  Updating GTF (adding gene_name)...\n")

  gtf_path <- paste0("genome and annotation/", sp, ".gtf")
  gtf_in  <- file(gtf_path, "r")
  gtf_out <- file(paste0("bat_annotation/", sp, ".gtf"), "w")

  while (TRUE) {
    lines <- readLines(gtf_in, n = 50000)
    if (length(lines) == 0) break

    for (l in lines) {
      if (grepl("^#", l)) {
        writeLines(l, gtf_out)
        next
      }

      parts <- strsplit(l, "\t")[[1]]
      if (length(parts) < 9) {
        writeLines(l, gtf_out)
        next
      }

      attrs <- parts[9]

      # Check if gene_name already exists
      if (!grepl("gene_name", attrs)) {
        # Extract gene_id
        gid_m <- regmatches(attrs, regexec('gene_id "([^"]+)"', attrs))[[1]]
        if (length(gid_m) > 1) {
          gid <- gid_m[2]

          # Also try without transcript suffix (.t1, .t2)
          gid_base <- sub("\\.t\\d+$", "", gid)

          gn <- gff_names[[gid]]
          if (is.null(gn)) gn <- gff_names[[gid_base]]

          if (!is.null(gn)) {
            # Add gene_name to attributes (before the final semicolon or at the end)
            if (grepl(";$", attrs)) {
              attrs <- paste0(attrs, ' gene_name "', gn, '";')
            } else {
              attrs <- paste0(attrs, '; gene_name "', gn, '";')
            }
            parts[9] <- attrs
            l <- paste(parts, collapse = "\t")
          }
        }
      }

      writeLines(l, gtf_out)
    }
  }
  close(gtf_in)
  close(gtf_out)
  cat(sprintf("  Wrote: bat_annotation/%s.gtf\n", sp))
}

# ---- 4. Copy rsi files (already complete) ----
cat("\n--- Copying rsi (already complete) ---\n")
file.copy("genome and annotation/rsi_2.gtf", "bat_annotation/rsi.gtf", overwrite = TRUE)
file.copy("genome and annotation/rsi_2.gff", "bat_annotation/rsi.gff3", overwrite = TRUE)
cat("  Copied: bat_annotation/rsi.gtf\n")
cat("  Copied: bat_annotation/rsi.gff3\n")

# ---- 5. Summary ----
cat("\n========== bat_annotation/ contents ==========\n")
files <- list.files("bat_annotation", full.names = TRUE)
for (f in files) {
  size_mb <- file.info(f)$size / 1e6
  cat(sprintf("  %s (%.1f MB)\n", basename(f), size_mb))
}
cat("\nDone!\n")
