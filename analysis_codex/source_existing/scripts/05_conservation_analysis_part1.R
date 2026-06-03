# ============================================================
# Step 1-2: Conservation Analysis after Annotation Rescue
# ============================================================
# 1. Update mfu/mpi gene_name with BLAST-rescued annotations
# 2. Build circ_id → gene_name lookup for all 6 species
# 3. Load CIRCexplorer3 expression + filter
# 4. Gene-level conservation analysis
# ============================================================

library(dplyr)
library(tidyr)

source(if (file.exists("New-analysis/00_project_setup.R")) "New-analysis/00_project_setup.R" else "00_project_setup.R", encoding = "UTF-8")
dir.create("New-analysis", showWarnings = FALSE)

cat("========== Step 1: Load circ_id → gene_name mappings ==========\n")

# ---- 1a: For hsa/mma/mmu/rsi — direct from *.cq.csv ----
load_cq_mapping <- function(cq_file, has_gene_name = TRUE) {
  cq <- read.csv(cq_file, stringsAsFactors = FALSE)
  if (has_gene_name) {
    # Has gene_name column
    mapping <- cq %>%
      select(circ_id, gene_name) %>%
      separate_rows(gene_name, sep = ",") %>%
      filter(!is.na(gene_name) & gene_name != "") %>%
      distinct()
    cat(sprintf("  %s: %d circ_ids → %d unique gene_names\n",
                cq_file, length(unique(mapping$circ_id)),
                length(unique(mapping$gene_name))))
  } else {
    # No gene_name — need GFF3 lookup (mfu/mpi before rescue)
    mapping <- cq %>%
      select(circ_id, gene_id) %>%
      separate_rows(gene_id, sep = ",") %>%
      distinct()
    cat(sprintf("  %s: %d circ_ids → %d unique gene_ids (no gene_name yet)\n",
                cq_file, length(unique(mapping$circ_id)),
                length(unique(mapping$gene_id))))
  }
  return(mapping)
}

# hsa, mma, mmu, rsi have gene_name
lookup_hsa <- load_cq_mapping("hsa.cq.csv")
lookup_mma <- load_cq_mapping("mma.cq.csv")
lookup_mmu <- load_cq_mapping("mmu.cq.csv")
lookup_rsi <- load_cq_mapping("rsi.cq.csv")

# ---- 1b: For mfu/mpi — use GFF3 + BLAST rescue ----
cat("\n--- Building mfu/mpi gene_name mappings (GFF3 + BLAST rescue) ---\n")

build_bat_mapping <- function(sp, cq_file, gff3_file, blast_mapping_file) {
  # Load CQ data (has gene_id but no gene_name)
  cq <- read.csv(cq_file, stringsAsFactors = FALSE)

  # Build gene_id → gene_name from GFF3
  cat(sprintf("  %s: Reading GFF3 for gene_id→gene_name...\n", sp))
  con <- file(gff3_file, "r")
  gff_names <- list()
  while (TRUE) {
    lines <- readLines(con, n = 50000)
    if (length(lines) == 0) break
    gene_lines <- grep("\tgene\t", lines, value = TRUE)
    for (l in gene_lines) {
      parts <- strsplit(l, "\t")[[1]]
      attrs <- parts[9]
      id_m <- regmatches(attrs, regexec("ID=([^;]+)", attrs))[[1]]
      name_m <- regmatches(attrs, regexec("Name=([^;]+)", attrs))[[1]]
      if (length(id_m) > 1 && length(name_m) > 1) {
        gff_names[[id_m[2]]] <- name_m[2]
      }
    }
  }
  close(con)
  cat(sprintf("  %d gene_id→Name mappings from GFF3\n", length(gff_names)))

  # Load BLAST rescue mappings if available
  if (file.exists(blast_mapping_file)) {
    blast_map <- read.csv(blast_mapping_file, stringsAsFactors = FALSE)
    # Override GFF3 names with BLAST names for novel genes
    cat(sprintf("  %d names rescued by BLAST\n", nrow(blast_map)))
    for (i in seq_len(nrow(blast_map))) {
      gid <- blast_map$gene_id[i]
      new_name <- blast_map$new_gene_name[i]
      if (!is.na(new_name) && new_name != "" && new_name != "UNKNOWN") {
        gff_names[[gid]] <- new_name
      }
    }
  }

  # Map circ_id → gene_name via gene_id
  mapping <- cq %>%
    select(circ_id, gene_id) %>%
    separate_rows(gene_id, sep = ",") %>%
    mutate(gene_id_clean = sub("\\.t\\d+$", "", gene_id)) %>%
    rowwise() %>%
    mutate(gene_name = ifelse(is.null(gff_names[[gene_id_clean]]),
                              NA, gff_names[[gene_id_clean]])) %>%
    ungroup() %>%
    filter(!is.na(gene_name)) %>%
    select(circ_id, gene_name) %>%
    distinct()

  cat(sprintf("  %s: %d circ_ids → %d unique gene_names\n",
              sp, length(unique(mapping$circ_id)),
              length(unique(mapping$gene_name))))
  return(mapping)
}

lookup_mfu <- build_bat_mapping("mfu", "mfu.cq.csv",
  "genome and annotation/mfu.gff3",
  "New-analysis/blast/mfu_novel_gene_mapping.csv")

lookup_mpi <- build_bat_mapping("mpi", "mpi.cq.csv",
  "genome and annotation/mpi.gff3",
  "New-analysis/blast/mpi_novel_gene_mapping.csv")

# ---- Combine all & normalize case ----
cat("\n========== Mapping summary (case-normalized) ==========\n")
gene_lookups <- list(hsa = lookup_hsa, mfu = lookup_mfu, mma = lookup_mma,
                     mmu = lookup_mmu, mpi = lookup_mpi, rsi = lookup_rsi)
for (sp in names(gene_lookups)) {
  n_before <- length(unique(gene_lookups[[sp]]$gene_name))
  gene_lookups[[sp]]$gene_name <- toupper(gene_lookups[[sp]]$gene_name)
  n_after <- length(unique(gene_lookups[[sp]]$gene_name))
  cat(sprintf("%s: %d circ_ids -> %d gene_names (merged %d by toupper)\n",
              sp, nrow(gene_lookups[[sp]]), n_after, n_before - n_after))
}

saveRDS(gene_lookups, "New-analysis/gene_lookups.rds")
cat("\nSaved: New-analysis/gene_lookups.rds\n")
