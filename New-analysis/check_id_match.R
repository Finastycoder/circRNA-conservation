source(if (file.exists("New-analysis/00_project_setup.R")) "New-analysis/00_project_setup.R" else "00_project_setup.R", encoding = "UTF-8")
novel <- read.csv("New-analysis/mfu_novel_gene_list.csv", stringsAsFactors = FALSE)
mapping <- read.csv("New-analysis/blast/mfu_novel_gene_mapping.csv", stringsAsFactors = FALSE)
mapping <- mapping[mapping$new_gene_name != "UNKNOWN" & !is.na(mapping$new_gene_name), ]

cat("=== GFF3 novel gene_ids (first 8) ===\n")
cat(head(novel$gene_id, 8), sep = "\n")

cat("\n=== BLAST mapping gene_ids (first 8) ===\n")
cat(head(unique(mapping$gene_id), 8), sep = "\n")

cat("\n=== Overlap ===\n")
cat("Mapping IDs in novel list:", sum(unique(mapping$gene_id) %in% novel$gene_id),
    "/", length(unique(mapping$gene_id)), "\n")
cat("Novel IDs in mapping:", sum(novel$gene_id %in% unique(mapping$gene_id)),
    "/", nrow(novel), "\n")

# Check a few IDs that are in mapping but NOT in novel
diff_ids <- setdiff(unique(mapping$gene_id), novel$gene_id)
cat("\n=== IDs in BLAST but NOT in novel list (first 10) ===\n")
cat(head(diff_ids, 10), sep = "\n")

# Check if these IDs exist in the GFF3 at all
cat("\n=== Checking if these are real gene IDs ===\n")
# Look at first diff ID in GFF3
test_id <- diff_ids[1]
cat("Searching GFF3 for:", test_id, "\n")
system(sprintf('findstr /c:"%s" "genome and annotation\\mfu.gff3"', test_id))
