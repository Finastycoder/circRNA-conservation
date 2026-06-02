setwd("D:/跨物种筛选保守circRNA")

for (sp in c("mfu", "mpi")) {
  cat("\n========== ", sp, " ==========\n", sep = "")

  novel <- read.csv(paste0("New-analysis/", sp, "_novel_gene_list.csv"), stringsAsFactors = FALSE)
  named <- read.csv(paste0("New-analysis/", sp, "_named_gene_list.csv"), stringsAsFactors = FALSE)
  mapping <- read.csv(paste0("New-analysis/blast/", sp, "_novel_gene_mapping.csv"), stringsAsFactors = FALSE)

  # Filter out UNKNOWN
  mapping <- mapping[!is.na(mapping$new_gene_name) & mapping$new_gene_name != "UNKNOWN", ]

  n_named   <- nrow(named)
  n_novel   <- nrow(novel)
  tot       <- n_named + n_novel
  n_mapping_rows <- nrow(mapping)
  n_unique_mapped <- length(unique(mapping$gene_id))
  rescued   <- length(intersect(unique(mapping$gene_id), novel$gene_id))

  cat("Total genes: ", tot, "\n")
  cat("Originally named: ", n_named, " (", round(100*n_named/tot, 1), "%)\n", sep = "")
  cat("Originally novel: ", n_novel, " (", round(100*n_novel/tot, 1), "%)\n", sep = "")
  cat("BLAST mapping rows: ", n_mapping_rows, "\n", sep = "")
  cat("Unique genes in BLAST mapping: ", n_unique_mapped, "\n", sep = "")
  cat("Rescued (overlap with novel list): ", rescued, " (", round(100*rescued/n_novel, 1), "% of novel)\n", sep = "")
  cat("Still novel: ", n_novel - rescued, "\n", sep = "")
  cat("FINAL named: ", n_named + rescued, " (", round(100*(n_named+rescued)/tot, 1), "%)\n", sep = "")
}
