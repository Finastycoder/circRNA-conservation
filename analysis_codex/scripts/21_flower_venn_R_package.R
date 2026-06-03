# ============================================================
# Six-species flower Venn plots using the R package `venn`
# ============================================================
# Recommended R runtime:
#   S:/R-4.2.2/bin/Rscript.exe analysis_codex/scripts/21_flower_venn_R_package.R
#
# Style source:
#   analysis_codex/multi_milk.R lines around the multispecies host-gene Venn:
#   library(venn); venn(venn_list, zcolor = "style", ilabels = FALSE,
#   ellipse = FALSE, opacity = 0.5, box = FALSE, borders = FALSE,
#   ilcs = 0.8, sncs = 1.5); get.venn.partitions(venn_list).
#
# Why not VennDiagram?
#   VennDiagram is available in R 4.2.2, but it exports up to
#   draw.quintuple.venn only. It cannot draw a 6-set Venn directly.
#   The `venn` package supports 6 sets with ellipse=TRUE and gives a
#   flower-like plot closer to 花瓣图示例.png.
# ============================================================

suppressPackageStartupMessages({
  library(venn)
  library(VennDiagram)
})

find_project_root <- function() {
  script_file <- tryCatch(sys.frame(1)$ofile, error = function(e) NA_character_)
  if (is.null(script_file) || !length(script_file) || is.na(script_file)) {
    file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
    script_file <- if (length(file_arg)) sub("^--file=", "", file_arg[1]) else NA_character_
  }
  script_dir <- if (!is.na(script_file) && nzchar(script_file)) {
    dirname(normalizePath(script_file, winslash = "/", mustWork = FALSE))
  } else {
    NA_character_
  }
  candidates <- c(
    getwd(),
    dirname(getwd()),
    if (!is.na(script_dir)) file.path(script_dir, "..", "..") else NA_character_
  )
  candidates <- unique(normalizePath(candidates[!is.na(candidates)], winslash = "/", mustWork = FALSE))
  for (candidate in candidates) {
    if (dir.exists(file.path(candidate, "analysis_codex")) &&
        file.exists(file.path(candidate, "New-analysis", "gene_presence_matrix.csv"))) {
      return(candidate)
    }
  }
  stop("Cannot find project root. Run from project root or analysis_codex/scripts.")
}

ROOT <- find_project_root()
FIG_DIR <- file.path(ROOT, "analysis_codex", "figures")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

DATA_SPECIES <- c("hsa", "mma", "mmu", "mfu", "mpi", "rsi")
PLOT_SPECIES <- c("hsa", "mfu", "mma", "mmu", "mpi", "rsi")
PLOT_NAMES <- c(hsa = "HSA", mfu = "MFU", mma = "MMA", mmu = "MMU", mpi = "MPI", rsi = "RSI")
PLOT_COLORS <- c(
  hsa = "#F8766D",
  mfu = "#FDAE61",
  mma = "#56B4E9",
  mmu = "#66C2A5",
  mpi = "#1F78B4",
  rsi = "#7B61FF"
)

parse_bool <- function(x) {
  toupper(as.character(x)) %in% c("TRUE", "T", "YES", "1")
}

strict_milk_host_gene_sets <- function() {
  presence <- read.csv(file.path(ROOT, "New-analysis", "gene_presence_matrix.csv"),
                       stringsAsFactors = FALSE, check.names = FALSE)
  presence$gene_name_upper <- toupper(presence$gene_name)
  sets <- lapply(PLOT_SPECIES, function(sp) {
    unique(presence$gene_name_upper[parse_bool(presence[[sp]])])
  })
  names(sets) <- PLOT_NAMES[PLOT_SPECIES]
  sets
}

extract_circatlas_host_genes <- function(path) {
  x <- read.delim(path, stringsAsFactors = FALSE, check.names = FALSE)
  id_col <- grep("circ", names(x), ignore.case = TRUE, value = TRUE)[1]
  ids <- x[[id_col]]
  genes <- sub("^[^-]+-", "", ids)
  genes <- sub("_[0-9]+$", "", genes)
  genes <- toupper(genes)
  unique(genes[genes != "" & genes != "INTERGENIC" & !is.na(genes)])
}

integrated_circRNA_host_gene_sets <- function() {
  sets <- strict_milk_host_gene_sets()
  atlas_files <- c(
    HSA = file.path(ROOT, "circAtlas_circRNA", "human_bed_v3.0.txt"),
    MMA = file.path(ROOT, "circAtlas_circRNA", "macaca_bed_v3.0.txt"),
    MMU = file.path(ROOT, "circAtlas_circRNA", "mouse_bed_v3.0.txt")
  )
  for (nm in names(atlas_files)) {
    sets[[nm]] <- unique(c(sets[[nm]], extract_circatlas_host_genes(atlas_files[[nm]])))
  }

  evidence <- read.csv(file.path(ROOT, "analysis_codex", "data", "output", "circACACA_integrated_evidence.csv"),
                       stringsAsFactors = FALSE, check.names = FALSE)
  for (i in seq_len(nrow(evidence))) {
    if (parse_bool(evidence$milk_sequencing_circACACA_detected[i])) {
      nm <- PLOT_NAMES[evidence$species[i]]
      sets[[nm]] <- unique(c(sets[[nm]], "ACACA"))
    }
  }
  sets[PLOT_NAMES[PLOT_SPECIES]]
}

save_venn <- function(sets, output_stem, main = NULL) {
  output_base <- file.path(FIG_DIR, output_stem)

  draw_one <- function() {
    par(mar = c(0, 0, 0, 0), bg = "white")
    venn(
      sets,
      snames = names(sets),
      zcolor = "style",
      ilabels = FALSE,
      ellipse = FALSE,
      opacity = 0.5,
      box = FALSE,
      borders = FALSE,
      plotsize = 15,
      ilcs = 0.8,
      sncs = 1.5,
      par = FALSE
    )
    if (!is.null(main)) {
      title(main = main, cex.main = 0.8, font.main = 2, line = -1.2)
    }
  }

  png(paste0(output_base, ".png"), width = 1800, height = 1800, res = 300)
  draw_one()
  dev.off()

  pdf(paste0(output_base, ".pdf"), width = 6, height = 6)
  draw_one()
  dev.off()

  tiff(paste0(output_base, ".tiff"), width = 1800, height = 1800, res = 300, compression = "lzw")
  draw_one()
  dev.off()

  message("Wrote: ", output_base, ".png/.pdf/.tiff")
}

combo_to_venn_code <- function(combo) {
  combo <- strsplit(combo, ";", fixed = TRUE)[[1]]
  idx <- match(combo, PLOT_SPECIES)
  idx <- sort(idx[!is.na(idx)])
  as.integer(paste(idx, collapse = ""))
}

overlay_zero_regions <- function(region_file) {
  regions <- read.csv(region_file, stringsAsFactors = FALSE, check.names = FALSE)
  zero_regions <- regions[regions$exact_region_count == 0 & regions$k >= 1, , drop = FALSE]
  if (nrow(zero_regions) == 0) return(invisible(NULL))

  icoords <- get("icoords", asNamespace("venn"))
  coords <- icoords[icoords$s == 6 & icoords$v == 0, c("x", "y", "l")]

  for (i in seq_len(nrow(zero_regions))) {
    code <- combo_to_venn_code(zero_regions$species_combo[i])
    hit <- coords[coords$l == code, , drop = FALSE]
    if (nrow(hit) == 1) {
      text(hit$x, hit$y, labels = "0", cex = 0.8, col = "black")
    }
  }
}

save_venn_with_zero_regions <- function(sets, region_file, output_stem, main = NULL) {
  output_base <- file.path(FIG_DIR, output_stem)

  draw_one <- function() {
    par(mar = c(0, 0, 0, 0), bg = "white")
    venn(
      sets,
      snames = names(sets),
      zcolor = "style",
      ilabels = FALSE,
      ellipse = FALSE,
      opacity = 0.5,
      box = FALSE,
      borders = FALSE,
      plotsize = 15,
      ilcs = 0.8,
      sncs = 1.5,
      par = FALSE
    )
    overlay_zero_regions(region_file)
    if (!is.null(main)) {
      title(main = main, cex.main = 0.8, font.main = 2, line = -1.2)
    }
  }

  png(paste0(output_base, ".png"), width = 1800, height = 1800, res = 300)
  draw_one()
  dev.off()

  pdf(paste0(output_base, ".pdf"), width = 6, height = 6)
  draw_one()
  dev.off()

  tiff(paste0(output_base, ".tiff"), width = 1800, height = 1800, res = 300, compression = "lzw")
  draw_one()
  dev.off()

  message("Wrote: ", output_base, ".png/.pdf/.tiff")
}

write_partitions <- function(sets, output_name) {
  inter <- get.venn.partitions(sets)
  for (i in seq_len(nrow(inter))) {
    inter[i, "values"] <- paste(inter[[i, "..values.."]], collapse = ",")
  }
  inter <- subset(inter, select = -..values..)
  inter <- subset(inter, select = -..set..)
  write.csv(inter, file.path(ROOT, "analysis_codex", "data", "output", output_name),
            row.names = FALSE, quote = FALSE)
}

strict_sets <- strict_milk_host_gene_sets()
integrated_sets <- integrated_circRNA_host_gene_sets()

write_partitions(strict_sets, "strict_milk_host_gene_vennPkg_partitions.csv")
write_partitions(integrated_sets, "integrated_circRNA_host_gene_vennPkg_partitions.csv")

save_venn_with_zero_regions(
  strict_sets,
  file.path(ROOT, "analysis_codex", "data", "output", "strict_milk_host_gene_6species_flower_regions.csv"),
  "FigS_6species_strict_milk_host_gene_flower_venn_vennPkg",
  main = NULL
)

save_venn_with_zero_regions(
  integrated_sets,
  file.path(ROOT, "analysis_codex", "data", "output", "integrated_circRNA_host_gene_6species_flower_regions.csv"),
  "FigS_6species_integrated_circRNA_host_gene_flower_venn_vennPkg",
  main = NULL
)
