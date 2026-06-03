# Shared setup for scripts in New-analysis.
# Run scripts from the project root or from New-analysis; set CIRCRNA_PROJECT_ROOT
# if you need to launch them from another working directory.

find_project_root <- function() {
  candidates <- c(
    Sys.getenv("CIRCRNA_PROJECT_ROOT", unset = NA_character_),
    getwd(),
    dirname(getwd())
  )
  candidates <- unique(candidates[!is.na(candidates) & nzchar(candidates)])

  for (candidate in candidates) {
    candidate <- normalizePath(candidate, winslash = "/", mustWork = FALSE)
    if (
      dir.exists(file.path(candidate, "New-analysis")) &&
        dir.exists(file.path(candidate, "genome and annotation"))
    ) {
      return(candidate)
    }
  }

  stop(
    "Cannot find project root. Run from the project root/New-analysis, ",
    "or set CIRCRNA_PROJECT_ROOT.",
    call. = FALSE
  )
}

PROJECT_ROOT <- find_project_root()
setwd(PROJECT_ROOT)
