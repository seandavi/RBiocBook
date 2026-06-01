#!/usr/bin/env Rscript
# Dependency-completeness lint for "The RBioc Book".
#
# Scans every chapter listed in _quarto.yml for the R packages it uses
# (library/require/requireNamespace/pkg::) and fails if any are NOT declared in
# DESCRIPTION's `Depends:` field. Base and recommended packages (always shipped
# with R) are exempt.
#
# Why: the book renders in CI by reusing committed `_freeze/`, so a chapter that
# library()s an undeclared package still "passes" — until the chapter is edited
# and re-executes, at which point CI dies with "there is no package called X".
# Running this BEFORE the render step turns that 10-minute render failure into a
# 5-second lint failure with a precise message.
#
# Usage (from the repo root):  Rscript tools/check-deps.R

quarto_yml  <- "_quarto.yml"
description <- "DESCRIPTION"
stopifnot(file.exists(quarto_yml), file.exists(description))

## 1. In-book chapter files from _quarto.yml (commented-out entries excluded) --
yml <- readLines(quarto_yml, warn = FALSE)
yml <- yml[!grepl("^\\s*#", yml)]                      # drop full-line comments
qmd <- unlist(regmatches(yml, gregexpr("[A-Za-z0-9_./-]+\\.qmd", yml)))
chapters <- unique(qmd)
chapters <- chapters[file.exists(chapters)]
if (!length(chapters)) stop("No chapter .qmd files resolved from _quarto.yml")

## 2. Packages referenced in the R code chunks of one chapter -----------------
pkgs_in_file <- function(path) {
  lines <- readLines(path, warn = FALSE)
  in_chunk <- FALSE
  skip    <- FALSE          # current chunk is eval=FALSE -> ignore its code
  code <- character()
  for (ln in lines) {
    if (!in_chunk && grepl("^```+\\{[rR][ ,}]", ln)) {
      in_chunk <- TRUE
      skip <- grepl("eval\\s*=\\s*FALSE", ln)          # old-style header option
      next
    }
    if (in_chunk && grepl("^```+\\s*$", ln)) { in_chunk <- FALSE; skip <- FALSE; next }
    if (in_chunk) {
      if (grepl("^\\s*#\\|\\s*eval:\\s*false", ln, ignore.case = TRUE)) { skip <- TRUE; next }
      if (skip) next                                    # non-executing chunk
      if (grepl("^\\s*#", ln)) next                     # comment / chunk-option line
      code <- c(code, ln)
    }
  }
  txt  <- paste(code, collapse = "\n")
  grab <- function(re) {
    m <- regmatches(txt, gregexpr(re, txt, perl = TRUE))[[1]]
    if (!length(m)) return(character())
    unique(sub(re, "\\1", m, perl = TRUE))
  }
  unique(c(
    grab("library\\(\\s*[\"']?([A-Za-z][A-Za-z0-9._]*)"),
    grab("require\\(\\s*[\"']?([A-Za-z][A-Za-z0-9._]*)"),
    grab("requireNamespace\\(\\s*[\"']([A-Za-z][A-Za-z0-9._]*)"),
    grab("\\b([A-Za-z][A-Za-z0-9._]*):::?")
  ))
}

ref <- lapply(chapters, function(f)
  tryCatch(pkgs_in_file(f), error = function(e) character()))
names(ref) <- chapters
all_ref <- sort(unique(unlist(ref)))

## 3. Declared deps + base/recommended packages ------------------------------
dep_field <- read.dcf(description)[1, "Depends"]
declared  <- trimws(gsub("[(].*?[)]", "", strsplit(dep_field, ",")[[1]]))
declared  <- declared[declared != "" & declared != "R"]
standard  <- rownames(installed.packages(priority = c("base", "recommended")))

## 4. Gap + report -----------------------------------------------------------
missing <- sort(setdiff(all_ref, c(declared, standard)))

if (!length(missing)) {
  cat(sprintf("OK  dep-completeness: all %d packages used across %d chapters are declared.\n",
              length(all_ref), length(chapters)))
  quit(status = 0)
}

cat("FAIL  dep-completeness: packages used by a chapter but missing from\n")
cat("      DESCRIPTION Depends: (and not base/recommended):\n\n")
for (p in missing) {
  where <- names(ref)[vapply(ref, function(v) p %in% v, logical(1))]
  cat(sprintf("  - %-24s used in: %s\n", p, paste(where, collapse = ", ")))
}
cat("\nFix: add them to DESCRIPTION Depends: (or drop the usage).\n")
quit(status = 1)
