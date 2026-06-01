# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

"The RBioc Book" — a [Quarto](https://quarto.org/) book by Sean Davis teaching introductory R, statistics, and Bioconductor, published to GitHub Pages at https://seandavi.github.io/RBiocBook/. Content is prose plus executable R code chunks; there is no application or library to build. A `DESCRIPTION` and `NAMESPACE` exist but only to declare the R package dependencies the chapters need (there is no `R/` source directory) — CI installs them via `devtools::install_deps()`.

## Building and previewing

```bash
quarto render            # render all formats (html, pdf, epub) into _book/
quarto render <file>.qmd --to html  # single chapter, fast iteration → _book/<file>.html
quarto preview           # live-reload server while editing
```

The rendering engine is `knitr`, so rendering executes the R in each chapter — a working R install with the packages from `DESCRIPTION` (`Depends:`) is required. View output by opening `_book/index.html`.

For single-chapter iteration, **always pass `--to html`**. The bare `quarto render <file>.qmd` (no `--to`) renders every configured format and writes per-chapter `*.html`, `*_files/`, and `site_libs/` artifacts to the repo root (not `_book/`) — a mess to clean up. The `--to html` form writes only `_book/<file>.html` and updates that chapter's `_freeze/<file>/execute-results/html.json`; CI re-executes the pdf/epub freeze on its own when source changes, so an html-only local render is enough to verify a prose/code edit. Bioconductor chapters self-install missing packages (e.g. `airway`) via `BiocManager` on first render, matching what `DESCRIPTION` declares.

## Architecture

- **`_quarto.yml`** is the single source of truth for book structure. The `chapters:` and `appendices:` lists control which `.qmd` files are part of the book **and their order** — a `.qmd` at the repo root that is not listed there is excluded (many such orphan/draft `.qmd` files exist). When adding a chapter, create the `.qmd` and register it here.
- **Each chapter is one `.qmd` file** at the repo root (or in a subfolder like `ggplot2/`, `machine_learning/`, `atac-seq/`, `single_cell/`). Chapters set their title via YAML frontmatter (`title: "..."`).
- **`bibliography.bib`** is the bibliography (set in `_quarto.yml`); cite with `@key`. `references.qmd` renders the reference list.
- **Three output formats** are configured: html (cosmo theme), pdf (`scrbook` via the `nmfs-opensci` titlepage extension in `_extensions/`), and epub. Changes should render cleanly in all three since CI publishes all of them.

## The `_freeze/` directory matters

`execute: freeze: auto` in `_quarto.yml` means Quarto caches each chapter's executed results under `_freeze/<chapter>/`, and **`_freeze/` is committed to git**. A chapter re-executes only when its source changes; otherwise frozen results are reused. So after editing a chapter and rendering, the updated `_freeze/<chapter>/` files (including `figure-*/` images) are expected changes — commit them alongside the `.qmd` edit. This keeps CI fast and deterministic.

Note: per-chapter `*_cache/` and `*_files/` directories and `figure/` are gitignored — only `_freeze/` carries committed results.

## Publishing (CI)

`.github/workflows/R-CMD-check.yaml` runs on push/PR to `main`: on macOS it sets up R, installs deps via `devtools::install_deps(dependencies=TRUE)`, sets up Quarto with TinyTeX (for PDF), installs the `nmfs-opensci/quarto_titlepages` extension, then renders and publishes (to `gh-pages`). There is no `renv` despite the package-style scaffold — dependencies come from `DESCRIPTION`.

## Conventions

- `_book/` and the `index.*` LaTeX intermediates at the repo root (`.tex`, `.aux`, `.log`, `.pdf`, `.toc`, etc.) are build artifacts — don't hand-edit them.
- This is a teaching book: keep prose accessible to R beginners and keep code chunks runnable and self-contained within a chapter.
