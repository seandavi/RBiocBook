# Chapter Review Findings

**Reviewed:** 2026-05-31 · **Chapters reviewed:** 40 (all `chapters:` + `appendices:` in `_quarto.yml`) · **Findings:** 82 across 28 chapters · **Clean:** 12 chapters.

Every finding below was produced by a per-chapter reviewer **and then grep-verified against the source file** — each one quotes text that literally exists at the cited line. No changes were made. Line numbers are as of this date.

Severity guide:
- **high** — would error at render, breaks meaning, or is a clearly wrong fact/value.
- **medium** — real defect (typo that garbles a word, malformed markup) but renders.
- **low** — minor typo or style nit.

## Summary

| Chapter | High | Med | Low | Total |
|---|---:|---:|---:|---:|
| index.qmd | 1 | 1 | 0 | 2 |
| intro.qmd | 2 | 1 | 0 | 3 |
| intro_to_rstudio.qmd | 0 | 1 | 0 | 1 |
| r_intro_mechanics.qmd | 0 | 1 | 0 | 1 |
| r_basics.qmd | 0 | 1 | 0 | 1 |
| packages_and_dice.qmd | 0 | 0 | 3 | 3 |
| reading_and_writing.qmd | 1 | 0 | 0 | 1 |
| vectors.qmd | 1 | 1 | 0 | 2 |
| dataframes_intro.qmd | 0 | 0 | 3 | 3 |
| factors.qmd | 1 | 0 | 0 | 1 |
| dplyr_intro_msleep.qmd | 0 | 0 | 3 | 3 |
| eda_and_univariate_brfss.qmd | 1 | 1 | 0 | 2 |
| ggplot2/intro.qmd | 4 | 0 | 1 | 5 |
| visualization_guide.qmd | 2 | 1 | 1 | 4 |
| norm.qmd | 2 | 2 | 0 | 4 |
| t-stats-and-tests.qmd | 7 | 0 | 0 | 7 |
| kmeans.qmd | 2 | 0 | 0 | 2 |
| machine_learning/intro.qmd | 1 | 0 | 0 | 1 |
| machine_learning/machine_learning_mlr3.qmd | 5 | 1 | 0 | 6 |
| geoquery.qmd | 1 | 0 | 0 | 1 |
| bioc-summarizedexperiment.qmd | 2 | 2 | 0 | 4 |
| 310_microbiome.qmd | 0 | 2 | 0 | 2 |
| genomic_ranges_tutorial.qmd | 2 | 1 | 0 | 3 |
| ranges_and_signals.qmd | 2 | 3 | 1 | 6 |
| ranges_exercises.qmd | 2 | 1 | 1 | 4 |
| atac-seq/atac-seq.qmd | 2 | 1 | 0 | 3 |
| single_cell/setup.qmd | 2 | 3 | 0 | 5 |
| appendix.qmd | 1 | 1 | 0 | 2 |

## Highest-impact items (would break rendering or teach something wrong)

- **Code that errors at render:** `eda_and_univariate_brfss.qmd:64` (`path` undefined — only set in an `eval=FALSE`/commented chunk), `ggplot2/intro.qmd` `label =` should be `labels =` (4×), `kmeans.qmd` chunk labels written as `{r #fig-...}`, `geoquery.qmd:92` stray space in `pca_results $sdev`, `visualization_guide.qmd:114` non-existent `colorblindFriendly=` arg, `machine_learning_mlr3.qmd` (`title=` on `boxplot`, `pred_train` used before defined), `single_cell/setup.qmd:217` trailing-comma syntax error + `:171` malformed callout, `factors.qmd`/`appendix.qmd` R Markdown `output:` YAML in a Quarto book.
- **Wrong facts/values:** `genomic_ranges_tutorial.qmd:111` ("1-based" contradicts the chapter's own 0-based explanation), `ranges_and_signals.qmd:568` (says hg19 for an hg38 package), `bioc-summarizedexperiment.qmd:208` (range `100000:1100000` ≠ the "100,000 to 110,000" in the comment), `ranges_exercises.qmd:114` (computes on `dnase`, question asks `dnase2`).
- **Broken cross-ref / markup:** `norm.qmd:36` (`@fig-pnorm-1` — label is `fig-pnorm`), `atac-seq/atac-seq.qmd:204` (`\@ref{...}` curly braces), `t-stats-and-tests.qmd:116,131` (mismatched LaTeX braces).

---

## Findings by chapter

### index.qmd
*Two issues: a malformed sentence and a duplicated word.*

| Sev | Cat | Line | Problem |
|---|---|---|---|
| high | typo | 141 | Duplicated word: "we will will not be prescriptive". |
| medium | grammar | 76 | Run-on/dangling clause: "…approach that focuses on the individual needs, preferences, and interests of each participant can greatly enhance…" — relative clause never resolves to a main predicate. |

### intro.qmd
*Two typos and one stray markup fragment.*

| Sev | Cat | Line | Problem |
|---|---|---|---|
| high | typo | 50 | "ara" should be "are" ("…extensively tested and ara available"). |
| high | typo | 66 | "interfact" should be "interface" ("R can interfact with FORTRAN"). |
| medium | formatting | 86 | Stray `]{}` at end of line ("…packages are unreliable.]{}") — leftover span syntax. |

### intro_to_rstudio.qmd
| Sev | Cat | Line | Problem |
|---|---|---|---|
| medium | typo | 24 | Missing space: "(no parentheses)in the R console" → "parentheses) in". |

### r_intro_mechanics.qmd
| Sev | Cat | Line | Problem |
|---|---|---|---|
| medium | typo | 56 | "fillowing" should be "following". |

### r_basics.qmd
| Sev | Cat | Line | Problem |
|---|---|---|---|
| medium | typo | 398 | "favlue" should be "value" ("does not change the favlue"). |

### packages_and_dice.qmd
*Minor style/formatting nits only.*

| Sev | Cat | Line | Problem |
|---|---|---|---|
| low | grammar | 12 | Semicolon before a list ("using two powerful and general tools;") — should be a colon. |
| low | formatting | 48 | `install.packages("dplyr")` in prose not wrapped in backticks. |
| low | grammar | 94 | "eg.," should be "e.g.,". |

### reading_and_writing.qmd
| Sev | Cat | Line | Problem |
|---|---|---|---|
| high | typo | 182 | "CSV riles" should be "CSV files". |

### vectors.qmd
| Sev | Cat | Line | Problem |
|---|---|---|---|
| high | broken-code | 331 | Exercise answer has 7 values `c(72, 75, 78, 81, 76, 73, 93)` but the prompt asks for the 6 listed values — extra `93`. |
| medium | formatting | 334 | Exercise list re-starts at "1." (should be "2."). |

### dataframes_intro.qmd
*Three typos; code is fine.*

| Sev | Cat | Line | Problem |
|---|---|---|---|
| low | typo | 18 | "examing" should be "examining". |
| low | typo | 193 | "look a the" should be "look at the". |
| low | typo | 257 | "varibles" should be "variables". |

### factors.qmd
| Sev | Cat | Line | Problem |
|---|---|---|---|
| high | broken-code | 3 | R Markdown frontmatter `output: html_document` in a Quarto book chapter — format is set globally in `_quarto.yml`; other chapters omit this. |

### dplyr_intro_msleep.qmd
| Sev | Cat | Line | Problem |
|---|---|---|---|
| low | typo | 28 | Heading "Why Is dplyr userful?" → "useful". |
| low | formatting | 165 | Heading has a doubled quote: `## "Piping"" with` `\|>`. |
| low | typo | 180 | Backslash-escaped underscore in prose: `sleep\_total` (escaping only needed in tables). Also at 201, 219, 253. |

### eda_and_univariate_brfss.qmd
| Sev | Cat | Line | Problem |
|---|---|---|---|
| high | broken-code | 64 | `stopifnot(file.exists(path))` runs in an executed chunk, but `path` is only set in an `eval=FALSE` chunk (60) and a commented-out chunk (55–57) — errors on render. |
| medium | formatting | 198 | Exercise numbered "2." duplicates the "2." at line 187 (should be "3."). |

### ggplot2/intro.qmd
| Sev | Cat | Line | Problem |
|---|---|---|---|
| high | broken-code | 145 | `scale_y_continuous(label = scales::dollar)` — arg is `labels=`, not `label=`. |
| high | broken-code | 170 | Same `label=` → `labels=`. |
| high | broken-code | 197 | Same `label=` → `labels=`. |
| high | broken-code | 229 | Same `label=` → `labels=`. |
| low | typo | 160 | Comment: "obsese" should be "obese". |

*(Note: ggplot tolerates `label` as a partial-match for `labels` in some versions, so these may render today — but they are not the documented argument and are fragile.)*

### visualization_guide.qmd
| Sev | Cat | Line | Problem |
|---|---|---|---|
| high | broken-code | 114 | `display.brewer.all(colorblindFriendly=TRUE)` — RColorBrewer's function has no `colorblindFriendly` argument. |
| high | grammar | 96 | "They not imply" — missing verb ("They do not imply"). |
| medium | formatting | 86 | `$e.g.$` rendered in math mode; should be plain "e.g.". |
| low | grammar | 7 | "there are a some truly amazing" — "a some". |

### norm.qmd
| Sev | Cat | Line | Problem |
|---|---|---|---|
| high | broken-link | 36 | `@fig-pnorm-1` cross-ref — the figure label is `fig-pnorm` (defined ~line 149). |
| high | missing-content | 140 | Lines ~140–143 duplicate the `pnorm` template/explanation from lines ~35–43 (the template paragraph appears twice). |
| medium | formatting | 237 | Missing space: "the inverse of the`pnorm` function" → "the `pnorm`". |
| medium | missing-content | 158 | `lower.tail` paragraph duplicates the explanation at lines ~145–146. |

### t-stats-and-tests.qmd
| Sev | Cat | Line | Problem |
|---|---|---|---|
| high | typo | 11 | "based no the data" → "based on". |
| high | typo | 25 | "distriubution" → "distribution". |
| high | typo | 37 | "is a the population mean" → "is the". |
| high | broken-code | 116 | LaTeX `\sigma = \sqrt{\frac{1}{N}\sum...({xi - \mu)^2}}` — mismatched braces around the squared term. |
| high | broken-code | 131 | LaTeX `s = \sqrt{...({x_{i} - \bar{x})^2}}` — same mismatched-brace pattern. |
| high | factual-error | 288 | `hist(pt(t10k,5))` — samples are size 5, so df should be 4 (n−1), not 5. |
| high | typo | 459 | "is now easily calculated" → "is not easily calculated" (reverses the meaning). |

### kmeans.qmd
| Sev | Cat | Line | Problem |
|---|---|---|---|
| high | broken-code | 186 | Chunk header `{r #fig-sds-derisi, ...}` — knitr label must not have `#`; should be `{r fig-sds-derisi, ...}`. |
| high | broken-code | 216 | Chunk header `{r #fig-kmeans-derisi, ...}` — same `#`-prefix problem. |

### machine_learning/intro.qmd
| Sev | Cat | Line | Problem |
|---|---|---|---|
| high | typo | 54 | "neaural networks" → "neural". |

### machine_learning/machine_learning_mlr3.qmd
| Sev | Cat | Line | Problem |
|---|---|---|---|
| high | broken-code | 313 | Sentence cuts off mid-word: "to get a vector of stan". |
| high | grammar | 440 | "Lets use our trained model works to predict" — garbled. |
| high | broken-code | 898 | `boxplot(fullFeatureSet, title='Original data')` — `boxplot` has no `title=`; use `main=`. |
| high | broken-code | 1098 | `pred_train$score(measures)` uses `pred_train` before it is created (line ~1115). |
| high | broken-code | 1093 | `learner$train(exp_pred_task)` lacks the split/partition used by the other training calls in the chapter. |
| medium | typo | 604 | "ont he ranger package" → "on the". |

### geoquery.qmd
| Sev | Cat | Line | Problem |
|---|---|---|---|
| high | broken-code | 92 | `sum(pca_results $sdev^2)` — stray space breaks `$` access; should be `pca_results$sdev^2`. |

### bioc-summarizedexperiment.qmd
| Sev | Cat | Line | Problem |
|---|---|---|---|
| high | factual-error | 208 | `GRanges(seqnames="1", ranges=100000:1100000)` spans 1,000,000 bp, but the comment (line 206) says "100,000 to 110,000" — likely meant `100000:110000`. |
| high | grammar | 68 | Broken fragment: "experimental and assays 64,102 gene transcripts." |
| medium | typo | 47 | Doubled quote: `the "child"" of the`. |
| medium | typo | 262 | Table cell "Green luorescence" → "fluorescence". |

### 310_microbiome.qmd
| Sev | Cat | Line | Problem |
|---|---|---|---|
| medium | typo | 135 | "rows or columns it you like" → "if you like". |
| medium | typo | 303 | Chunk option `{r messasge=FALSE}` → `message=FALSE` (so the option is silently ignored). |

### genomic_ranges_tutorial.qmd
| Sev | Cat | Line | Problem |
|---|---|---|---|
| high | broken-code | 78–79 | `readr::read_table(...)` but `readr` is never loaded (setup loads rtracklayer, GenomicRanges, ggplot2, dplyr). |
| high | factual-error | 111 | Comment "End positions (1-based)" contradicts the chapter's own 0-based BED explanation (lines 29–30, 46, 54). |
| medium | missing-content | 298 | Key Takeaways lists "Finding overlaps and performing grouped operations" as covered, but the tutorial never demonstrates it. |

### ranges_and_signals.qmd
| Sev | Cat | Line | Problem |
|---|---|---|---|
| high | factual-error | 568 | Text says "hg19 build" for what the surrounding code uses as the hg38 knownGene package. |
| high | broken-code | 499 | `findOverlaps(gr1, gr2, select="first")` appears twice in a row (lines 499 & 502). |
| medium | typo | 354 | "are are comprised of exons" — doubled "are". |
| medium | typo | 356 | "we thing of each transcript" → "think". |
| medium | grammar | 359 | Fragment: "…list of `GRanges` objects that. Continuing with…". |
| low | typo | 61 | Table: "single nicleotide locations" → "nucleotide". |

### ranges_exercises.qmd
| Sev | Cat | Line | Problem |
|---|---|---|---|
| high | broken-code | 114 | Question asks about `dnase2` but code computes `sum(width(dnase))/sum(seqlengths(dnase))`. |
| high | broken-code | 202 | `prop_dnase = sum(width(dnase))/sum(seqlengths(prom_regions))` — denominator uses the promoter regions' seqlengths, inconsistent normalization. |
| medium | typo | 98 | "Are there are similar number?" — doubled "are". |
| low | typo | 203 | "# Iff the dnase…" — "Iff" (math: iff) where plain "If" is meant. |

### atac-seq/atac-seq.qmd
| Sev | Cat | Line | Problem |
|---|---|---|---|
| high | broken-code | 204 | Cross-ref `\@ref{fig:flgreenleaf}` uses braces; Quarto/bookdown needs `\@ref(fig:flgreenleaf)`. |
| high | broken-code | 21 | `Biocpkg('heatmaps')` — likely wrong package name (verify against what the code chunks actually load; common one is ComplexHeatmap). |
| medium | typo | 74 | Table: "genomic reagions" → "regions". |

### single_cell/setup.qmd
| Sev | Cat | Line | Problem |
|---|---|---|---|
| high | formatting | 171 | Malformed callout: `::: {.callout-tip)` — `)` should be `}`. |
| high | missing-content | 276 | Dangling fragment "See the" with no continuation. |
| medium | broken-code | 217 | `plotReducedDim(sce, "TSNE", ..., shape_by="tissue", )` — trailing comma. |
| medium | grammar | 149 | "PCA is a widely _linear_ dimensionality reduction technique" — garbled (missing "used"). |
| medium | broken-code | 242 | Comment ends with an unmatched "(": "…in a violin plot(". |

### appendix.qmd
| Sev | Cat | Line | Problem |
|---|---|---|---|
| high | broken-code | 4 | R Markdown `output: { html_document, pdf_document }` YAML in a Quarto book chapter (format is controlled by `_quarto.yml`). |
| medium | missing-content | 22 | swirl blockquote ends abruptly ("…have no fear.") — looks truncated. |

---

## Clean chapters (no obvious problems found)

`data_structures_overview.qmd` · `matrices.qmd` · `lists.qmd` · `eda_overview.qmd` · `machine_learning/models.qmd` · `machine_learning/mlr3verse_intro.qmd` · `references.qmd` · `git_and_github.qmd` · `additional_resources.qmd` · `ai_tools.qmd` · `dataviz.qmd` · `matrix_exercises.qmd`

---

*Scope: obvious, verifiable defects only (broken code, wrong facts, malformed markup, clear typos). Stylistic/subjective issues were intentionally excluded. This was a static read — chapters were not rendered, so runtime errors that only surface with live data or specific package versions may not all be captured here. The most reliable next step for the code findings is `quarto render <chapter>.qmd` on the flagged chapters.*
