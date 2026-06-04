# Chapter harmonization ranking

Worst-first priority list for the pedagogical harmonization effort. Scores are the
**mean of the 7-dimension 1–5 rubric** in [`STYLE_GUIDE.md`](STYLE_GUIDE.md) §4
(5 = matches the class-tested intro chapters). Full per-dimension scorecards for
the three pilot chapters are in [`PILOT_REVIEWS.md`](PILOT_REVIEWS.md).

*Scored 2026-06-01 by a panel of parallel reviewers, calibrated to the three
pilots (`dataframes_intro` 3.3, `t-stats-and-tests` 2.4, `bioc-summarizedexperiment`
2.0). A snapshot — re-score if chapters change substantially.*

## Done (merged to main)

| Chapter | Was | Now (self-est.) |
|---|:--:|:--:|
| `bioc-summarizedexperiment` | 2.0 | harmonized |
| `t-stats-and-tests` | 2.4 | harmonized |
| `norm` | (cold) | harmonized |
| `factors` | 1.7 | ~4.4 |
| `vectors` | 2.4 | ~4.5 |
| `matrices` | 2.4 | ~4.4 |
| `lists` | 3.1 | ~4.4 |
| `machine_learning/models` | 1.4 | ~4.7 |
| `visualization_guide` | 1.6 | ~4.7 |
| `ranges_and_signals` | 1.6 | ~4.3 |
| `intro` | 1.9 | ~4.7 |
| `kmeans` | 1.9 | ~4.5 |
| `eda_and_univariate_brfss` | 2.0 | ~4.5 |
| `atac-seq/atac-seq` | 1.6 | ~4.5 |
| `geoquery` | 1.7 | ~4.4 |
| `genomic_ranges_tutorial` | 2.1 | ~4.6 |
| `machine_learning/intro` | 2.3 | ~4.5 |
| `machine_learning/mlr3verse_intro` | 2.1 | ~4.9 |
| `machine_learning/machine_learning_mlr3` | 2.9 | ~4.5 (split — see below) |
| `dataframes_intro` | 3.3 | ~4.7 |
| `reading_and_writing` | 2.3 | ~4.7 |
| `dplyr_intro_msleep` | 2.4 | ~4.7 |
| `ggplot2/intro` | 2.7 | ~4.7 |
| `310_microbiome` | 2.3 | ~4.7 |
| `single_cell/setup` | 2.6 | ~4.3 |
| `ranges_exercises` | 2.4 | ~4.9 |

## Remaining — worst first

**None — the harmonization sweep is complete.** All instructional chapters tracked
here have been raised to the house style (motivating open, single objectives
block, callout taxonomy, exercises with collapsible solutions, summary, `<-`
throughout). Re-score if chapters change substantially, or fold in any chapters
not previously tracked.

## Natural batches

- **ML part — DONE.** All four chapters harmonized (`intro`, `models`,
  `mlr3verse_intro`, plus the worked-examples capstone). The capstone
  `machine_learning_mlr3.qmd` (992 lines after harmonization) was **split into
  three self-contained worked-example chapters**: `ml_cancer.qmd` (cancer
  classification from gene expression), `ml_methylation_age.qmd` (age from DNA
  methylation), and `ml_expression.qmd` (expression from histone marks).
- **Ranges / Bioc — DONE.** `ranges_exercises`, `ranges_and_signals`,
  `genomic_ranges_tutorial`, `atac-seq`, and `geoquery` are all harmonized.
- **Final sweep — DONE.** `dataframes_intro`, `reading_and_writing`,
  `dplyr_intro_msleep`, `ggplot2/intro`, `310_microbiome`, `single_cell/setup`,
  and `ranges_exercises` were harmonized together (one author per chapter,
  parallel across chapters) and each rendered cleanly to HTML with its
  `_freeze/` regenerated.

## Process reminders (from `STYLE_GUIDE.md` / `PILOT_REVIEWS.md`)

- One author per chapter (don't split a chapter's prose — it fractures voice);
  parallelize *across* chapters. Anchor voice to `r_basics.qmd`, structure to an
  already-harmonized chapter (e.g. `norm.qmd`).
- Before editing a chapter, run `Rscript tools/check-deps.R` and add any missing
  package to `DESCRIPTION` (editing a chapter forces CI re-execution that needs
  its deps).
- Audience is **not** very stats-savvy — bias toward accessible prose; prefer
  inline `r ...` over hardcoded numbers; seed random chunks for reproducible
  freeze.
